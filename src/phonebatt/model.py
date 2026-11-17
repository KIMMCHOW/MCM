"""Fast L2/L3 system x = (z, v1, T) and the time-to-empty solvers (FINAL_SPEC §2.13, §6).

* ``fast_state`` / ``rhs_batch`` / ``event_batch`` are vectorised over a trailing sample axis: ``x`` has shape
  ``(3, n)`` and the input arrays and any ``Params`` field broadcast to ``(n,)``.
* ``Segment`` wraps one constant input for the reference LSODA solver (§6.2).
* ``run`` integrates one trajectory (LSODA reference, or RK4); ``run_batch`` integrates many parameter samples
  with ``integrate.rk4_batch`` (§6.3) under (a) a cyclic deterministic schedule or (b) CTMC usage paths.

Events (§6.1): charge (z), cutoff (V − V_cut), brownout (V_pk − V_sys, only while ΔP_b > 0 and the envelope
is enabled), infeasible (D ≤ 0) and thermal (T_shut − T).  The batch solver uses the combined event function
g = min(z/0.01, (V − V_cut)/0.1, (V_pk − V_sys)/0.1, (T_shut − T)/10, D/Ũ²), NaN when D ≤ 0, and identifies
the event type at the crossing with :func:`classify` ('infeasible' whenever D ≤ 0; simultaneous violations at a
switch resolved in the order charge, cutoff, brownout, thermal; otherwise the argmin).
"""
from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, fields, replace
from typing import NamedTuple

import numpy as np

from . import cell
from .cell import DUDT_V, DUDT_Z, current, peak_voltage
from .integrate import rk4_batch, solve_ivp_segmented
from .params import K0, RG, RP63, T_REF, Params, with_theta
from .power import power_parts, rho_perf, throttle
from .usage import CTMC, PRIMITIVES, Paths, Schedule, UsageState

EVENTS = ("charge", "cutoff", "brownout", "thermal", "infeasible")


class Inputs(NamedTuple):
    """Per-sample segment inputs: throttling-independent power, SoC demand, burst amplitude [W], T_amb [K]."""

    Pns: np.ndarray
    dem: np.ndarray
    dPb: np.ndarray
    Tamb: np.ndarray


# ------------------------------------------------------------------------------------------- vectorised core
def fast_state(x, inp: Inputs, p: Params) -> dict:
    """Algebraic quantities of the fast system at state x (§2.2, §2.5–2.9).  Works on (3,) or (3, n).

    A fused evaluation of ``cell.ocv``, ``cell.R0``, ``cell.R1``, ``cell.tau1``, ``power.throttle``,
    ``cell.current`` and ``cell.peak_voltage`` that shares the Arrhenius and entropic sub-expressions
    (tests/test_core.py checks it against the separate functions).
    """
    z, v1, T = x[0], x[1], x[2]
    zc = np.minimum(np.maximum(z, 0.0), 1.0)
    if p.dudt_table is None:
        dudt = np.interp(zc, DUDT_Z, DUDT_V)
    else:
        dudt = np.interp(zc, *p.dudt_table)
    U = cell.ocv_ref(zc, p, clamped=True) + (T - T_REF) * p.k_S * dudt
    ar = (1.0 / T - 1.0 / T_REF) / RG
    g = 1.0 + p.rg
    rinf = p.R0_ref / (1.0 + p.A_R * math.exp(-0.5 / p.z0)) if np.ndim(p.A_R) == 0 and np.ndim(p.z0) == 0 \
        else p.R0_ref / (1.0 + p.A_R * np.exp(-0.5 / p.z0))
    R0_ = rinf * (1.0 + p.A_R * np.exp(-np.maximum(z, 0.0) / p.z0)) * np.exp(p.Ea_R0 * ar) * g
    R1_ = p.R1_ref * np.exp(p.Ea_R1 * ar) * g
    tau = p.tau1_ref * np.exp(p.Ea_tau * ar)
    Ts = inp.Tamb + p.R_skin / p.R_th * (T - inp.Tamb)
    Psoc = throttle(inp.dem, Ts, p)
    P = inp.Pns + Psoc
    Ue = U - v1
    I, V, D = current(Ue, R0_, P)
    Vpk = peak_voltage(Ue, R0_, P + inp.dPb)
    return dict(U=U, R0=R0_, R1=R1_, tau1=tau, dUdT=dudt, Ts=Ts, Psoc=Psoc, P=P, Ue=Ue, I=I, V=V, D=D, Vpk=Vpk)


def rhs_from_state(x, st: dict, inp: Inputs, p: Params):
    """Right-hand side of §2.13 given the algebraic state ``st`` (E3, E4, E7)."""
    z, v1, T = x[0], x[1], x[2]
    I, R0_, R1_ = st["I"], st["R0"], st["R1"]
    q = (I * st["V"] if p.load_heat else 0.0) + I * I * R0_ + v1 * v1 / R1_ - I * T * p.k_S * st["dUdT"]
    dz = -I / (3600.0 * p.Q_max)
    dv1 = (R1_ * I - v1) / st["tau1"]
    dT = (q - (T - inp.Tamb) / (p.k_env * p.R_th)) / p.C_th
    return dz, dv1, dT


def _scalar_state(z, v1, T, inp: Inputs, p: Params, c: dict) -> dict:
    """Pure-Python scalar version of :func:`fast_state` (LSODA speed); ``c`` holds per-Params constants."""
    zc = min(max(z, 0.0), 1.0)
    i = min(int(zc / 0.05), 19)
    dudt = DUDT_V[i] + (DUDT_V[i + 1] - DUDT_V[i]) * (zc - DUDT_Z[i]) / 0.05
    U = (-c["b0"] * math.exp(-c["b1"] * zc) + c["a0"] + zc * (c["a1"] + zc * (c["a2"] + zc * c["a3"]))
         + c["dV_hv"] * zc ** c["p_hv"] + c["dU"] + (T - T_REF) * c["k_S"] * dudt)
    ar = (1.0 / T - 1.0 / T_REF) / RG
    R0_ = c["rinf"] * (1.0 + c["A_R"] * math.exp(-max(z, 0.0) / c["z0"])) * math.exp(c["Ea_R0"] * ar) * c["g"]
    R1_ = c["R1_ref"] * math.exp(c["Ea_R1"] * ar) * c["g"]
    tau = c["tau1_ref"] * math.exp(c["Ea_tau"] * ar)
    Tamb = inp.Tamb
    Ts = Tamb + c["skin"] * (T - Tamb)
    dem = inp.dem
    if c["throttling"] and dem > c["P_fl"]:
        y = (c["T_thr"] - Ts) / c["w_thr"]
        Psoc = c["P_fl"] + (dem - c["P_fl"]) / (1.0 + math.exp(-y)) if y > -700 else c["P_fl"]
    else:
        Psoc = dem
    P = inp.Pns + Psoc
    Ue = U - v1
    D = Ue * Ue - 4.0 * R0_ * P
    if D > 0:
        sq = math.sqrt(D)
        I = 2.0 * P / (Ue + sq)
    else:
        I = Ue / (2.0 * R0_)
    V = Ue - I * R0_
    Dk = Ue * Ue - 4.0 * R0_ * (P + inp.dPb)
    Vpk = 0.5 * (Ue + math.sqrt(Dk)) if Dk > 0 else 0.0
    return dict(U=U, R0=R0_, R1=R1_, tau1=tau, dUdT=dudt, Ts=Ts, Psoc=Psoc, P=P, Ue=Ue, I=I, V=V, D=D, Vpk=Vpk)


def _scalar_constants(p: Params) -> dict | None:
    """Constants for :func:`_scalar_state`, or None if p has array fields or tables (use the vector path)."""
    if p.n_samples() > 1 or p.ocv_table is not None or p.dudt_table is not None:
        return None
    names = ("b0", "b1", "a0", "a1", "a2", "a3", "dV_hv", "p_hv", "dU", "k_S", "A_R", "z0", "Ea_R0", "R1_ref",
             "tau1_ref", "Ea_R1", "Ea_tau", "P_fl", "w_thr")
    c = {k: float(getattr(p, k)) for k in names}
    c["rinf"] = float(p.R0_ref) / (1.0 + c["A_R"] * math.exp(-0.5 / c["z0"]))
    c["g"] = 1.0 + float(p.rg)
    c["skin"] = float(p.R_skin) / float(p.R_th)
    c["T_thr"] = float(p.T_lim) - float(p.dT_thr)
    c["throttling"] = bool(p.throttling)
    return c


def rhs_batch(x, inp: Inputs, p: Params, st: dict | None = None) -> np.ndarray:
    """dx/dt for a batch, shape (3, n)."""
    st = fast_state(x, inp, p) if st is None else st
    out = np.empty(np.shape(x))
    out[0], out[1], out[2] = rhs_from_state(x, st, inp, p)
    return out


def event_components(x, inp: Inputs, p: Params, st: dict | None = None) -> np.ndarray:
    """Normalised event functions [z/0.01, (V−V_cut)/0.1, (V_pk−V_sys)/0.1, (T_shut−T)/10, D/Ũ²], (5, n)."""
    st = fast_state(x, inp, p) if st is None else st
    z, T = x[0], x[2]
    active = (inp.dPb > 0) & bool(p.burst)
    g3 = np.where(active, (st["Vpk"] - p.V_sys) / 0.1, np.inf)
    comps = np.broadcast_arrays(z / 0.01, (st["V"] - p.V_cut) / 0.1, g3, (p.T_shut - T) / 10.0,
                                st["D"] / (st["Ue"] * st["Ue"]))
    return np.stack(comps)


def event_batch(x, inp: Inputs, p: Params, st: dict | None = None) -> np.ndarray:
    """Combined event function of §6.3 (NaN where D ≤ 0)."""
    c = event_components(x, inp, p, st)
    g = np.min(c, axis=0)
    return np.where(c[4] > 0, g, np.nan)


def classify(comps: np.ndarray) -> np.ndarray:
    """Event name per sample from the event components (5, n) evaluated at the crossing.

    * 'infeasible' if D ≤ 0 (no operating point; the clamped voltage is meaningless);
    * else the first violated (≤ 0) function in the order charge, cutoff, brownout, thermal — this decides
      simultaneous violations at a switching instant (§6.1) and matches the prototype's cut-off-first rule;
    * else (an interpolated smooth crossing, all components ≈ 0⁺) the argmin (§6.3).
    """
    names = np.array(EVENTS[:4], dtype=object)
    g = comps[:4]
    viol = g <= 0
    first = np.argmax(viol, axis=0)
    out = np.where(viol.any(axis=0), names[first], names[np.argmin(g, axis=0)]).astype(object)
    out[~(comps[4] > 0)] = "infeasible"
    return out


# ------------------------------------------------------------------------------------------- scalar segment
class Segment:
    """One constant input u at ambient T_amb [K] for the LSODA reference solver (scalar state)."""

    def __init__(self, u: UsageState, Tamb_K: float, p: Params):
        self.u, self.Tamb, self.p = u, float(Tamb_K), p
        Pns, dem, dPb, br = power_parts(u, p)
        self.breakdown = {k: float(v) for k, v in br.items()}
        self.inp = Inputs(float(Pns), float(dem), float(dPb), self.Tamb)
        self.brownout_active = bool(p.burst) and float(dPb) > 0
        self._c = _scalar_constants(p)
        self._key = None
        self._st = None
        self._Q = float(p.Q_max)
        self._kS, self._Rth_env, self._Cth = float(p.k_S), float(p.k_env) * float(p.R_th), float(p.C_th)
        self._load_heat = bool(p.load_heat)

    def state(self, x) -> dict:
        key = (float(x[0]), float(x[1]), float(x[2]))
        if key != self._key:
            if self._c is not None:
                self._st = _scalar_state(key[0], key[1], key[2], self.inp, self.p, self._c)
            else:
                self._st = {k: float(v) for k, v in fast_state(np.array(key), self.inp, self.p).items()}
            self._key = key
        return self._st

    def rhs(self, t, x):
        st = self.state(x)
        z, v1, T = float(x[0]), float(x[1]), float(x[2])
        I = st["I"]
        q = (I * st["V"] if self._load_heat else 0.0) + I * I * st["R0"] + v1 * v1 / st["R1"] \
            - I * T * self._kS * st["dUdT"]
        return [-I / (3600.0 * self._Q), (st["R1"] * I - v1) / st["tau1"],
                (q - (T - self.Tamb) / self._Rth_env) / self._Cth]

    def gvals(self, x) -> dict:
        """Raw event functions g1..g5 at x (brownout only while active)."""
        st = self.state(x)
        p = self.p
        g = dict(charge=float(x[0]), cutoff=st["V"] - float(p.V_cut))
        if self.brownout_active:
            g["brownout"] = st["Vpk"] - float(p.V_sys)
        g["thermal"] = float(p.T_shut) - float(x[2])
        g["infeasible"] = st["D"]
        return g

    def first_violated(self, x) -> str | None:
        """Event at a switching instant: None if all g > 0, else :func:`classify` of the components."""
        g = self.gvals(x)
        if all(v > 0 for v in g.values()):
            return None
        comps = event_components(np.asarray(x, float)[:, None], Inputs(*[np.atleast_1d(v) for v in self.inp]),
                                 self.p)
        return str(classify(comps)[0])

    def events(self, z_e: float):
        evs = []
        for name in self.gvals(np.array([0.5, 0.0, self.Tamb])):
            f = (lambda nm: (lambda t, x: self.gvals(x)[nm]))(name)
            f.terminal, f.direction = True, -1
            f.__name__ = name
            evs.append(f)
        f = lambda t, x: x[0] - z_e  # displayed 0 % (non-terminal)
        f.terminal, f.direction = False, -1
        f.__name__ = "disp0"
        evs.append(f)
        return evs


# ------------------------------------------------------------------------------------------- results
@dataclass
class Result:
    """Outcome of one TTE run (§1.5, §6.1)."""

    tte_h: float
    event: str
    s_off: float  # displayed SOC (D1) at shutdown
    soe_off: float  # energy-fraction SOC (D2) at shutdown
    t_disp0_h: float  # time when D1 reaches 0 (NaN if shutdown happens first)
    traj: dict | None = None
    T_end_C: float = math.nan
    Ts_end_C: float = math.nan
    P_end_W: float = math.nan
    Psoc_end_W: float = math.nan
    V_end: float = math.nan
    z_end: float = math.nan
    P_mean_W: float = math.nan  # delivered energy / TTE
    rho_perf: float = math.nan  # time-averaged delivered SoC performance over SoC-demanding segments
    segments: int = 0

    def as_dict(self) -> dict:
        d = {f.name: getattr(self, f.name) for f in fields(self) if f.name != "traj"}
        d["s_off_pct"] = 100 * self.s_off
        d["soe_off_pct"] = 100 * self.soe_off
        return d


def _as_schedule(schedule) -> Schedule:
    if isinstance(schedule, Schedule):
        return schedule
    if isinstance(schedule, str):
        return Schedule.constant(PRIMITIVES[schedule])
    if isinstance(schedule, UsageState):
        return Schedule.constant(schedule)
    return Schedule(list(schedule))


def _trapz(y, t):
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(t))) if len(t) > 1 else 0.0


def run(schedule, Tamb_C: float, s0: float = 1.0, T0_C: float | None = None, p: Params = RP63,
        solver: str = "lsoda", h: float = 2.0, t_max_h: float = 2000.0, record: bool = False,
        rtol: float = 1e-7, atol=(1e-9, 1e-7, 1e-5), max_step: float = 60.0) -> Result:
    """Time-to-empty of one trajectory (§6.1–6.2).

    ``schedule``: a ``Schedule``, a list of ``(duration_s, UsageState)`` (cyclic), a ``UsageState`` or a primitive
    name.  ``Tamb_C``: ambient (°C; a segment's ``T_amb_C`` overrides it).  ``s0``: displayed SOC (D1) at start.
    ``T0_C``: initial temperature (default ambient of the first segment).  ``solver``: 'lsoda' (reference: segment
    by segment, events terminal with direction −1, rtol 1e-7, atol [1e-9, 1e-7, 1e-5], max_step 60 s) or 'rk4'
    (``run_batch`` with one sample and step ``h``).
    """
    sch = _as_schedule(schedule)
    if solver == "rk4":
        b = run_batch(sch, Tamb_C, None, s0=s0, h=h, T0_C=T0_C, p=p, t_max_h=t_max_h, record=record)
        return Result(float(b["tte_h"][0]), str(b["event"][0]), float(b["s_off"][0]), float(b["soe_off"][0]),
                      float(b["t_disp0_h"][0]), traj=b.get("traj"), T_end_C=float(b["T_end_C"][0]),
                      Ts_end_C=float(b["Ts_end_C"][0]), P_end_W=float(b["P_end_W"][0]), V_end=float(b["V_end"][0]),
                      z_end=float(b["z_end"][0]))
    if solver != "lsoda":
        raise ValueError(solver)
    ze = float(cell.z_empty(p))
    seg_T = [(Tamb_C if u.T_amb_C is None else u.T_amb_C) + K0 for _, u in sch]
    T0 = seg_T[0] if T0_C is None else T0_C + K0
    x = np.array([ze + s0 * (1.0 - ze), 0.0, T0])
    cache: dict[int, tuple[Segment, list]] = {}
    t_disp0 = math.nan
    E_out = 0.0
    rho_num = rho_den = 0.0
    tr = {k_: [] for k_ in ("t", "z", "v1", "T", "V", "I", "P", "Ts", "Psoc", "seg")} if record else None
    t_max = t_max_h * 3600.0

    def make_traj():
        if not record:
            return None
        traj = {k_: np.concatenate(v) if v else np.array([]) for k_, v in tr.items()}
        traj["t_h"] = traj["t"] / 3600.0
        traj["s"] = (traj["z"] - ze) / (1.0 - ze)
        traj["T_C"] = traj["T"] - K0
        traj["Ts_C"] = traj["Ts"] - K0
        return traj

    def finish(t_ev, ev, xe, seg, k):
        st = seg.state(xe)
        s_off = (xe[0] - ze) / (1.0 - ze)
        soe = float(cell.display_soc(xe[0], p, "D2"))
        traj = make_traj()
        return Result(t_ev / 3600.0, ev, float(s_off), soe, t_disp0, traj=traj, T_end_C=float(xe[2] - K0),
                      Ts_end_C=st["Ts"] - K0, P_end_W=st["P"], Psoc_end_W=st["Psoc"], V_end=st["V"],
                      z_end=float(xe[0]), P_mean_W=E_out / t_ev if t_ev > 0 else st["P"],
                      rho_perf=rho_num / rho_den if rho_den > 0 else math.nan, segments=k + 1)

    def problem(j, u, t_sw, x_sw):
        if j not in cache:
            sg = Segment(u, seg_T[j], p)
            cache[j] = (sg, sg.events(ze))
        seg, evs = cache[j]
        if seg.first_violated(x_sw) is not None:
            return None
        return seg.rhs, evs

    x_last, k_last = x, 0
    for k, j, t_sw, x_sw, sol in solve_ivp_segmented(sch.segments, x, t_max, problem, method="LSODA",
                                                      rtol=rtol, atol=list(atol), max_step=max_step):
        seg, evs = cache[j]
        k_last = k
        if sol is None:  # an event is violated at the switching instant (§6.1)
            return finish(t_sw, seg.first_violated(x_sw), x_sw, seg, k)
        u = sch[j][1]
        # statistics on the solver's own points (vectorised evaluation of the algebraic state)
        Y = sol.y
        stv = fast_state(Y, seg.inp, p)
        E_out += _trapz(np.broadcast_to(stv["P"], sol.t.shape), sol.t)
        if float(seg.inp.dem) > 0:
            rp = rho_perf(stv["Psoc"], u.x_c, p)
            rho_num += _trapz(np.broadcast_to(rp, sol.t.shape), sol.t)
            rho_den += sol.t[-1] - sol.t[0]
        if record:
            for k_, v in (("t", sol.t), ("z", Y[0]), ("v1", Y[1]), ("T", Y[2]), ("V", stv["V"]), ("I", stv["I"]),
                          ("P", stv["P"]), ("Ts", stv["Ts"]), ("Psoc", stv["Psoc"])):
                tr[k_].append(np.broadcast_to(v, sol.t.shape).astype(float))
            tr["seg"].append(np.full(sol.t.shape, j))
        hit = None
        for ev_f, te, ye in zip(evs, sol.t_events, sol.y_events):
            if len(te) == 0:
                continue
            if ev_f.__name__ == "disp0":
                if math.isnan(t_disp0):
                    t_disp0 = float(te[0]) / 3600.0
                continue
            if hit is None or te[0] < hit[0]:
                hit = (float(te[0]), ev_f.__name__, ye[0])
        if hit is not None:
            return finish(hit[0], hit[1], np.asarray(hit[2], float), seg, k)
        x_last = sol.y[:, -1]
    return Result(math.inf, "none", math.nan, math.nan, t_disp0, traj=make_traj(), z_end=float(x_last[0]),
                  T_end_C=float(x_last[2] - K0), P_mean_W=E_out / t_max, segments=k_last + 1)


# ------------------------------------------------------------------------------------------- batch driver
def _take(p: Params, idx: np.ndarray, n: int) -> Params:
    """Restrict array-valued parameter fields of length n to the samples ``idx``."""
    kw = {}
    for f in fields(p):
        v = getattr(p, f.name)
        if isinstance(v, np.ndarray) and v.ndim >= 1 and v.shape[-1] == n and f.name not in ("ocv_table",
                                                                                              "dudt_table"):
            kw[f.name] = v[..., idx]
    return replace(p, **kw) if kw else p


def _bcast(v, n: int) -> np.ndarray:
    return np.broadcast_to(np.asarray(v, dtype=float), (n,)).astype(float)


class _Driver:
    """Holds the per-step inputs for rk4_batch.  The event call at the end of each step advances the step
    index, so all four RK4 stages of a step use the input at the step start (§6.3).

    At a step whose end is an input switch the event function is evaluated twice (§6.1): first with the input
    the step was integrated with (a smooth crossing inside the step, e.g. the cut-off at the very end of a burst
    followed by a lighter segment), then with the new input (an event at the switching instant, whose time is
    then the switch time itself: the returned g is 0, so rk4_batch's interpolation gives frac = 1).
    [REVIEW 2026-11-17] Before this, only the new input was tested, so a crossing in the last step of a heavy
    segment was missed whenever the next segment was lighter (S9: +1–2.4 % TTE on ≈ 10 % of parameter samples).
    """

    def __init__(self, tables, step_index, h: float, p: Params):
        self.tables = tables  # dict of (m, n) arrays: Pns, dem, dPb, Tamb
        self.step_index = step_index  # callable k -> segment/state index, scalar or (n,)
        self.h = h
        self.p = p
        self.t0 = 0.0
        self.k = None
        self.idx = None
        self.cols = None
        self.inp = None
        self._ev_x = None
        self._ev_st = None

    def set_cols(self, cols, p_sub):
        self.cols = cols
        self.p = p_sub
        self.k = None
        self.idx = None
        self.inp = None
        self._ev_x = None
        self._p_vec = p_sub.n_samples() > 1  # array-valued parameters must be subset with the columns
        m = len(cols)
        self.alive = np.ones(m, dtype=bool)  # mirrors rk4_batch's crossing rule to record detection states
        self.det = np.full((8, m), np.nan)  # 5 event components + P, V, T_s at the detection step

    def inputs_for(self, idx, cols):
        T = self.tables
        if np.ndim(idx) == 0:
            return Inputs(T["Pns"][idx][cols], T["dem"][idx][cols], T["dPb"][idx][cols], T["Tamb"][idx][cols])
        return Inputs(T["Pns"][idx, cols], T["dem"][idx, cols], T["dPb"][idx, cols], T["Tamb"][idx, cols])

    def set_step(self, k: int):
        if k == self.k:
            return
        self.k = k
        idx = self.step_index(k)
        if np.ndim(idx) > 0:
            idx = idx[self.cols]
        self.idx = idx
        self.inp = self.inputs_for(idx, self.cols)

    def f(self, t, x):
        # rk4_batch evaluates k1 of the next step on the very array it just passed to the event function:
        # reuse that algebraic state (same x object, same step inputs).
        if x is self._ev_x:
            st, self._ev_x = self._ev_st, None
            return rhs_batch(x, self.inp, self.p, st)
        return rhs_batch(x, self.inp, self.p)

    def event(self, t, x):
        k_old, idx_old, inp_old = self.k, self.idx, self.inp
        self.set_step(int(round((self.t0 + t) / self.h)))
        st = fast_state(x, self.inp, self.p)
        self._ev_x, self._ev_st = x, st
        comps = event_components(x, self.inp, self.p, st)
        g = np.where(comps[4] > 0, np.min(comps, axis=0), np.nan)
        m = g.shape[0]
        extra = [np.broadcast_to(st[key], (m,)).copy() for key in ("P", "V", "Ts")]
        if k_old is not None and self.k != k_old:
            changed = np.broadcast_to(np.asarray(self.idx) != np.asarray(idx_old), (m,))
            if changed.any():
                cols = np.nonzero(changed)[0]
                full = cols.size == m
                xs = x if full else x[:, cols]
                inp_s = inp_old if full else Inputs(*(np.broadcast_to(a, (m,))[cols] for a in inp_old))
                ps = self.p if (full or not self._p_vec) else _take(self.p, cols, m)
                st_o = fast_state(xs, inp_s, ps)
                c_o = event_components(xs, inp_s, ps, st_o)
                g_o = np.where(c_o[4] > 0, np.min(c_o, axis=0), np.nan)
                smooth = ~(g_o > 0)  # crossed during the step, under the input it was integrated with
                at_switch = ~smooth & ~(g[cols] > 0)  # violated only by the new input: TTE = switch time
                g[cols[at_switch]] = 0.0
                cs = cols[smooth]
                if cs.size:
                    g[cs] = g_o[smooth]
                    comps[:, cs] = c_o[:, smooth]
                    for e, key in zip(extra, ("P", "V", "Ts")):
                        e[cs] = np.broadcast_to(st_o[key], (cols.size,))[smooth]
        new = self.alive & ~(g > 0)  # same test as rk4_batch (g <= 0 or non-finite)
        if new.any():
            self.det[:5, new] = comps[:, new]
            for i, e in enumerate(extra, 5):
                self.det[i, new] = e[new]
            self.alive &= ~new
        return g


def _input_tables(us: list[UsageState], p: Params, Tamb_C, n: int) -> dict:
    Pns, dem, dPb, Ta = [], [], [], []
    for u in us:
        a, b, c, _ = power_parts(u, p)
        Pns.append(_bcast(a, n))
        dem.append(_bcast(b, n))
        dPb.append(_bcast(c, n))
        Ta.append(_bcast((Tamb_C if u.T_amb_C is None else u.T_amb_C), n) + K0)
    return dict(Pns=np.stack(Pns), dem=np.stack(dem), dPb=np.stack(dPb), Tamb=np.stack(Ta))


def run_batch(schedules_or_ctmc, Tamb_C, theta: dict | None = None, s0=1.0, h: float = 2.0, T0_C=None,
              p: Params = RP63, t_max_h: float | None = None, n: int | None = None, paths: Paths | None = None,
              seed: int = 7, chunk_s: float = 3600.0, record: bool = False) -> dict:
    """Batch TTE with ``rk4_batch`` over parameter samples (§6.3).

    Parameters
    ----------
    schedules_or_ctmc : ``Schedule`` (or anything ``run`` accepts) shared by all samples, or a ``CTMC`` (each
        sample gets its own path; pass ``paths`` to reuse sampled paths).
    Tamb_C, s0, T0_C : scalars or per-sample arrays (n,).
    theta : dict name -> array (n,) of parameter values replacing fields of ``p``.
    h : RK4 step [s]; deterministic segment durations should be multiples of h.
    t_max_h : horizon (default 2000 h for schedules, 240 h for CTMC); samples without event get TTE = inf.
    n : number of samples when it cannot be inferred from theta/arrays (e.g. CTMC with no theta).
    chunk_s : samples that have died are dropped from the integration every ``chunk_s`` seconds.
    record : keep the full RK4 history (``traj = {t, x (steps, 3, n)}``); intended for small n.

    Returns dict of arrays: ``tte_h``, ``event`` (str), ``s_off`` (D1), ``soe_off`` (D2), ``z_end``,
    ``T_end_C``, ``Ts_end_C``, ``P_end_W``, ``V_end``, ``t_disp0_h`` (NaN: not tracked in batch) and, for CTMC,
    ``paths``.
    """
    pb = with_theta(p, theta)
    # RK4 is stable on the RC mode only for h < 2.785*tau1.  tau1 shrinks when hot (Arrhenius) and UQ samples
    # tau1_ref down to 10 s, so check the smallest tau1 at a conservative hot-cell temperature.
    T_hot = float(np.max(Tamb_C)) + 25.0 + 273.15
    tau_min = float(np.min(cell.tau1(T_hot, pb)))
    if h > 2.5 * tau_min:
        warnings.warn(f"run_batch: h = {h} s exceeds 2.5*tau1_min = {2.5 * tau_min:.1f} s at "
                      f"{T_hot - 273.15:.0f} C; RK4 may be inaccurate or unstable for some samples")
    sizes = [pb.n_samples(), np.size(Tamb_C), np.size(s0), np.size(T0_C) if T0_C is not None else 1,
             n or 1, paths.n if paths is not None else 1]
    N = max(sizes)
    is_ctmc = isinstance(schedules_or_ctmc, CTMC)
    if t_max_h is None:
        t_max_h = 240.0 if is_ctmc else 2000.0
    if is_ctmc:
        ctmc = schedules_or_ctmc
        if paths is None:
            paths = ctmc.sample_paths(N, t_max_h, h=h, seed=seed)
        elif paths.n != N:
            raise ValueError(f"{paths.n} paths supplied for {N} samples")
        paths.h = h
        us = ctmc.primitives()
        tables = _input_tables(us, pb, Tamb_C, N)
        cursor = paths.cursor()
        step_index = cursor.state_for_step
        T_first = tables["Tamb"][paths.states[:, 0], np.arange(N)]
    else:
        sch = _as_schedule(schedules_or_ctmc)
        rem = np.mod(sch.durations / h, 1.0)
        if np.any(np.minimum(rem, 1 - rem) > 1e-9):
            warnings.warn("segment durations are not multiples of h; switches are rounded to the step grid")
        us = [u for _, u in sch]
        tables = _input_tables(us, pb, Tamb_C, N)
        step_index = lambda k: int(sch.index_at(k * h + 1e-9 * h))
        T_first = tables["Tamb"][0]
    ze = _bcast(cell.z_empty(pb), N)
    z0 = ze + _bcast(s0, N) * (1.0 - ze)
    T0 = T_first.copy() if T0_C is None else _bcast(T0_C, N) + K0
    x = np.stack([z0, np.zeros(N), T0])

    tte = np.full(N, np.inf)
    xf = x.copy()
    det = np.full((8, N), np.nan)
    alive_idx = np.arange(N)
    drv = _Driver(tables, step_index, h, pb)
    # dead samples are compacted away every chunk; with record=True keep one chunk so histories align
    chunk = h * max(1, int(round(chunk_s / h))) if not record else t_max_h * 3600.0
    t0 = 0.0
    t_end = t_max_h * 3600.0
    hist_t, hist_x = [], []
    while alive_idx.size and t0 < t_end - 1e-9:
        drv.set_cols(alive_idx, _take(pb, alive_idx, N))
        drv.t0 = t0
        L = min(chunk, t_end - t0)
        te, xe, hist = rk4_batch(drv.f, x[:, alive_idx], L, h, event=drv.event, record_every=1 if record else 0)
        if record:
            hist_t.append(t0 + hist[0])
            hist_x.append(hist[1])
        died = np.isfinite(te)
        di = alive_idx[died]
        tte[di] = t0 + te[died]
        xf[:, di] = xe[:, died]
        det[:, di] = drv.det[:, died]
        x[:, alive_idx] = xe
        alive_idx = alive_idx[~died]
        t0 += L
    xf[:, alive_idx] = x[:, alive_idx]
    dead = np.isfinite(tte)
    ev = np.full(N, "none", dtype=object)
    if dead.any():
        ev[dead] = classify(det[:5, dead])
    out = dict(tte_h=tte / 3600.0, event=ev.astype(str), z_end=xf[0], T_end_C=xf[2] - K0)
    out["s_off"] = np.where(dead, (xf[0] - ze) / (1.0 - ze), np.nan)
    out["soe_off"] = np.where(dead, cell.display_soc(xf[0], pb, "D2") + 0.0 * ze, np.nan)
    out["P_end_W"], out["V_end"], out["Ts_end_C"] = det[5], det[6], det[7] - K0
    out["t_disp0_h"] = np.full(N, np.nan)
    if is_ctmc:
        out["paths"] = paths
    if record:
        out["traj"] = dict(t=np.concatenate(hist_t), x=np.concatenate(hist_x, axis=0))
    return out


__all__ = ["EVENTS", "Inputs", "fast_state", "rhs_from_state", "rhs_batch", "event_components", "event_batch",
           "classify", "Segment", "Result", "run", "run_batch"]
