"""Structural variants X2–X14, X18 of FINAL_SPEC §10.4 as one vectorised "union" model (AN10 helper).

The frozen core (``src/phonebatt``) integrates the baseline L2/L3 system x = (z, v1, T).  The structural variants
change the model form (state dimension, current law, resistance law, radio law, ...), so they are implemented here
as a superset of the core equations whose extra terms vanish for the baseline X0.  Every constitutive function is
the core's own (``phonebatt.cell``: OCV incl. entropic term, R0/R1/τ1 Arrhenius forms, constant-power root,
burst peak; ``phonebatt.power``: E12–E15 loads and the E9 throttle); only the variant terms are new.

State per column: y = (z, v1, T, v2, δ)
* v2  second RC branch (X5; R2 = 0 elsewhere, so v2 ≡ 0),
* δ   KiBaM bound charge [A h] (X9; ≡ 0 elsewhere),
* the RRC tail occupancy r of X10 is not an RK4 state: its ODE dr/dt = λ_c(1 − r) − r/T_tail does not depend on x,
  so it is advanced exactly (exponential solution) inside each step; the RK4 stages see r(t) exactly.

Per-column switches (``Variant``):
* ``series``  L1-type quasi-static ECM: Ũ = U, series R0 + R1, no RC state (X2, and the series resistance of X9);
* ``iso``     isothermal at T_amb (X2, X3); ``throttle=False`` disables E9 (T_lim → 10⁶ °C); ``burst=False`` (X7);
* ``cc``      constant-current load I = P/V_nom, V_nom = 3.87 V (X4);
* ``pcm_frac`` temperature-independent fraction of R0 (X6: 35/60);
* ``rc2``     2RC (R_a, τ_a, R_b, τ_b) at 25 °C replacing (R1, τ1) (X5);
* ``kibam``   KiBaM bound charge (c, κ_K, E_a,K) replacing the RC state, R1 kept as series resistance (X9);
* ``tail_ode`` cellular RRC tail as the ODE state r (X10);
* ``k_lk``    SoC leakage 0.1 W·[exp(k_lk(T_s − 25 °C)) − 1] while awake (X13);
* ``params``  Params overrides (X8 d_sys = 0, X11 γ_b, X18 OCV stretch);
* ``input_map`` a map on every UsageState (X14 audio on a DSP: awake = 0 in AU/AU-L);
* ``soc_mode`` 'D1' or 'D2' for s0 and s_off (X12).

Inputs: every column has its own piecewise-constant input given as a ``phonebatt.usage.Paths`` (jump times and the
state after each jump) over its own list of UsageStates; deterministic cyclic schedules are expanded into such a
path.  The step convention and the event logic are those of ``phonebatt.model.run_batch`` (§6.3, including the
two evaluations of g at a switching step, [REVIEW 2026-11-17]); with all switches off the union model reproduces
``run_batch`` (checked in run_structural.py).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import Callable

import numpy as np

from phonebatt import cell
from phonebatt.cell import DUDT_V, DUDT_Z, current, peak_voltage
from phonebatt.integrate import rk4_batch
from phonebatt.model import classify
from phonebatt.params import K0, RG, RP63, T_REF, Params, aged, with_theta
from phonebatt.power import m_sig, power_parts, throttle
from phonebatt.usage import Schedule, UsageState

V_NOM = 3.87  # X4 constant-current conversion voltage (mean OCV)
X9_DEFAULT = dict(c=0.166, kappa_h=40.0, Ea=35e3)  # FINAL_SPEC §10.4 (same as analysis/validation_cells.py)
RT_RATE = 1.0e3  # 1/s: real-time traffic keeps the cellular radio connected (r -> 1 within ms)


@dataclass
class Variant:
    id: str
    label: str = ""
    series: bool = False
    iso: bool = False
    throttle: bool = True
    burst: bool = True
    cc: bool = False
    pcm_frac: float = 0.0
    rc2: tuple | None = None  # (R_a, tau_a, R_b, tau_b) at 25 °C [Ohm, s]
    kibam: dict | None = None
    tail_ode: bool = False
    k_lk: float = 0.0
    params: dict = field(default_factory=dict)
    input_map: Callable[[UsageState], UsageState] | None = None
    soc_mode: str = "D1"


@dataclass
class Column:
    """One (variant, scenario) run: inputs, path, environment and initial condition."""

    variant: Variant
    key: str  # scenario key
    inputs: list  # list[UsageState] (state space of the path)
    paths_times: np.ndarray  # (M,) jump times [s] (times[0] = 0), +inf padded
    paths_states: np.ndarray  # (M,) state index after each jump
    Tamb_C: float = 25.0
    soh: float = 1.0
    T0_C: float | None = None
    s0: float = 1.0
    p: Params = RP63
    table_key: str | None = None  # columns with equal keys share input tables / initial state (CTMC paths)


def schedule_path(sch: Schedule, t_max_s: float) -> tuple[np.ndarray, np.ndarray]:
    """Expand a cyclic schedule into jump times and segment indices up to t_max."""
    if len(sch) == 1:
        return np.array([0.0]), np.array([0])
    n_cyc = int(math.ceil(t_max_s / sch.period)) + 1
    starts = np.concatenate([[0.0], sch.ends[:-1]])
    t = (starts[None, :] + sch.period * np.arange(n_cyc)[:, None]).ravel()
    s = np.tile(np.arange(len(sch)), n_cyc)
    keep = t <= t_max_s
    return t[keep], s[keep]


# ------------------------------------------------------------------------------------------- column setup
def _column_params(col: Column) -> Params:
    v = col.variant
    p = col.p
    if col.soh < 1.0:
        p = aged(p, 1.0 - col.soh)
    if v.params:
        p = replace(p, **v.params)
    return p


def _tail_aux(u: UsageState, p: Params) -> tuple:
    """X10 per-input constants for the cellular radio: (lam_c, Tt_c, dPt_c, lam_pro, Epro_c, qs_c).

    lam_c: rate driving r (bursts on the cellular data path; RT_RATE during real-time traffic; 0 otherwise);
    Tt_c: cellular tail timer; dPt_c = m (P_tail − P_idle) is the tail power per unit occupancy (0 during real-time
    traffic, whose m·P_rt already is the connected-radio power; after the call r ≈ 1 decays as a tail);
    lam_pro·(1 − r)·Epro_c is the promotion power; qs_c: the quasi-steady tail + promotion power already contained
    in P_nonsoc (E14), which X10 replaces.
    """
    if u.c == "none":
        return 0.0, 1.0, 0.0, 0.0, 0.0, 0.0
    r = p.radio(u.c)
    m = float(m_sig(u.RSRP, p))
    Tt = float(r["T_tail"])
    dPt = m * float(r["P_tail"] - r["P_idle"])
    Ep = m * float(r["E_pro"])
    if u.rho == "cell":
        if u.kappa == "rt":
            return RT_RATE, Tt, 0.0, 0.0, Ep, 0.0
        lam = float(u.lam)
        qs = dPt * (1.0 - math.exp(-lam * Tt)) + lam * math.exp(-lam * Tt) * Ep
        return lam, Tt, dPt, lam, Ep, qs
    return 0.0, Tt, dPt, 0.0, Ep, 0.0


class UnionBatch:
    """Vectorised union model over columns (variants × scenarios) with an rk4_batch driver."""

    FIELDS = ("Q_n", "q_loss", "R0_ref", "A_R", "z0", "Ea_R0", "R1_ref", "tau1_ref", "Ea_R1", "Ea_tau", "V_cut",
              "d_sys", "C_th", "R_th", "R_skin", "T_lim_C", "dT_thr", "w_thr", "P_fl", "k_env", "k_S", "dU", "b0",
              "b1", "a0", "a1", "a2", "a3", "dV_hv", "p_hv", "T_shut_C", "gamma_R")

    def __init__(self, columns: list[Column], h: float):
        self.cols = columns
        self.h = float(h)
        n = len(columns)
        self.n = n
        ps = [_column_params(c) for c in columns]
        # array-valued Params over the columns (fields that may differ); r_growth from gamma_R*q_loss
        th = {f: np.array([float(getattr(pp, f)) for pp in ps]) for f in self.FIELDS}
        for pp in ps:
            if pp.r_growth is not None or pp.ocv_table is not None or pp.dudt_table is not None:
                raise ValueError("union model: explicit r_growth / tables not supported")
        thr = np.array([c.variant.throttle for c in columns])
        th["T_lim_C"] = np.where(thr, th["T_lim_C"], 1.0e6)
        self.p = with_theta(RP63, th)
        V = [c.variant for c in columns]
        f = lambda a: np.array(a, dtype=float)
        self.cfg = dict(
            series=f([v.series or v.kibam is not None for v in V]) > 0,
            kib=f([v.kibam is not None for v in V]) > 0,
            iso=f([v.iso for v in V]) > 0,
            cc=f([v.cc for v in V]) > 0,
            pcm=f([v.pcm_frac for v in V]),
            R2=f([0.0 if v.rc2 is None else v.rc2[2] for v in V]),
            tau2=f([1.0 if v.rc2 is None else v.rc2[3] for v in V]),
            kc=f([(v.kibam or X9_DEFAULT)["c"] for v in V]),
            kkap=f([(v.kibam or X9_DEFAULT)["kappa_h"] for v in V]),
            kEa=f([(v.kibam or X9_DEFAULT)["Ea"] for v in V]),
            tail=f([v.tail_ode for v in V]) > 0,
            klk=f([v.k_lk for v in V]),
        )
        # X5 replaces R1, tau1 by the first branch
        R1 = self.p.R1_ref.copy()
        t1 = self.p.tau1_ref.copy()
        for j, v in enumerate(V):
            if v.rc2 is not None:
                R1[j], t1[j] = v.rc2[0], v.rc2[1]
        self.p = replace(self.p, R1_ref=R1, tau1_ref=t1)
        self.any = {k: bool(np.any(self.cfg[k])) for k in ("kib", "tail", "cc")}
        self.any["R2"] = bool(np.any(self.cfg["R2"] > 0))
        self.any["klk"] = bool(np.any(self.cfg["klk"] != 0))
        self.any["pcm"] = bool(np.any(self.cfg["pcm"] != 0))
        # input tables (S_max, n)
        S = max(len(c.inputs) for c in columns)
        keys = ("Pns", "dem", "dPb", "Tamb", "awake", "lam_c", "Tt_c", "dPt_c", "lam_pro", "Epro_c", "qs_c")
        T = {k: np.zeros((S, n)) for k in keys}
        T["Tt_c"][:] = 1.0
        done: dict = {}
        for j, (c, pp) in enumerate(zip(columns, ps)):
            v = c.variant
            tk = None if c.table_key is None else (v.id, c.table_key)
            if tk is not None and tk in done:
                jj = done[tk]
                for k in keys:
                    T[k][:, j] = T[k][:, jj]
                continue
            if tk is not None:
                done[tk] = j
            for s, u in enumerate(c.inputs):
                uu = v.input_map(u) if v.input_map is not None else u
                a, b, d, _ = power_parts(uu, pp)
                T["Pns"][s, j] = float(a)
                T["dem"][s, j] = float(b)
                T["dPb"][s, j] = float(d) if v.burst else 0.0
                T["Tamb"][s, j] = (c.Tamb_C if uu.T_amb_C is None else uu.T_amb_C) + K0
                T["awake"][s, j] = float(np.maximum(uu.s_scr, uu.awake))
                if v.tail_ode:
                    (T["lam_c"][s, j], T["Tt_c"][s, j], T["dPt_c"][s, j], T["lam_pro"][s, j], T["Epro_c"][s, j],
                     T["qs_c"][s, j]) = _tail_aux(uu, pp)
        self.tab = T
        # paths
        M = max(c.paths_times.size for c in columns)
        self.times = np.full((n, M + 1), np.inf)
        self.states = np.zeros((n, M + 1), dtype=np.int64)
        for j, c in enumerate(columns):
            self.times[j, :c.paths_times.size] = c.paths_times
            self.states[j, :c.paths_states.size] = c.paths_states
            self.states[j, c.paths_states.size:] = c.paths_states[-1]
        self.M = M + 1
        # initial conditions
        ze = cell.z_empty(self.p)
        self.ze = np.broadcast_to(ze, (n,)).astype(float)
        z0 = np.empty(n)
        T0 = np.empty(n)
        zd2: dict = {}
        for j, c in enumerate(columns):
            if c.variant.soc_mode == "D2":
                kk = (c.variant.id, c.table_key if c.table_key is not None else j, c.s0)
                if kk not in zd2:
                    zd2[kk] = float(cell.soc_to_z(c.s0, ps[j], "D2"))
                z0[j] = zd2[kk]
            else:
                z0[j] = self.ze[j] + c.s0 * (1.0 - self.ze[j])
            Tfirst = T["Tamb"][self.states[j, 0], j]
            T0[j] = Tfirst if (c.T0_C is None or c.variant.iso) else c.T0_C + K0
        self.x0 = np.stack([z0, np.zeros(n), T0, np.zeros(n), np.zeros(n)])

    # --------------------------------------------------------------------------------- algebra
    def algebra(self, x, inp: dict, p: Params, cfg: dict, r=None) -> dict:
        z, v1, T, v2, dK = x[0], x[1], x[2], x[3], x[4]
        zc = np.minimum(np.maximum(z, 0.0), 1.0)
        dudt = np.interp(zc, DUDT_Z, DUDT_V)
        U = cell.ocv_ref(zc, p, clamped=True) + (T - T_REF) * p.k_S * dudt
        kib = cfg["kib"]
        if self.any["kib"] and kib.any():
            zs = np.minimum(np.maximum(z - (1.0 - cfg["kc"]) * dK / p.Q_max, 0.0), 1.0)
            Us = np.where(kib, cell.ocv_ref(zs, p, clamped=True) + (T - T_REF) * p.k_S * np.interp(zs, DUDT_Z,
                                                                                                   DUDT_V), U)
        else:
            Us = U
        ar = (1.0 / T - 1.0 / T_REF) / RG
        g = 1.0 + p.rg
        rinf = p.R0_ref / (1.0 + p.A_R * np.exp(-0.5 / p.z0))
        aR0 = np.exp(p.Ea_R0 * ar)
        if self.any["pcm"]:
            aR0 = cfg["pcm"] + (1.0 - cfg["pcm"]) * aR0
        R0_ = rinf * (1.0 + p.A_R * np.exp(-np.maximum(z, 0.0) / p.z0)) * aR0 * g
        aR1 = np.exp(p.Ea_R1 * ar)
        aT = np.exp(p.Ea_tau * ar)
        R1_ = p.R1_ref * aR1 * g
        tau = p.tau1_ref * aT
        Ts = inp["Tamb"] + p.R_skin / p.R_th * (T - inp["Tamb"])
        Psoc = throttle(inp["dem"], Ts, p)
        P = inp["Pns"] + Psoc
        if self.any["klk"]:
            P = P + inp["awake"] * 0.1 * (np.exp(cfg["klk"] * (Ts - T_REF)) - 1.0)
        if self.any["tail"] and r is not None:
            P = P + np.where(cfg["tail"], inp["dPt_c"] * r + inp["lam_pro"] * (1.0 - r) * inp["Epro_c"]
                             - inp["qs_c"], 0.0)
        ser = cfg["series"]
        vsum = v1 + v2
        Ue = np.where(ser, Us, U - vsum)
        Rin = np.where(ser, R0_ + R1_, R0_)
        I, V, D = current(Ue, Rin, P)
        Ue_pk = np.where(ser, Us - R1_ * I, Ue)
        Ppk = P + inp["dPb"]
        Vpk = peak_voltage(Ue_pk, R0_, Ppk)
        if self.any["cc"]:
            cc = cfg["cc"]
            Icc = P / V_NOM
            I = np.where(cc, Icc, I)
            V = np.where(cc, Ue - Rin * Icc, V)
            D = np.where(cc, Ue * Ue, D)
            Vpk = np.where(cc, Ue_pk - R0_ * Ppk / V_NOM, Vpk)
        out = dict(U=U, Us=Us, R0=R0_, R1=R1_, tau1=tau, dUdT=dudt, Ts=Ts, Psoc=Psoc, P=P, Ue=Ue, I=I, V=V, D=D,
                   Vpk=Vpk, aR1=aR1, aT=aT)
        return out

    def rhs(self, x, st, inp, p, cfg):
        z, v1, T, v2, dK = x
        I = st["I"]
        ser = cfg["series"]
        out = np.empty_like(x)
        out[0] = -I / (3600.0 * p.Q_max)
        out[1] = np.where(ser, 0.0, (st["R1"] * I - v1) / st["tau1"])
        q = I * st["V"] + I * I * st["R0"] - I * T * p.k_S * st["dUdT"]
        rc_heat = v1 * v1 / st["R1"]
        if self.any["R2"]:
            R2 = cfg["R2"] * st["aR1"] * (1.0 + p.rg)
            tau2 = cfg["tau2"] * st["aT"]
            has2 = (cfg["R2"] > 0) & ~ser
            out[3] = np.where(has2, (R2 * I - v2) / tau2, 0.0)
            rc_heat = rc_heat + np.where(has2, v2 * v2 / np.where(has2, R2, 1.0), 0.0)
        else:
            out[3] = 0.0
        q = q + np.where(ser, I * I * st["R1"], rc_heat)
        if self.any["kib"]:
            kib = cfg["kib"]
            kp = cfg["kkap"] * np.exp(-cfg["kEa"] / RG * (1.0 / T - 1.0 / T_REF)) / 3600.0 * (1.0 - cfg["kc"]) / \
                cfg["kc"]
            out[4] = np.where(kib, I / (3600.0 * cfg["kc"]) - kp * dK, 0.0)
            q = q + np.where(kib, I * (st["U"] - st["Us"]), 0.0)
        else:
            out[4] = 0.0
        dT = (q - (T - st["Tamb_"]) / (p.k_env * p.R_th)) / p.C_th
        out[2] = np.where(cfg["iso"], 0.0, dT)
        return out

    def components(self, x, st, inp, p):
        z, T = x[0], x[2]
        g3 = np.where(inp["dPb"] > 0, (st["Vpk"] - p.V_sys) / 0.1, np.inf)
        return np.stack(np.broadcast_arrays(z / 0.01, (st["V"] - p.V_cut) / 0.1, g3, (p.T_shut - T) / 10.0,
                                            st["D"] / (st["Ue"] * st["Ue"])))

    # --------------------------------------------------------------------------------- driver
    def run(self, t_max_h: float, chunk_s: float = 3600.0, progress: bool = False) -> dict:
        from phonebatt.model import _take
        n, h = self.n, self.h
        x = self.x0.copy()
        tte = np.full(n, np.inf)
        xf = x.copy()
        det = np.full((8, n), np.nan)
        alive = np.arange(n)
        ptr = np.zeros(n, dtype=np.int64)
        r_all = np.zeros(n)
        t0 = 0.0
        t_end = t_max_h * 3600.0
        chunk = h * max(1, int(round(chunk_s / h)))
        tab = self.tab
        model = self

        class Drv:
            pass

        while alive.size and t0 < t_end - 1e-9:
            cols = alive
            m = cols.size
            pc = _take(self.p, cols, n)
            cfg = {k: (v[cols] if isinstance(v, np.ndarray) else v) for k, v in self.cfg.items()}
            times = self.times[cols]
            states = self.states[cols]
            rows = np.arange(m)
            d = Drv()
            d.k = None
            d.idx = None
            d.inp = None
            d.r0 = r_all[cols].copy()
            d.ptr = ptr[cols].copy()
            d.alive = np.ones(m, bool)
            d.det = np.full((8, m), np.nan)
            d.ev_x = None
            d.ev_st = None

            def idx_for_step(k):
                tk = (k + 1.0) * h
                while True:
                    nxt = np.minimum(d.ptr + 1, self.M - 1)
                    adv = (d.ptr + 1 < self.M) & (times[rows, nxt] < tk)
                    if not adv.any():
                        break
                    d.ptr = np.where(adv, d.ptr + 1, d.ptr)
                return states[rows, d.ptr]

            def inputs(idx):
                return {k: tab[k][idx, cols] for k in tab}

            def r_at(s, inp):
                if not model.any["tail"]:
                    return None
                beta = inp["lam_c"] + 1.0 / inp["Tt_c"]
                rinf = inp["lam_c"] / beta
                return rinf + (d.r0 - rinf) * np.exp(-beta * s)

            def set_step(k):
                if k == d.k:
                    return
                d.k = k
                d.idx = idx_for_step(k)
                d.inp = inputs(d.idx)

            def state(xx, inp, s):
                st = model.algebra(xx, inp, pc, cfg, r_at(s, inp))
                st["Tamb_"] = inp["Tamb"]
                return st

            def f(t, xx):
                s = (t0 + t) - d.k * h
                if xx is d.ev_x and abs(s) < 1e-9:
                    st, d.ev_x = d.ev_st, None
                else:
                    st = state(xx, d.inp, s)
                return model.rhs(xx, st, d.inp, pc, cfg)

            def event(t, xx):
                k_new = int(round((t0 + t) / h))
                k_old, idx_old, inp_old = d.k, d.idx, d.inp
                if k_old is not None and model.any["tail"]:
                    d.r0 = r_at(h * (k_new - k_old), inp_old)  # exact r at the step end (old input)
                set_step(k_new)
                st = state(xx, d.inp, 0.0)
                d.ev_x, d.ev_st = xx, st
                comps = model.components(xx, st, d.inp, pc)
                g = np.where(comps[4] > 0, np.min(comps, axis=0), np.nan)
                extra = [np.broadcast_to(st[key], (m,)).copy() for key in ("P", "V", "Ts")]
                if k_old is not None and d.k != k_old:
                    changed = np.asarray(d.idx) != np.asarray(idx_old)
                    if changed.any():
                        st_o = state(xx, inp_old, 0.0)
                        c_o = model.components(xx, st_o, inp_old, pc)
                        g_o = np.where(c_o[4] > 0, np.min(c_o, axis=0), np.nan)
                        smooth = changed & ~(g_o > 0)
                        at_switch = changed & ~smooth & ~(g > 0)
                        g[at_switch] = 0.0
                        if smooth.any():
                            g[smooth] = g_o[smooth]
                            comps[:, smooth] = c_o[:, smooth]
                            for e, key in zip(extra, ("P", "V", "Ts")):
                                e[smooth] = np.broadcast_to(st_o[key], (m,))[smooth]
                new = d.alive & ~(g > 0)
                if new.any():
                    d.det[:5, new] = comps[:, new]
                    for i, e in enumerate(extra, 5):
                        d.det[i, new] = e[new]
                    d.alive &= ~new
                return g

            L = min(chunk, t_end - t0)
            te, xe, _ = rk4_batch(f, x[:, cols], L, h, event=event)
            died = np.isfinite(te)
            di = cols[died]
            tte[di] = t0 + te[died]
            xf[:, di] = xe[:, died]
            det[:, di] = d.det[:, died]
            x[:, cols] = xe
            ptr[cols] = d.ptr
            if model.any["tail"]:
                r_all[cols] = d.r0  # r at the chunk end (each chunk ends on a step boundary)
            alive = cols[~died]
            t0 += L
            if progress:
                print(f"  t = {t0 / 3600:.0f} h, alive {alive.size}", flush=True)
        xf[:, alive] = x[:, alive]
        dead = np.isfinite(tte)
        ev = np.full(n, "none", dtype=object)
        if dead.any():
            ev[dead] = classify(det[:5, dead])
        zf = xf[0]
        s1 = (zf - self.ze) / (1.0 - self.ze)
        s2 = np.asarray(cell.display_soc(zf, self.p, "D2"), float)
        mode = np.array([c.variant.soc_mode for c in self.cols])
        s_off = np.where(mode == "D2", s2, s1)
        return dict(tte_h=tte / 3600.0, event=ev.astype(str), s_off=np.where(dead, s_off, np.nan),
                    s_off_D1=np.where(dead, s1, np.nan), T_end_C=xf[2] - K0, z_end=zf)


__all__ = ["Variant", "Column", "UnionBatch", "schedule_path", "V_NOM", "X9_DEFAULT"]
