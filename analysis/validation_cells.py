"""Cell-test drivers for the validation group (AN6: FINAL_SPEC §8 V2, V3, V4, V5).

The frozen core (``src/phonebatt``) integrates a phone as a constant-*power* load (A1).  The open cell datasets
are constant-current (HG2, NASA), a logged current profile (Oxford) or constant power (PyBaMM runs), so this module
integrates the same cell equations under a prescribed load.  Every constitutive function is the core's own
(``phonebatt.cell``: OCV with the entropic term, R0(z, T), R1(T), tau1(T), a_R, a_deg); only the load closure and
the structural variants are local.

Variants (FINAL_SPEC §1.1, §10.4)
---------------------------------
* ``L1``  quasi-static ECM: V = U(z, T) − I·(R0 + R1)                                   (X2 when isothermal)
* ``L2``  1RC Thevenin (E4): V = U − v1 − I·R0, dv1/dt = (R1·I − v1)/tau1               (baseline)
* ``X5``  2RC: V = U − v1 − v2 − I·R0, each branch with its own R_k, tau_k and Arrhenius energies
* ``X9``  KiBaM bound charge instead of the RC state (PS M2-K): dδ/dt = I/(3600 c) − k'(T)·δ  [δ in A h],
          z_surf = z − (1 − c)·δ/Q, V = U(z_surf, T) − I·(R0 + R1).  k' = κ_K(T)·(1 − c)/c with
          κ_K(T) = κ_K·a_deg(T; E_a,K) (a rate constant: slower when colder).  R1 is kept as a quasi-static series
          resistance so that the d.c. resistance budget equals L2's; X9 then differs from L2 only by the
          diffusion-limited availability (the question V2 asks).  ``x9_r1=False`` drops R1 (sensitivity).

Heat (E7, Bernardi): q = I²R0 + Σ v_k²/R_k (dynamic variants) or I²(R0 + R1) (L1, X9) + I·(U(z) − U(z_surf)) (X9)
− I·T·k_S·∂U/∂T(z), plus the delivered power I·V only if ``p.load_heat`` (phones; False for cells on a tester).
``isothermal=True`` holds T at the ambient (the PyBaMM DFN reference is isothermal).

Loads: a list of segments ``(duration_s, kind, value)`` with kind 'I' (A, > 0 discharge) or 'P' (W), integrated
segment by segment with LSODA (rtol 1e-7, as the core's reference solver), or a logged current profile with a
fixed-step RK4 and zero-order hold (``simulate_profile``).  Events: V ≤ V_cut ('cutoff'), z ≤ 0 ('charge'),
D ≤ 0 ('infeasible', constant power only).  A segment list can end with a long segment to run to the event.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

from phonebatt import cell
from phonebatt.params import K0, Params

# X9 (FINAL_SPEC §10.4): c = 0.166 (JH09), κ_K = 40 h^-1, E_a,K = 35 kJ/mol
X9_DEFAULT = dict(c=0.166, kappa_h=40.0, Ea=35e3)


@dataclass
class SimResult:
    t: np.ndarray
    V: np.ndarray
    I: np.ndarray
    z: np.ndarray
    T_C: np.ndarray
    tte_s: float
    event: str
    Ah: float  # delivered charge (integral of I dt) up to the event / end
    Wh: float  # delivered energy (integral of I V dt)
    sols: list = field(default_factory=list)  # (t0, t1, dense solution, kind, value) per segment
    model: "CellModel | None" = None

    @property
    def V_mean(self) -> float:
        return self.Wh / self.Ah if self.Ah > 0 else math.nan

    def V_at(self, t) -> np.ndarray:
        """Terminal voltage at times t (s) from the dense segment solutions (NaN outside the simulated range)."""
        t = np.atleast_1d(np.asarray(t, float))
        out = np.full(t.shape, np.nan)
        for t0, t1, sol, kind, val in self.sols:
            m = (t >= t0 - 1e-9) & (t <= t1 + 1e-9) & np.isnan(out)  # a switch instant belongs to the earlier segment
            if m.any():
                y = sol(np.clip(t[m], t0, t1))
                out[m] = self.model.algebra(y, kind, val)["V"]
        return out


class CellModel:
    """One cell (or phone pack) under a prescribed load; see the module docstring."""

    def __init__(self, p: Params, variant: str = "L2", isothermal: bool = False, x5: dict | None = None,
                 x9: dict | None = None, x9_r1: bool = True):
        if variant not in ("L1", "L2", "X5", "X9"):
            raise ValueError(variant)
        self.p, self.variant, self.isothermal = p, variant, isothermal
        self.x5 = x5
        self.x9 = dict(X9_DEFAULT, **(x9 or {}))
        self.x9_r1 = x9_r1
        if variant == "X5" and x5 is None:
            raise ValueError("X5 needs branch parameters x5=dict(R1, tau1, R2, tau2, Ea_R1, Ea_tau1, Ea_R2, Ea_tau2)")
        self.Q = float(p.Q_max)
        self.Vcut = float(p.V_cut)
        self.rg = 1.0 + float(p.rg)
        self.Ta = None  # set by simulate (K)

    # --------------------------------------------------------------------------------------- algebra
    def _resist(self, zc, T):
        p = self.p
        R0 = cell.R0(zc, T, p)
        if self.variant == "X5" and "R0" in self.x5:
            R0 = R0 * self.x5["R0"] / float(p.R0_ref)
        return R0

    def _branches(self, T):
        """(R, tau) of the RC branches at T."""
        p = self.p
        if self.variant == "X5":
            x = self.x5
            return [(x["R1"] * cell.a_R(T, x["Ea_R1"]) * self.rg, x["tau1"] * cell.a_R(T, x["Ea_tau1"])),
                    (x["R2"] * cell.a_R(T, x["Ea_R2"]) * self.rg, x["tau2"] * cell.a_R(T, x["Ea_tau2"]))]
        return [(cell.R1(T, p), cell.tau1(T, p))]

    def _kprime(self, T):
        c = self.x9["c"]
        kap = self.x9["kappa_h"] * cell.a_deg(T, self.x9["Ea"]) / 3600.0  # s^-1
        return kap * (1.0 - c) / c

    def algebra(self, y, kind: str, val: float) -> dict:
        """Algebraic state for state array y = (z, a1, a2, T) of shape (4,) or (4, m)."""
        z, a1, a2, T = y[0], y[1], y[2], y[3]
        if self.isothermal:
            T = np.full_like(np.asarray(z, float), self.Ta)
        zc = np.clip(z, 0.0, 1.0)
        U = cell.ocv(zc, T, self.p)
        R0 = self._resist(zc, T)
        br = self._branches(T)
        v = self.variant
        if v == "L1":
            Ue, Rin = U, R0 + br[0][0]
        elif v == "L2":
            Ue, Rin = U - a1, R0
        elif v == "X5":
            Ue, Rin = U - a1 - a2, R0
        else:  # X9
            zs = np.clip(z - (1.0 - self.x9["c"]) * a1 / self.Q, 0.0, 1.0)
            Ue = cell.ocv(zs, T, self.p)
            Rin = R0 + (br[0][0] if self.x9_r1 else 0.0)
        if kind == "I":
            I = np.broadcast_to(np.asarray(val, float), np.shape(Ue)).astype(float)
            D = np.ones_like(Ue)
            V = Ue - I * Rin
        else:
            I, V, D = cell.current(Ue, Rin, val)
        return dict(z=z, zc=zc, U=U, Ue=Ue, R0=R0, Rin=Rin, br=br, I=I, V=V, D=D, T=T)

    def rhs(self, t, y, kind, val):
        a = self.algebra(y, kind, val)
        I, T = a["I"], a["T"]
        v = self.variant
        dz = -I / (3600.0 * self.Q)
        da1 = da2 = 0.0
        if v == "L1":
            q = I * I * a["Rin"]
        elif v == "L2":
            R1, tau1 = a["br"][0]
            da1 = (R1 * I - y[1]) / tau1
            q = I * I * a["R0"] + y[1] ** 2 / R1
        elif v == "X5":
            (R1, t1), (R2, t2) = a["br"]
            da1 = (R1 * I - y[1]) / t1
            da2 = (R2 * I - y[2]) / t2
            q = I * I * a["R0"] + y[1] ** 2 / R1 + y[2] ** 2 / R2
        else:
            da1 = I / (3600.0 * self.x9["c"]) - self._kprime(T) * y[1]
            q = I * I * a["Rin"] + I * (a["U"] - a["Ue"])
        q = q - I * T * float(self.p.k_S) * cell.dUdT(a["zc"], self.p)
        if bool(self.p.load_heat):
            q = q + I * a["V"]
        if self.isothermal:
            dT = 0.0
        else:
            dT = (q - (T - self.Ta) / (float(self.p.k_env) * float(self.p.R_th))) / float(self.p.C_th)
        return [float(dz), float(da1), float(da2), float(dT)]

    # --------------------------------------------------------------------------------------- integration
    def simulate(self, segments, T_amb_C: float, z0: float = 1.0, T0_C: float | None = None,
                 rtol: float = 1e-7, max_step: float = 60.0, dt_rec: float | None = None) -> SimResult:
        """Integrate a list of load segments ``(duration_s, kind, value)`` from a rested state at z0."""
        self.Ta = T_amb_C + K0
        y = np.array([z0, 0.0, 0.0, (T_amb_C if T0_C is None else T0_C) + K0], float)
        t0 = 0.0
        ts, Vs, Is, zs, Ts, sols = [], [], [], [], [], []
        Ah = Wh = 0.0
        event = "none"
        atol = [1e-9, 1e-7, 1e-7, 1e-5]
        for dur, kind, val in segments:
            a0 = self.algebra(y, kind, val)
            if a0["V"] <= self.Vcut or (kind == "P" and a0["D"] <= 0):  # violated at the switch (§6.1)
                event = "infeasible" if (kind == "P" and a0["D"] <= 0) else "cutoff"
                break

            def ev_cut(t, yy, kind, val):  # solve_ivp passes ``args`` to the events as well
                return float(self.algebra(yy, kind, val)["V"]) - self.Vcut

            def ev_z(t, yy, kind, val):
                return float(yy[0])

            def ev_D(t, yy, kind, val):
                return float(self.algebra(yy, kind, val)["D"]) if kind == "P" else 1.0

            for e in (ev_cut, ev_z, ev_D):
                e.terminal, e.direction = True, -1
            sol = solve_ivp(self.rhs, (t0, t0 + dur), y, method="LSODA", rtol=rtol, atol=atol,
                            max_step=max_step, events=(ev_cut, ev_z, ev_D), args=(kind, val), dense_output=True)
            if not sol.success:
                raise RuntimeError(sol.message)
            t_end = float(sol.t[-1])
            if dt_rec:
                tg = np.unique(np.concatenate([np.arange(t0, t_end, dt_rec), [t_end]]))
                Y = sol.sol(tg)
            else:
                tg, Y = sol.t, sol.y
            al = self.algebra(Y, kind, val)
            I_, V_ = np.broadcast_to(al["I"], tg.shape), np.broadcast_to(al["V"], tg.shape)
            if len(tg) > 1:
                Ah += float(np.sum(0.5 * (I_[1:] + I_[:-1]) * np.diff(tg))) / 3600.0
                Wh += float(np.sum(0.5 * (I_[1:] * V_[1:] + I_[:-1] * V_[:-1]) * np.diff(tg))) / 3600.0
            ts.append(tg), Vs.append(V_), Is.append(I_), zs.append(Y[0]), Ts.append(np.broadcast_to(al["T"], tg.shape))
            sols.append((t0, t_end, sol.sol, kind, val))
            y = sol.y[:, -1].copy()
            t0 = t_end
            hit = [i for i, te in enumerate(sol.t_events) if len(te)]
            if hit:
                event = ("cutoff", "charge", "infeasible")[min(hit, key=lambda i: sol.t_events[i][0])]
                break
        cat = lambda L: np.concatenate(L) if L else np.array([])
        return SimResult(cat(ts), cat(Vs), cat(Is), cat(zs), cat(Ts) - K0, t0 if event != "none" else math.inf,
                         event, Ah, Wh, sols=sols, model=self)

    def simulate_profile(self, t_data, I_data, T_amb_C: float, z0: float = 1.0, T_data_C=None,
                         h: float = 0.25) -> dict:
        """Logged current profile (zero-order hold between samples) with fixed-step RK4.

        If ``T_data_C`` is given the temperature is prescribed from the log (thermal ODE bypassed); the voltage is
        returned at the sample times, computed with the sample's own current (simultaneous V and I logging).
        """
        self.Ta = T_amb_C + K0
        y = np.array([z0, 0.0, 0.0, (T_amb_C if T_data_C is None else T_data_C[0]) + K0], float)
        n = len(t_data)
        V = np.empty(n)
        zz = np.empty(n)
        iso = self.isothermal
        if T_data_C is not None:  # prescribed temperature: hold T at the logged value over each interval
            self.isothermal = True
        try:
            for k in range(n):
                if T_data_C is not None:
                    self.Ta = T_data_C[k] + K0
                    y[3] = self.Ta
                a = self.algebra(y, "I", I_data[k])
                V[k], zz[k] = float(a["V"]), y[0]
                if k == n - 1:
                    break
                dt = t_data[k + 1] - t_data[k]
                m = max(1, int(round(dt / h)))
                hh = dt / m
                f = lambda yy: np.array(self.rhs(0.0, yy, "I", I_data[k]))
                for _ in range(m):
                    k1 = f(y)
                    k2 = f(y + 0.5 * hh * k1)
                    k3 = f(y + 0.5 * hh * k2)
                    k4 = f(y + hh * k3)
                    y = y + hh / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        finally:
            self.isothermal = iso
        return dict(t=np.asarray(t_data, float), V=V, z=zz)


def run_to_event(model: CellModel, kind: str, value: float, T_amb_C: float, z0: float = 1.0,
                 T0_C: float | None = None, t_max_s: float = 3.0e5, **kw) -> SimResult:
    """Constant current ('I') or constant power ('P') until the first event."""
    return model.simulate([(t_max_s, kind, value)], T_amb_C, z0=z0, T0_C=T0_C, **kw)


def rmse(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    return float(np.sqrt(np.mean((a[m] - b[m]) ** 2))) if m.any() else math.nan


__all__ = ["CellModel", "SimResult", "run_to_event", "rmse", "X9_DEFAULT"]
