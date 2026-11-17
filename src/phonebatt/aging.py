"""L3 slow aging layer: squared-loss calendar, SEI-cycle and plating states (FINAL_SPEC §2.12, §11).

Each √-law loss q = k·√(stress-weighted exposure) is integrated as d(q²)/dt = k²·stress² (E17–E19), which is
exact over the piecewise-constant segments of the quasi-static day profile (§11.2).

State: ``AgingState(Q_cal, Q_H, Q_L)`` (squared losses).  q_loss = √Q_cal + √Q_H + √Q_L, SOH = 1 − q_loss,
r_growth = γ_R q_loss (E20).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, replace
from typing import Callable, NamedTuple

import numpy as np

from . import cell
from .cell import a_R, a_deg
from .params import K0, RP63, T_REF, Params


class AgingState(NamedTuple):
    Q_cal: float = 0.0
    Q_H: float = 0.0
    Q_L: float = 0.0

    @property
    def q_cal(self):
        return math.sqrt(self.Q_cal)

    @property
    def q_H(self):
        return math.sqrt(self.Q_H)

    @property
    def q_L(self):
        return math.sqrt(self.Q_L)

    @property
    def q_loss(self):
        return self.q_cal + self.q_H + self.q_L

    @property
    def soh(self):
        return 1.0 - self.q_loss


class Seg(NamedTuple):
    """One quasi-static aging segment: duration [h], mean SOC, cell temperature [K], FCE throughput (charge +
    discharge halves), charge throughput [FCE], charge C-rate [1/h]."""

    dt_h: float
    z: float
    T_K: float
    dn: float
    dn_ch: float
    c_ch: float


def _sig(y):
    return 1.0 / (1.0 + np.exp(-y))


def g_stress(z, p: Params = RP63):
    """Calendar SOC stress g(z) (E16): KEIL16 graphite step × SCH14 high-voltage factor h_HV^κ_HV.

    g = [0.7 + 0.3 σ((z − 0.2)/0.02) + (g_hi − 1) σ((z − 0.6)/0.03)] · h_HV(z)^κ_HV,
    h_HV = max{1, (7.543 U(z) − 23.75)/(7.543 U(0.8) − 23.75)} with U at T_ref.
    """
    z = np.asarray(z, dtype=float)
    step = 0.7 + 0.3 * _sig((z - 0.2) / 0.02) + (p.g_hi - 1.0) * _sig((z - 0.6) / 0.03)
    U = cell.ocv_ref(z, p)
    U08 = cell.ocv_ref(0.8, p)
    hv = np.maximum(1.0, (7.543 * U - 23.75) / (7.543 * U08 - 23.75))
    return step * hv ** p.k_HV


def jeita_rate(c, T_cell_C):
    """JEITA rule (TI-JEITA): c ≤ 0.5C at 0–10 °C cell temperature; no charging below 0 °C or above 60 °C."""
    if T_cell_C < 0.0 or T_cell_C > 60.0:
        return 0.0
    if T_cell_C < 10.0:
        return min(c, 0.5)
    return c


def day_profile(strategy: str, Tamb_C: float, dsoc: float, p: Params = RP63, c_cc: float = 1.0,
                c_cv: float = 0.3, z_cc: float = 0.8, P_day: float = 0.440, heat_h: float = 0.0,
                heat_C: float = 45.0, cal_mode: bool = False, Tch_amb_C: float | None = None,
                jeita: bool = False, dt_charge_h: float = 0.25) -> list[Seg]:
    """Quasi-static day (§11.2), starting 07:00 unplugged.

    * 07:00–23:00: SOC falls linearly by ``dsoc`` in 16 one-hour segments at T_use = T_amb + k_env R_th P_day
      (the first ``heat_h`` hours at ``heat_C``: car dashboard, A5).
    * 23:00: plug in; CC at ``c_cc`` up to ``z_cc`` then ``c_cv`` up to the target (E21) at
      T_ch = T_ch,amb + ΔT_ch; hold at the target at T_amb until 07:00.
    * Strategies: 'full' (target 1.0), 'limit80' (target and unplug SOC 0.8), 'obc' (charge to 0.8, hold,
      resume at ``c_cv`` 1 h before 07:00 to 1.0).
    * ``cal_mode``: protocol C0 (cell held at T_amb throughout: no use or charge heating).
    * ``jeita``: apply the JEITA current limits to the charging temperature (off by default; A6c sets 0.5C).
    Throughput: Δn = ΔSOC/2 (discharge) and c Δt/2 (charge); Δn_ch = c Δt.
    """
    seg: list[Seg] = []
    ztop = 0.8 if strategy == "limit80" else 1.0
    nd = 16
    T_use = Tamb_C + (0.0 if cal_mode else float(p.k_env * p.R_th) * P_day) + K0
    for i in range(nd):
        T = heat_C + K0 if (heat_h and i < heat_h) else T_use
        seg.append(Seg(1.0, ztop - dsoc * (i + 0.5) / nd, T, (dsoc / nd) / 2.0, 0.0, 0.0))
    z = max(ztop - dsoc, 0.0)
    t = 0.0
    Tch_C = (Tamb_C if Tch_amb_C is None else Tch_amb_C) + (0.0 if cal_mode else float(p.dT_ch))
    Tch = Tch_C + K0
    target = 0.8 if strategy in ("limit80", "obc") else 1.0
    cc, cv = (jeita_rate(c_cc, Tch_C), jeita_rate(c_cv, Tch_C)) if jeita else (c_cc, c_cv)

    def charge_to(z, t, ztarget, t_stop):
        while z < ztarget - 1e-9 and t < t_stop - 1e-12:
            c = cc if z < z_cc - 1e-9 else cv
            if c <= 0:
                break
            zstop = min(ztarget, z_cc) if z < z_cc - 1e-9 else ztarget
            dt = min(dt_charge_h, (zstop - z) / c, t_stop - t)
            seg.append(Seg(dt, z + c * dt / 2.0, Tch, c * dt / 2.0, c * dt, c))
            z += c * dt
            t += dt
        return z, t

    z, t = charge_to(z, t, target, 8.0)
    hold_end = 8.0 - (1.0 if strategy == "obc" else 0.0)
    if t < hold_end:
        seg.append(Seg(hold_end - t, z, Tamb_C + K0, 0.0, 0.0, 0.0))
        t = hold_end
    if strategy == "obc":
        cc_save = cc
        cc = cv  # the resumed top-up runs at the CV rate
        z, t = charge_to(z, t, 1.0, 8.0)
        cc = cc_save
        if t < 8.0:
            seg.append(Seg(8.0 - t, z, Tamb_C + K0, 0.0, 0.0, 0.0))
    return seg


def advance_day(state: AgingState, profile: list[Seg], p: Params = RP63, k_H=None, k_L=None) -> AgingState:
    """Integrate E17–E19 exactly over the piecewise-constant segments of one day."""
    kH = p.k_H if k_H is None else k_H
    kL = p.k_L if k_L is None else k_L
    QC, QH, QL = state
    for s in profile:
        QC += float(p.k_c * a_deg(s.T_K, p.E_cal) * g_stress(s.z, p)) ** 2 * s.dt_h / 24.0
        QH += float(kH * a_deg(s.T_K, p.E_H)) ** 2 * s.dn
        if s.dn_ch > 0:
            QL += float(kL * a_R(s.T_K, p.E_L) * math.exp(p.beta_L * (s.c_ch - p.c_ref))) ** 2 * s.dn_ch
    return AgingState(QC, QH, QL)


def calibrate(p: Params = RP63, N_80: float | None = None, soh_target: float = 0.80) -> tuple[float, float]:
    """Protocol C0 (§11.1): (k_H, k_L) such that SOH(N_80 FCE) = soh_target with plating share w_L."""
    c = calibrate_details(p, N_80, soh_target)
    return c["k_H"], c["k_L"]


def calibrate_details(p: Params = RP63, N_80: float | None = None, soh_target: float = 0.80) -> dict:
    """Protocol C0 (§11.1) with diagnostics.

    C0: cell at a constant 25 °C; each day discharge 100 → 0 % over 16 h, charge per E21, hold at 100 %
    (1 FCE/day).  Returns dict(k_H, k_L, q_cal_N, mean_g), where mean_g = q_cal/(k_c √N) is the
    root-mean-square of g(z) over the protocol (the quantity the spec calls "mean g").
    """
    N = int(round(p.N_80 if N_80 is None else N_80))
    prof = day_profile("full", 25.0, 1.0, p, cal_mode=True)
    st = AgingState()
    for _ in range(N):
        st = advance_day(st, prof, p, k_H=1.0, k_L=1.0)
    qc = st.q_cal
    need = (1.0 - soh_target) - qc
    kH = (1.0 - p.w_L) * need / math.sqrt(st.Q_H)
    kL = p.w_L * need / math.sqrt(st.Q_L)
    return dict(k_H=kH, k_L=kL, q_cal_N=qc, mean_g=qc / (float(p.k_c) * math.sqrt(N)))


AGING_SCENARIOS = {
    "A1": dict(strategy="full", Tamb_C=25.0),
    "A2": dict(strategy="limit80", Tamb_C=25.0),
    "A3": dict(strategy="obc", Tamb_C=25.0),
    "A4": dict(strategy="full", Tamb_C=35.0),
    "A5": dict(strategy="full", Tamb_C=25.0, heat_h=2.0),
    "A6a": dict(strategy="full", Tamb_C=25.0, c_cc=1.5),
    "A6b": dict(strategy="full", Tamb_C=25.0, Tch_amb_C=0.0),
    "A6c": dict(strategy="full", Tamb_C=25.0, Tch_amb_C=0.0, c_cc=0.5, c_cv=0.3),
    "A7": dict(strategy="full", Tamb_C=25.0, dsoc=0.95),
}


def run_years(strategy: str = "full", Tamb_C: float = 25.0, years: float = 3.0, p: Params = RP63,
              P_day: float = 0.440, dsoc: float | Callable | None = None, recalibrate: bool = False,
              **opts) -> dict:
    """Daily aging driver (§11.2) for ``years``; returns per-day arrays.

    ``dsoc``: fixed daily SOC swing, a callable of q_loss, or None for
    ΔSOC = min(0.95, 24 h·P_day/(E_n (1 − q_loss))).  ``recalibrate``: recompute k_H, k_L with protocol C0 for
    this ``p`` (needed when aging parameters are sampled).  Other options go to :func:`day_profile`.
    Returns dict(day, t_yr, soh, q_cal, q_H, q_L, r_ratio = R/R_new, years_to_80, k_H, k_L).
    """
    if recalibrate:
        kH, kL = calibrate(p)
        p = replace(p, k_H=kH, k_L=kL)
    En = float(cell.nominal_energy(p))
    if dsoc is None:
        dfn = lambda q: min(0.95, 24.0 * P_day / (En * (1.0 - q)))
    elif callable(dsoc):
        dfn = dsoc
    else:
        dfn = lambda q, d=float(dsoc): d
    ndays = int(round(years * 365))
    st = AgingState()
    out = np.zeros((ndays, 4))
    cache: dict[float, list[Seg]] = {}
    for d in range(ndays):
        ds = round(dfn(st.q_loss), 12)
        if ds not in cache:
            cache = {ds: day_profile(strategy, Tamb_C, ds, p, P_day=P_day, **opts)}
        st = advance_day(st, cache[ds], p)
        out[d] = (st.q_cal, st.q_H, st.q_L, st.q_loss)
    soh = 1.0 - out[:, 3]
    below = np.nonzero(soh <= 0.8)[0]
    return dict(day=np.arange(1, ndays + 1), t_yr=np.arange(1, ndays + 1) / 365.0, soh=soh, q_cal=out[:, 0],
                q_H=out[:, 1], q_L=out[:, 2], r_ratio=1.0 + float(p.gamma_R) * out[:, 3],
                years_to_80=(below[0] + 1) / 365.0 if below.size else math.inf, k_H=float(p.k_H),
                k_L=float(p.k_L))


def closed_form_losses(t_days=0.0, n_fce=0.0, n_ch=0.0, T_K=T_REF, z=0.5, c_ch=1.0, p: Params = RP63) -> dict:
    """Constant-condition closed forms (F11): q_cal = k_c a_deg g √t_d, q_H = k_H a_deg √n, q_L = ... √n_ch."""
    q_cal = p.k_c * a_deg(T_K, p.E_cal) * g_stress(z, p) * np.sqrt(t_days)
    q_H = p.k_H * a_deg(T_K, p.E_H) * np.sqrt(n_fce)
    q_L = p.k_L * a_R(T_K, p.E_L) * np.exp(p.beta_L * (c_ch - p.c_ref)) * np.sqrt(n_ch)
    return dict(q_cal=q_cal, q_H=q_H, q_L=q_L, q_loss=q_cal + q_H + q_L)


__all__ = ["AgingState", "Seg", "g_stress", "jeita_rate", "day_profile", "advance_day", "calibrate",
           "calibrate_details",
           "AGING_SCENARIOS", "run_years", "closed_form_losses"]
