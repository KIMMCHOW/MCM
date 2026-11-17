"""Component power sub-models at the battery terminals (FINAL_SPEC §2.1, §2.8, §3).

``power_parts(u, p)`` returns the throttling-independent power, the SoC demand, the burst amplitude and a
component breakdown.  Everything broadcasts over array-valued ``Params`` / ``UsageState`` fields.
"""
from __future__ import annotations

import numpy as np

from .cell import bisect
from .params import K0, Params
from .usage import PRIMITIVES, UsageState, wifi_cliff

BREAKDOWN_KEYS = ("susp", "awake", "wake", "disp", "net", "gnss", "cam", "aud", "sens", "extra")


def soc_demand(x_c, p: Params):
    """SoC demand P_SoC,dem(x_c) = P_SoC,max · x · [v0 + (1 − v0) x]² (E13) [W]."""
    x = np.asarray(x_c, dtype=float)
    return p.P_SoC_max * x * (p.v0 + (1.0 - p.v0) * x) ** 2


def m_sig(RSRP, p: Params):
    """Cellular weak-signal multiplier m_sig = min{m_max, 1 + γ_sig·max(0, (RSRP_ref − RSRP)/10 dB)} (E14)."""
    return np.minimum(p.m_max, 1.0 + p.gamma_sig * np.maximum(0.0, (p.RSRP_ref - np.asarray(RSRP, float)) / 10.0))


def radio_power(u: UsageState, p: Params):
    """P_net (E14): idle of the associated radios plus the traffic term Φ_ρ of the data path."""
    net = u.w * p.wifi_P_idle
    if u.c != "none":
        net = net + getattr(p, f"{u.c}_P_idle")
    if u.rho == "none":
        return net + 0.0 * p.P_susp
    rat = "wifi" if u.rho == "wifi" else u.c
    if rat == "none":
        raise ValueError(f"{u.name}: cellular data path without a camping RAT")
    r = p.radio(rat)
    if u.rho == "wifi":
        m = 1.0 if u.RSSI is None else wifi_cliff(u.RSSI)
    else:
        m = m_sig(u.RSRP, p)
    if u.kappa == "rt":
        phi = m * r["P_rt"]
    else:
        lam = u.lam
        phi = m * ((r["P_tail"] - r["P_idle"]) * (1.0 - np.exp(-lam * r["T_tail"]))
                   + (r["P_act"] - r["P_tail"]) * (1.0 - np.exp(-lam * r["T_act"]))
                   + lam * np.exp(-lam * r["T_tail"]) * r["E_pro"]
                   + r["e_b"] * (u.d_dl + p.k_ul * u.d_ul))
    return net + phi


def display_power(u: UsageState, p: Params):
    """OLED display power (E12)."""
    return u.s_scr * (p.P_d0 + p.p_c * p.A_disp * (np.asarray(u.b, float) / 0.5) ** p.gamma_b
                      * (u.m_APL * u.APL) / p.APL_ref + p.k_RR * np.maximum(np.asarray(u.f_r, float) - 60.0, 0.0))


def power_parts(u: UsageState, p: Params):
    """Battery-side power demand of input u (E1).

    Returns ``(P_nonsoc, P_soc_dem, dPb, breakdown)``: the throttling-independent power [W], the SoC demand
    before throttling [W], the burst-envelope amplitude ΔP_b·m_b [W] (0 when the envelope is disabled) and a
    dict of component powers (keys ``BREAKDOWN_KEYS`` plus ``soc_dem``).
    """
    a = np.maximum(u.s_scr, u.awake)
    br = {
        "susp": p.P_susp + 0.0 * a,
        "awake": p.P_awake * a,
        "wake": (1.0 - a) * u.lam * p.E_wake,
        "disp": display_power(u, p),
        "net": radio_power(u, p),
        "gnss": u.g * p.P_G,
        "cam": u.c_r * p.P_cam + u.c_f * p.P_camf,
        "aud": u.au * p.P_aud + u.bt * p.P_BT,
        "sens": u.P_sens + 0.0 * p.P_susp,
        "extra": u.P_extra + 0.0 * p.P_susp,
    }
    dem = soc_demand(u.x_c, p)
    br["soc_dem"] = dem
    Pns = sum(br[k] for k in BREAKDOWN_KEYS)
    dPb = u.dPb * p.m_b if p.burst else 0.0 * p.m_b
    return Pns, dem, dPb, br


def demand_power(u: UsageState, p: Params):
    """Total demand before throttling, P_nonsoc + P_SoC,dem [W] (the 'P [PROTO]' column of §5.1)."""
    Pns, dem, _, _ = power_parts(u, p)
    return Pns + dem


def throttle(P_dem, Ts, p: Params):
    """Throttled SoC power (E9); ``Ts`` is the skin temperature in K."""
    if not p.throttling:
        return P_dem + 0.0 * Ts
    sig = 1.0 / (1.0 + np.exp(-(p.T_lim - p.dT_thr - Ts) / p.w_thr))
    return np.where(P_dem > p.P_fl, p.P_fl + (P_dem - p.P_fl) * sig, P_dem)


def x_eff(P_soc, p: Params):
    """Delivered SoC throughput x_eff solving P_SoC,dem(x_eff) = P_SoC (§3.6), by bisection on [0, 1]."""
    P_soc = np.asarray(P_soc, dtype=float)
    return bisect(lambda x: soc_demand(x, p) - P_soc, np.zeros_like(P_soc + 0.0 * p.v0),
                  np.ones_like(P_soc + 0.0 * p.v0), n_iter=60)


def rho_perf(P_soc, x_c, p: Params):
    """Performance ratio ρ_perf = x_eff/x_c (1 when x_c = 0)."""
    x_c = np.asarray(x_c, dtype=float)
    return np.where(x_c > 0, x_eff(P_soc, p) / np.where(x_c > 0, x_c, 1.0), 1.0)


def skin_temperature(T, Tamb, p: Params):
    """Quasi-steady skin hot spot T_s = T_amb + (R_skin/R_th)(T − T_amb) (E8), K."""
    return Tamb + p.R_skin / p.R_th * (T - Tamb)


def steady_skin_temperature(u: UsageState, p: Params, Tamb_C=25.0, q_batt=0.0):
    """Thermal fixed point T_s* = T_amb + k_env R_skin [P(T_s*) + q_batt] (F5), K."""
    Pns, dem, _, _ = power_parts(u, p)
    Ta = np.asarray(Tamb_C, dtype=float) + K0

    def f(Ts):
        return Ts - Ta - p.k_env * p.R_skin * (Pns + throttle(dem, Ts, p) + q_batt)

    shape = np.broadcast(np.asarray(Pns), np.asarray(dem), Ta, np.asarray(p.R_skin), np.asarray(p.k_env)).shape
    lo = np.broadcast_to(Ta - 1.0, shape)
    hi = np.broadcast_to(Ta + 200.0, shape)
    return bisect(f, lo, hi, n_iter=60)


def steady_power(u: UsageState | str, p: Params, Tamb_C=25.0, q_batt=0.0):
    """Quasi-steady (throttled) battery power P* at the F5 fixed point [W] (q_batt neglected by default)."""
    u = PRIMITIVES[u] if isinstance(u, str) else u
    Pns, dem, _, _ = power_parts(u, p)
    Ts = steady_skin_temperature(u, p, Tamb_C, q_batt)
    return Pns + throttle(dem, Ts, p)


__all__ = [
    "PRIMITIVES", "BREAKDOWN_KEYS", "soc_demand", "m_sig", "radio_power", "display_power", "power_parts",
    "demand_power", "throttle", "x_eff", "rho_perf", "skin_temperature", "steady_skin_temperature",
    "steady_power",
]
