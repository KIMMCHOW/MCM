"""Cell and pack electrical model (FINAL_SPEC §2.2, §2.5, §2.6, §2.9, §2.11).

All functions broadcast: ``z``/``T`` may be scalars or arrays, and any ``Params`` field may be an array
over a trailing sample axis.  Temperatures are in kelvin.
"""
from __future__ import annotations

import numpy as np

from .params import RG, T_REF, Params

# ∂U/∂T table (V/K): docs/research/data/cell_ecm_pybamm/pybamm_ocv_fit.json, key dUdT_cell_mV_per_K (§2.5).
DUDT_Z = np.round(np.arange(0.0, 1.0001, 0.05), 2)
DUDT_V = 1e-3 * np.array([-0.4485, -0.7057, -0.5343, -0.4982, -0.4787, -0.4213, -0.3722, -0.3439, -0.3068,
                          -0.2655, -0.2482, -0.2412, -0.2482, -0.2281, -0.2047, -0.1798, -0.1528, -0.1246,
                          -0.0966, -0.0705, -0.0459])


# ------------------------------------------------------------------------------------------- Arrhenius
def a_R(T, E):
    """Resistance-like Arrhenius factor exp[+(E/R_g)(1/T − 1/T_ref)]: > 1 when colder than 25 °C."""
    return np.exp(E / RG * (1.0 / T - 1.0 / T_REF))


def a_deg(T, E):
    """Degradation-like Arrhenius factor exp[−(E/R_g)(1/T − 1/T_ref)]: > 1 when hotter than 25 °C."""
    return np.exp(-E / RG * (1.0 / T - 1.0 / T_REF))


# ------------------------------------------------------------------------------------------- OCV (E5)
def dUdT(z, p: Params | None = None):
    """Entropic coefficient ∂U/∂T(z) [V/K]: linear interpolation in the 21-point table, clamped to [0, 1]."""
    zt, vt = (DUDT_Z, DUDT_V) if p is None or p.dudt_table is None else p.dudt_table
    return np.interp(np.minimum(np.maximum(z, 0.0), 1.0), zt, vt)


def ocv_ref(z, p: Params, clamped: bool = False):
    """OCV at T_ref, U(z, T_ref) [V] (E5 without the entropic term).  z is clamped to [0, 1]
    (pass ``clamped=True`` if it already is)."""
    zc = z if clamped else np.minimum(np.maximum(z, 0.0), 1.0)
    if p.ocv_table is not None:
        zt, ut = p.ocv_table
        return np.interp(zc, zt, ut) + p.dU
    return (-p.b0 * np.exp(-p.b1 * zc) + p.a0 + zc * (p.a1 + zc * (p.a2 + zc * p.a3))
            + p.dV_hv * zc ** p.p_hv + p.dU)


def ocv(z, T, p: Params):
    """Open-circuit voltage U(z, T) [V] (E5), including (T − T_ref)·k_S·∂U/∂T(z)."""
    return ocv_ref(z, p) + (T - T_REF) * p.k_S * dUdT(z, p)


def docv_dz(z, p: Params):
    """dU/dz at T_ref [V per unit SOC] (analytic for the CRM06 form, central difference for a table)."""
    if p.ocv_table is not None:
        dz = 1e-5
        return (ocv_ref(z + dz, p) - ocv_ref(z - dz, p)) / (2 * dz)
    return (p.b0 * p.b1 * np.exp(-p.b1 * z) + p.a1 + 2 * p.a2 * z + 3 * p.a3 * z * z
            + p.dV_hv * p.p_hv * z ** (p.p_hv - 1.0))


def ocv_integral(z, p: Params):
    """Antiderivative F(z) = ∫_0^z U(ζ, T_ref) dζ [V] (closed form for E5; trapezoid for a table)."""
    zc = np.clip(z, 0.0, 1.0)
    if p.ocv_table is not None:
        zt, ut = (np.asarray(a, dtype=float) for a in p.ocv_table)
        zf = np.linspace(0.0, 1.0, 4001)
        uf = np.interp(zf, zt, ut)
        F = np.concatenate([[0.0], np.cumsum(0.5 * (uf[1:] + uf[:-1]) * np.diff(zf))])
        return np.interp(zc, zf, F) + p.dU * zc
    return (p.b0 / p.b1 * (np.exp(-p.b1 * zc) - 1.0) + p.a0 * zc + p.a1 * zc ** 2 / 2 + p.a2 * zc ** 3 / 3
            + p.a3 * zc ** 4 / 4 + p.dV_hv * zc ** (p.p_hv + 1.0) / (p.p_hv + 1.0) + p.dU * zc)


# ------------------------------------------------------------------------------------------- resistances (E6)
def r_growth(p: Params):
    """Resistance growth fraction (E20): γ_R·q_loss unless set explicitly."""
    return p.rg


def R0(z, T, p: Params):
    """Pack ohmic resistance R0(z, T) [Ω] (E6a); R0(0.5, T_ref) = R0_ref (fresh)."""
    rinf = p.R0_ref / (1.0 + p.A_R * np.exp(-0.5 / p.z0))
    return rinf * (1.0 + p.A_R * np.exp(-np.maximum(z, 0.0) / p.z0)) * a_R(T, p.Ea_R0) * (1.0 + p.rg)


def R1(T, p: Params):
    """RC-branch resistance R1(T) [Ω] (E6b)."""
    return p.R1_ref * a_R(T, p.Ea_R1) * (1.0 + p.rg)


def tau1(T, p: Params):
    """RC time constant τ1(T) [s] (E6c)."""
    return p.tau1_ref * a_R(T, p.Ea_tau)


def R_sigma(z, T, p: Params):
    """Quasi-static series resistance R0 + R1 [Ω]."""
    return R0(z, T, p) + R1(T, p)


# ------------------------------------------------------------------------------------------- constant-power root
def current(Ue, R0_, P):
    """Constant-power operating point (E2).

    Returns ``(I, V, D)`` with D = Ũ² − 4 R0 P.  For D > 0, I = 2P/(Ũ + √D) (cancellation-free physical root)
    and V = (Ũ + √D)/2.  For D ≤ 0 the current is clamped to the maximum-power current Ũ/(2 R0) and V = Ũ/2
    (the 'infeasible' event fires there).
    """
    D = Ue * Ue - 4.0 * R0_ * P
    sq = np.sqrt(np.maximum(D, 0.0))
    ok = D > 0
    I = np.where(ok, 2.0 * P / np.where(ok, Ue + sq, 1.0), Ue / (2.0 * R0_))
    V = Ue - I * R0_
    return I, V, D


def peak_voltage(Ue, R0_, Ppk):
    """Terminal voltage under a sub-second burst P_pk with v1 frozen (E10); 0 if D_pk ≤ 0."""
    Dk = Ue * Ue - 4.0 * R0_ * Ppk
    return np.where(Dk > 0, 0.5 * (Ue + np.sqrt(np.maximum(Dk, 0.0))), 0.0)


def p_max(Ue, R0_):
    """Largest deliverable power Ũ²/(4 R0) [W]."""
    return Ue * Ue / (4.0 * R0_)


def p_cut(Ue, R, Vstar):
    """Power at which the terminal voltage equals V*: P_cut(V*) = V*(Ũ − V*)/R [W] (thermal_aging T3)."""
    return Vstar * (Ue - Vstar) / R


# ------------------------------------------------------------------------------------------- roots / SOC maps
def bisect(f, lo, hi, n_iter: int = 64):
    """Vectorised bisection for an increasing function f on [lo, hi] (arrays broadcast)."""
    lo = np.asarray(lo, dtype=float) + 0.0
    hi = np.asarray(hi, dtype=float) + 0.0
    lo, hi = np.broadcast_arrays(lo, hi)
    lo, hi = lo.copy(), hi.copy()
    for _ in range(n_iter):
        mid = 0.5 * (lo + hi)
        pos = f(mid) > 0
        hi = np.where(pos, mid, hi)
        lo = np.where(pos, lo, mid)
    return 0.5 * (lo + hi)


def z_empty(p: Params, T=T_REF):
    """z_e = U⁻¹(V_cut; T): the coulombic SOC at which the open-circuit voltage equals the cut-off (§2.11)."""
    shape = np.broadcast(np.asarray(p.V_cut), np.asarray(p.dU), np.asarray(p.Q_n), np.asarray(T)).shape
    return bisect(lambda z: ocv(z, T, p) - p.V_cut, np.full(shape, 1e-9), np.full(shape, 0.6), n_iter=60)


def z_from_s0(s0, p: Params):
    """Initial coulombic SOC from the displayed SOC s0 (D1, E11): z0 = z_e + s0 (1 − z_e)."""
    ze = z_empty(p)
    return ze + np.asarray(s0, dtype=float) * (1.0 - ze)


def usable_energy(p: Params, z0=1.0):
    """E_u(z0) = Q_max ∫_{z_e}^{z0} U(z, T_ref) dz [Wh] (F1)."""
    ze = z_empty(p)
    return p.Q_max * (ocv_integral(z0, p) - ocv_integral(ze, p))


def nominal_energy(p: Params):
    """E_n = Q_n ∫_0^1 U(z, T_ref) dz [Wh]."""
    return p.Q_n * (ocv_integral(1.0, p) - ocv_integral(0.0, p))


def mean_ocv(p: Params):
    """Mean OCV over 0 ≤ z ≤ 1 at T_ref [V]."""
    return ocv_integral(1.0, p) - ocv_integral(0.0, p)


def display_soc(z, p: Params, mode: str = "D1"):
    """Displayed SOC from coulombic z (E11).

    D1: s = (z − z_e)/(1 − z_e) (coulombic, used for s0 and s_off).
    D2: ς = ∫_{z_e}^{z} U dz / ∫_{z_e}^{1} U dz (energy fraction of the problem glossary).
    """
    ze = z_empty(p)
    if mode == "D1":
        return (z - ze) / (1.0 - ze)
    if mode == "D2":
        Fe = ocv_integral(ze, p)
        return (ocv_integral(z, p) - Fe) / (ocv_integral(1.0, p) - Fe)
    raise ValueError(mode)


def soc_to_z(s, p: Params, mode: str = "D1"):
    """Inverse of :func:`display_soc` (D2 by bisection)."""
    ze = z_empty(p)
    if mode == "D1":
        return ze + np.asarray(s, dtype=float) * (1.0 - ze)
    if mode == "D2":
        s = np.asarray(s, dtype=float)
        return bisect(lambda z: display_soc(z, p, "D2") - s, np.zeros_like(s + ze), np.ones_like(s + ze))
    raise ValueError(mode)


def ios_display(s, h: float = 0.065):
    """iOS top-hold display map (Bonn comparison only): floor_1%(min(1, s/(1 − h)))."""
    return np.floor(100.0 * np.minimum(1.0, np.asarray(s) / (1.0 - h)) + 1e-9) / 100.0


__all__ = [
    "DUDT_Z", "DUDT_V", "a_R", "a_deg", "dUdT", "ocv_ref", "ocv", "docv_dz", "ocv_integral", "r_growth", "R0",
    "R1", "tau1", "R_sigma", "current", "peak_voltage", "p_max", "p_cut", "bisect", "z_empty", "z_from_s0",
    "usable_energy", "nominal_energy", "mean_ocv", "display_soc", "soc_to_z", "ios_display",
]
