"""Analytic and approximate results, the insight layer (FINAL_SPEC §7, F1–F11, Appendix A).

All functions broadcast over array-valued powers and ``Params`` fields (trailing sample axis), except the
matrix routines of the CTMC (one generator at a time; loop over θ samples, milliseconds each).
Units: power W, energy W·h, time h, generator rates 1/h.
"""
from __future__ import annotations

import math
from typing import Callable, Mapping, Sequence

import numpy as np

from . import cell
from .cell import R0, R1, bisect, ocv
from .params import K0, RP63, T_REF, Params
from .power import power_parts, steady_power, steady_skin_temperature, throttle, x_eff
from .usage import CTMC, LEVERS, PRIMITIVES, Schedule, UsageState, s12_schedule

_GL_CACHE: dict[int, tuple[np.ndarray, np.ndarray]] = {}


def _gauss_legendre(n: int):
    if n not in _GL_CACHE:
        _GL_CACHE[n] = np.polynomial.legendre.leggauss(n)
    return _GL_CACHE[n]


def _integrate(f: Callable, a, b, n: int = 200):
    """∫_a^b f(z) dz by n-point Gauss–Legendre, vectorised over the broadcast shape of a and b."""
    a, b = np.broadcast_arrays(np.asarray(a, float), np.asarray(b, float))
    x, w = _gauss_legendre(n)
    shp = (n,) + (1,) * a.ndim
    mid, half = 0.5 * (a + b), 0.5 * (b - a)
    z = mid + half * x.reshape(shp)
    return half * np.sum(w.reshape(shp) * f(z), axis=0)


# ------------------------------------------------------------------------------------------- F1 energy bucket
def tte_L0(P, p: Params = RP63, s0=1.0):
    """L0 time-to-empty E_u(z0)/P [h] (F1), z0 from the displayed s0."""
    return cell.usable_energy(p, cell.z_from_s0(s0, p)) / np.asarray(P, dtype=float)


def energy_shares(parts: Mapping[str, np.ndarray]) -> dict:
    """Energy shares s_i of components (= −elasticity of the L0 TTE w.r.t. P_i, F1)."""
    tot = sum(np.asarray(v, float) for v in parts.values())
    return {k: np.asarray(v, float) / tot for k, v in parts.items()}


# ------------------------------------------------------------------------------------------- F2 L1 quadrature
def v_qs(z, P, T, p: Params):
    """Quasi-static terminal voltage [U + √(U² − 4 R_Σ P)]/2 with the RC branch equilibrated."""
    U = ocv(z, T, p)
    Rs = R0(z, T, p) + R1(T, p)
    return 0.5 * (U + np.sqrt(np.maximum(U * U - 4.0 * Rs * P, 0.0)))


def cutoff_soc(P, T, p: Params, z0=1.0):
    """z_c solving U(z_c, T) = V_cut + P R_Σ(z_c, T)/V_cut, i.e. V_qs(z_c) = V_cut (bisection on [1e-6, z0])."""
    P = np.asarray(P, dtype=float)
    shape = np.broadcast(P, np.asarray(T), np.asarray(z0), np.asarray(p.V_cut), np.asarray(p.R0_ref)).shape
    return bisect(lambda z: v_qs(z, P, T, p) - p.V_cut, np.full(shape, 1e-6), np.broadcast_to(z0, shape) + 0.0,
                  n_iter=60)


def tte_L1(P, p: Params = RP63, T_C=25.0, s0=1.0, n_gl: int = 200, return_zc: bool = False):
    """Exact quasi-static TTE for constant P and T (F2): (Q_max/P) ∫_{z_c}^{z0} V_qs dz [h].

    Vectorised (Gauss–Legendre on [z_c, z0], bisection for z_c), so it serves as the Sobol surrogate (§10.3).
    """
    T = np.asarray(T_C, dtype=float) + K0
    P = np.asarray(P, dtype=float)
    z0 = cell.z_from_s0(s0, p)
    zc = cutoff_soc(P, T, p, z0)
    I = _integrate(lambda z: v_qs(z, P, T, p), zc, z0, n=n_gl)
    tte = p.Q_max / P * I
    return (tte, zc) if return_zc else tte


# ------------------------------------------------------------------------------------------- F3 resistive constant
def a_coefficient(p: Params = RP63, T_C=25.0, n_gl: int = 200) -> dict:
    """First-order resistive constant a of E_u(P) ≈ E_u0 (1 − a P) (F3, Appendix A.4) [1/W].

    a = [∫_{z_e}^1 R_Σ/U dz + R_Σ(z_e)/U′(z_e)] / ∫_{z_e}^1 U dz, with U, R_Σ at temperature T (z_e at T_ref).
    """
    T = np.asarray(T_C, dtype=float) + K0
    ze = cell.z_empty(p)
    EU = _integrate(lambda z: ocv(z, T, p), ze, 1.0, n_gl)
    ohm = _integrate(lambda z: (R0(z, T, p) + R1(T, p)) / ocv(z, T, p), ze, 1.0, n_gl) / EU
    dU = cell.docv_dz(ze, p)
    cut = (R0(ze, T, p) + R1(T, p)) / dU / EU
    return dict(ohm=ohm, cut=cut, total=ohm + cut, dUdz_ze=dU)


def elasticity_L1(P, a):
    """∂ln TTE/∂ln P = −1 − aP/(1 − aP) (F3)."""
    aP = np.asarray(a) * np.asarray(P)
    return -1.0 - aP / (1.0 - aP)


def jensen_L1(mean_L0, a, Pbar, cv_P):
    """Burstiness penalty E[TTE]_L1 ≈ E[TTE]_L0 (1 − a P̄ (1 + CV_P²)) (F10)."""
    return np.asarray(mean_L0) * (1.0 - np.asarray(a) * np.asarray(Pbar) * (1.0 + np.asarray(cv_P) ** 2))


# ------------------------------------------------------------------------------------------- F4 brown-out
def brownout_root(P, dPb, T_K, p: Params = RP63, v1=None) -> dict:
    """Brown-out threshold z* solving U(z*,T) = V_sys + v1 + R0(z*,T)(P + ΔP_b)/V_sys (F4, A.5).

    ``dPb`` is the burst amplitude in W (already including m_b); ``v1`` defaults to R1 P/U(z*) (RC at
    equilibrium under the mean load).  Returns dict(z_star, s_off) with s_off = (z* − z_e)/(1 − z_e).
    """
    P = np.asarray(P, dtype=float)
    Ppk = P + np.asarray(dPb, dtype=float)
    Vs = p.V_sys

    def f(z):
        U = ocv(z, T_K, p)
        vv = R1(T_K, p) * P / U if v1 is None else v1
        return U - Vs - vv - R0(z, T_K, p) * Ppk / Vs

    shape = np.broadcast(P, Ppk, np.asarray(T_K), np.asarray(p.R0_ref), np.asarray(Vs)).shape
    z = bisect(f, np.zeros(shape), np.ones(shape), n_iter=60)
    ze = cell.z_empty(p)
    return dict(z_star=z, s_off=(z - ze) / (1.0 - ze))


def safe_burst(z, T_K, P, p: Params = RP63, v1=None):
    """Largest burst above the mean load that avoids a brown-out: V_sys(Ũ − V_sys)/R0 − P [W] (F4)."""
    U = ocv(z, T_K, p)
    vv = R1(T_K, p) * np.asarray(P) / U if v1 is None else v1
    Ue = U - vv
    return p.V_sys * (Ue - p.V_sys) / R0(z, T_K, p) - np.asarray(P)


def pcut_table(p: Params = RP63, T_C=25.0, zs=(0.5, 0.2, 0.1)) -> dict:
    """Cut-off-limited power P_cut(V_cut, R_Σ) and burst-limited P_cut,sys(V_sys, R0) at given z [W] (F4)."""
    T = float(T_C) + K0
    z = np.asarray(zs, dtype=float)
    U = ocv(z, T, p)
    return dict(z=z, P_cut=cell.p_cut(U, R0(z, T, p) + R1(T, p), p.V_cut),
                P_cut_sys=cell.p_cut(U, R0(z, T, p), p.V_sys))


# ------------------------------------------------------------------------------------------- F5 thermal cap
def thermal_cap(u: UsageState | str = "GM", p: Params = RP63, Tamb_C=25.0, q_batt=0.0) -> dict:
    """Thermal fixed point of a constant activity (F5).

    Returns P_star (W), Ts_star_C, rho_perf, tte_cap_h = E_u/P* (L0 form, as used for the [PROTO] numbers),
    tte_cap_L1_h (F2 at P* and the steady battery temperature T_amb + k_env R_th P*), sharp-cap bound (W),
    dTTE/dT_amb (h/K, valid when the cap binds) and the non-SoC power.
    """
    u = PRIMITIVES[u] if isinstance(u, str) else u
    Pns, dem, _, _ = power_parts(u, p)
    Ts = steady_skin_temperature(u, p, Tamb_C, q_batt)
    Psoc = throttle(dem, Ts, p)
    Pstar = Pns + Psoc
    Eu = cell.usable_energy(p)
    xc = np.asarray(u.x_c, dtype=float)
    rho = np.where(xc > 0, x_eff(Psoc, p) / np.where(xc > 0, xc, 1.0), 1.0)
    sharp = np.maximum(Pns + p.P_fl, (p.T_lim_C - p.dT_thr - np.asarray(Tamb_C)) / (p.k_env * p.R_skin))
    T_batt = np.asarray(Tamb_C) + p.k_env * p.R_th * Pstar
    return dict(P_star=Pstar, Ts_star_C=Ts - K0, rho_perf=rho, tte_cap_h=Eu / Pstar,
                tte_cap_L1_h=tte_L1(Pstar, p, T_batt), sharp_cap_W=sharp,
                dTTE_dTamb=Eu / (p.k_env * p.R_skin * Pstar ** 2), P_nonsoc=Pns, P_demand=Pns + dem)


def thermal_step(t_s, T0, Tamb, q, p: Params = RP63):
    """Lumped thermal step response T(t) = T_amb + R q + (T0 − T_amb − R q) e^(−t/τ_th) (F5, HE19 Eq. 9)."""
    R = p.k_env * p.R_th
    return Tamb + R * q + (T0 - Tamb - R * q) * np.exp(-np.asarray(t_s) / (R * p.C_th))


# ------------------------------------------------------------------------------------------- F6–F7 CTMC moments
def stationary(L: np.ndarray) -> np.ndarray:
    """Stationary distribution π of a generator (πᵀΛ = 0, Σπ = 1)."""
    n = L.shape[0]
    A = np.vstack([L.T, np.ones(n)])
    b = np.zeros(n + 1)
    b[-1] = 1.0
    return np.linalg.lstsq(A, b, rcond=None)[0]


def ctmc_moments(L: np.ndarray, P: np.ndarray, Eu_Wh: float, x0: int = 0) -> dict:
    """PDMP mean (Dynkin, F6), CLT variance and inverse-Gaussian shape (F7) of the L0 TTE under a CTMC.

    Returns dict(pi, Pbar, sigma2 [W²h], h [Wh], Dm [h], mu_fluid, mu_dynkin, sd_clt [h], var_P, cv_P, tau_c [h],
    cv_tte (= sd_clt/mu_dynkin), lam_IG [h]).
    """
    L = np.asarray(L, dtype=float)
    P = np.asarray(P, dtype=float)
    pi = stationary(L)
    n = pi.size
    Pi = np.outer(np.ones(n), pi)
    Dm = np.linalg.inv(Pi - L) - Pi  # deviation matrix
    Pbar = float(pi @ P)
    hv = Dm @ (P - Pbar)
    s2 = float(2.0 * np.sum(pi * (P - Pbar) * hv))
    hP = float(np.sum(pi * P * hv) / Pbar)  # E[h(X_TTE)] ≈ Σ π_x P_x h_x / P̄ (inspection paradox)
    mu_f = Eu_Wh / Pbar
    mu_d = (Eu_Wh - hv[x0] + hP) / Pbar
    sd = math.sqrt(s2 * Eu_Wh / Pbar ** 3)
    varP = float(np.sum(pi * (P - Pbar) ** 2))
    return dict(pi=pi, Pbar=Pbar, sigma2=s2, h=hv, Dm=Dm, mu_fluid=mu_f, mu_dynkin=mu_d, sd_clt=sd, var_P=varP,
                cv_P=math.sqrt(varP) / Pbar, tau_c=s2 / (2.0 * varP) if varP > 0 else math.nan, cv_tte=sd / mu_d,
                lam_IG=Eu_Wh ** 2 / s2 if s2 > 0 else math.inf)


def two_state_sigma2(alpha, beta, P_on, P_off):
    """CLT variance rate of an on/off chain: 2αβ(P_on − P_off)²/(α + β)³ [W²h] (F7 check)."""
    return 2.0 * alpha * beta * (P_on - P_off) ** 2 / (alpha + beta) ** 3


def ig_quantiles(mu, lam, q=(0.05, 0.5, 0.95)):
    """Quantiles of IG(mean μ, shape λ) (F7)."""
    from scipy.stats import invgauss
    return invgauss.ppf(np.asarray(q, dtype=float), mu / lam, scale=lam)


def ig_cdf(t, mu, lam):
    """Pr(TTE < t) under IG(mean μ, shape λ)."""
    from scipy.stats import invgauss
    return invgauss.cdf(np.asarray(t, dtype=float), mu / lam, scale=lam)


def ig_pdf(t, mu, lam):
    from scipy.stats import invgauss
    return invgauss.pdf(np.asarray(t, dtype=float), mu / lam, scale=lam)


# ------------------------------------------------------------------------------------------- F8–F9 hub and spoke
def hub_closed_form(lam_on: float, p_j: Sequence[float], d_j_h: Sequence[float], P_idle: float,
                    P_j: Sequence[float], names: Sequence[str] | None = None) -> dict:
    """Hub-and-spoke closed forms (F8, Appendix A.7) for the activity chain alone.

    ``d_j_h`` are mean session durations in hours.  Returns E[L], P̄, σ², π, variance shares and the usage
    elasticities of μ_fluid with respect to d_j and λ_on.
    """
    p_j, d, Pj = (np.asarray(v, dtype=float) for v in (p_j, d_j_h, P_j))
    EL = 1.0 / lam_on + np.sum(p_j * d)
    Pbar = (P_idle / lam_on + np.sum(p_j * d * Pj)) / EL
    contrib = 2.0 * p_j * d ** 2 * (Pj - Pbar) ** 2 / EL
    s2 = float(np.sum(contrib))
    pi_j = p_j * d / EL
    pi_idle = (1.0 / lam_on) / EL
    P_sess = np.sum(p_j * d * Pj) / np.sum(p_j * d)
    el_d = -pi_j * (Pj - Pbar) / Pbar
    el_lam = -pi_idle * (1.0 - pi_idle) * (P_sess - P_idle) / Pbar
    names = list(names) if names is not None else [str(i) for i in range(len(p_j))]
    return dict(EL=EL, Pbar=Pbar, sigma2=s2, pi_idle=pi_idle, pi=dict(zip(names, pi_j)),
                share=dict(zip(names, contrib / s2)), elast_d=dict(zip(names, el_d)), elast_lam_on=el_lam,
                P_session=P_sess)


def semi_markov_sigma2(lam_on: float, p_j, d_j_h, P_idle: float, P_j, cv2_j=1.0, cv2_0=1.0) -> float:
    """CLT variance rate with general sojourn laws (F9): squared CVs ``cv2_j`` (sessions) and ``cv2_0`` (idle)."""
    p_j, d, Pj = (np.asarray(v, dtype=float) for v in (p_j, d_j_h, P_j))
    cv2_j = np.broadcast_to(np.asarray(cv2_j, dtype=float), p_j.shape)
    EL = 1.0 / lam_on + np.sum(p_j * d)
    Pbar = (P_idle / lam_on + np.sum(p_j * d * Pj)) / EL
    num = np.sum(p_j * (1.0 + cv2_j) * d ** 2 * (Pj - Pbar) ** 2) + (P_idle - Pbar) ** 2 * (cv2_0 - 1.0) / lam_on ** 2
    return float(num / EL)


# ------------------------------------------------------------------------------------------- Monte Carlo L0
def mc_L0(L: np.ndarray, P: np.ndarray, Eu_Wh: float, n: int = 20000, seed: int = 1, x0: int = 0) -> np.ndarray:
    """Exact L0 TTE samples [h] under the CTMC: energy E_u depleted at rate P_X(t) (vectorised over paths)."""
    rng = np.random.default_rng(seed)
    rates = -np.diag(L)
    J = (L - np.diag(np.diag(L))) / rates[:, None]
    cum = np.cumsum(J, axis=1)
    cum[:, -1] = 1.0
    P = np.asarray(P, dtype=float)
    s = np.full(n, x0, dtype=np.int64)
    E = np.full(n, float(Eu_Wh))
    t = np.zeros(n)
    out = np.full(n, np.nan)
    alive = np.ones(n, bool)
    while alive.any():
        idx = np.nonzero(alive)[0]
        si = s[idx]
        tau = rng.exponential(1.0, idx.size) / rates[si]
        e = P[si] * tau
        done = e >= E[idx]
        dI = idx[done]
        out[dI] = t[dI] + E[dI] / P[s[dI]]
        alive[dI] = False
        cI = idx[~done]
        E[cI] -= e[~done]
        t[cI] += tau[~done]
        u = rng.random(cI.size)
        s[cI] = (u[:, None] > cum[s[cI]]).sum(axis=1)
    return out


# ------------------------------------------------------------------------------------------- lever metrics (§12.1)
def schedule_mean_power(schedule: Schedule, p: Params = RP63, Tamb_C=25.0, steady: bool = True,
                        map_fn: Callable | None = None) -> float:
    """Time-averaged power of a schedule with quasi-steady state powers (F5 fixed points) [W]."""
    E = 0.0
    for d, u in schedule:
        uu = map_fn(u.name, u) if map_fn is not None else u
        E += d * float(steady_power(uu, p, Tamb_C) if steady else np.sum(power_parts(uu, p)[:2]))
    return E / schedule.period


def ctmc_mean_tte(ctmc: CTMC, p: Params = RP63, Tamb_C=25.0, map_fn: Callable | None = None, x0: int = 0):
    """Dynkin mean TTE [h] of a CTMC with quasi-steady state powers (F6)."""
    P = ctmc.state_powers(p, Tamb_C, capped=True, map_fn=map_fn)
    return ctmc_moments(ctmc.generator(), P, float(cell.usable_energy(p)), x0)["mu_dynkin"]


def lever_effects(p: Params = RP63, Tamb_C=25.0, levers: Mapping | None = None) -> dict:
    """Smooth lever metrics of §12.1: % gain of E_u/P̄ (S12 day) and of the S13 Dynkin mean.

    Also includes ``wifi_share_80``: connectivity sojourns 4 h Wi-Fi / 1 h LTE (CTMC only).
    """
    levers = LEVERS if levers is None else levers
    Eu = float(cell.usable_energy(p))
    sch = s12_schedule()
    ct = CTMC("typical")
    P0 = schedule_mean_power(sch, p, Tamb_C)
    mu0 = ctmc_mean_tte(ct, p, Tamb_C)
    out = {"base": dict(s12_Pbar_W=P0, s12_tte_inf_h=Eu / P0, ctmc_mu_h=mu0)}
    for k, (label, fn) in levers.items():
        Pb = schedule_mean_power(sch, p, Tamb_C, map_fn=fn)
        mu = ctmc_mean_tte(ct, p, Tamb_C, map_fn=fn)
        out[k] = dict(label=label, s12_pct=100 * (P0 / Pb - 1), ctmc_pct=100 * (mu / mu0 - 1))
    mu80 = ctmc_mean_tte(CTMC("typical", mean_w=4.0, mean_c=1.0), p, Tamb_C)
    out["wifi_share_80"] = dict(label="Prefer Wi-Fi: Wi-Fi share 50 -> 80 %", s12_pct=math.nan,
                                ctmc_pct=100 * (mu80 / mu0 - 1), ctmc_mu_h=mu80)
    return out


__all__ = [
    "tte_L0", "energy_shares", "v_qs", "cutoff_soc", "tte_L1", "a_coefficient", "elasticity_L1", "jensen_L1",
    "brownout_root", "safe_burst", "pcut_table", "thermal_cap", "thermal_step", "stationary", "ctmc_moments",
    "two_state_sigma2", "ig_quantiles", "ig_cdf", "ig_pdf", "hub_closed_form", "semi_markov_sigma2", "mc_L0",
    "schedule_mean_power", "ctmc_mean_tte", "lever_effects",
]
