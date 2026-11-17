"""AN9: global variance-based sensitivity analysis (Sobol), FINAL_SPEC §10.3.

Design
------
* Factors: the 24 groups of §10.3 (``phonebatt.params.SOBOL_GROUPS``); every parameter of a group moves jointly,
  each by its own §4 distribution (``UQ_SPEC``).  Unit-cube points come from SALib's Saltelli design
  (``SALib.sample.sobol.sample``, scrambled Sobol', seed 2026, second order: N(2k + 2) rows) and are mapped
  through the §4 inverse CDFs (``phonebatt.params.map_unit``).
* Estimators (SALib's pairing): S_i by Saltelli et al. (2010), S_i = E[f(B)(f(A_B^i) − f(A))]/V, and S_Ti by
  Jansen (1999), S_Ti = E[(f(A) − f(A_B^i))²]/(2V), V = Var(f(A) ∪ f(B)); 95 % percentile-bootstrap CIs over the
  N base rows (1000 resamples); ``*_conf`` is the half-width of that interval.  Jansen's first-order estimator
  is stored as ``S1_jansen`` (much noisier for small S_i).  Second-order S_ij (interaction study only) come from
  SALib's ``analyze`` (Saltelli 2002 estimator, its own 1000-resample bootstrap).
* Constant-load scenarios S3, S4, S5a/b/c, S6, S8a/b, S11a/c (25 °C, new, s0 = 100 %): vectorised surrogate,
  N = 1024 (k = 24: 51,200 evaluations each).  Surrogate = F2 (exact quasi-static quadrature, §7) at the thermal
  steady state T = T_amb + k_env R_th P*, with P* the F5 fixed point (throttling at the steady skin temperature)
  and the end point max(z_c, z*) where z* is the F4 brown-out root under P* + ΔP_b (v1 = R1 P/U).  Checked
  against the full L2 model (``run_batch``, rk4_batch, h = 2 s) on 200 random design points per scenario.
  S11c (video call at −110 dBm, up to 10 W) fails that check at the thermal-shutdown boundary (TTE jumps where the
  battery just reaches 68 °C), so its reported indices come from full L2 (N = 256) and the surrogate indices are
  kept as ``S11c_surrogate``.  The surrogate includes the battery heat q_batt and the thermal event (F5).
* Full L2 (``run_batch``, h = 2 s) for S7 (gaming, 25 °C, new), S9 (cold hike, −10 °C, SOH 0.8, T0 = 20 °C) and
  S10a (social, −10 °C, SOH 0.8): N = 256 (12,800 runs each).  Outputs TTE and the stranded charge s_off.
* Interaction study S10a-G (social S4): T_amb ~ U(−20, 10) °C, SOH ~ U(0.75, 1.0) and m_b ~ U(0.5, 2.0) are
  explicit factors next to the other 23 groups (k = 26), N = 512 (27,648 runs), T0 = T_amb.  Reports S_i, S_Ti,
  S_ij within {T_amb, SOH, m_b} and the interaction share 1 − Σ S_i, for TTE and for s_off.

Writes results/sobol.json (keys per FINAL_SPEC §14: ``{S}.{S1, ST, S1_conf, ST_conf}`` (dicts keyed by factor
group), ``S10a-G.{S2, interaction_share}``, ``surrogate_check_maxerr_pct``; extras: ``S1_ci``/``ST_ci``,
``S1_jansen``, ``s_off`` blocks for the full-L2 scenarios, ``surrogate_check``) and paper/figures/fig_sobol.pdf
(+ PNG preview).  The digest for the paper writer, results/global_summary.md, is written by run_structural.py
(AN10, which runs after AN9 in run_all.py) from both JSON files.

Run:  cd /home/user/MCM && PYTHONPATH=src:analysis python3 analysis/run_sobol.py [--jobs 2] [--figure-only]
"""
from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import replace
from multiprocessing import get_context

import numpy as np

from phonebatt import RP63, aged, analytic, cell, run_batch
from phonebatt.params import K0, SOBOL_GROUPS, UQ_SPEC, map_unit, with_theta
from phonebatt.power import power_parts, steady_power
from phonebatt.usage import PRIMITIVES, Schedule, s9_schedule
from basics_common import RESULTS, meta, write_json

SEED = 2026
N_SUR = 1024  # surrogate scenarios
N_L2 = 256  # full-L2 scenarios (S7, S9, S10a)
N_G = 512  # interaction study S10a-G (k = 26: 27,648 runs)
N_CHECK = 200  # surrogate-vs-L2 check points per surrogate scenario
N_BOOT = 1000
H = 2.0  # rk4 step for every full-L2 run (§6.3 default)
CHECK_TOL_PCT = 1.0  # §10.3: surrogate error must be <= 1 %

GROUPS = list(SOBOL_GROUPS)  # 24 group labels, §10.3 order
GROUP_LABEL = {  # short labels for the figure / table
    "Q_n": "Q_n", "R0_ref": "R0,ref", "A_R": "A_R", "Ea_R0": "Ea,R0", "RC": "RC (R1,τ1,Ea)", "OCV": "OCV (δU,k_S)",
    "V_cut": "V_cut", "d_sys": "Δ_sys", "C_th": "C_th", "R_th": "R_th", "R_skin": "R_skin",
    "throttle": "throttle (T_lim,w,P_fl)", "k_env": "k_env", "floor": "floor (P_susp,P_awake,E_wake)",
    "display": "display (P_d0,p_c,γ_b)", "k_RR": "κ_RR", "SoC": "SoC (P_max,v0)", "wifi": "Wi-Fi radio",
    "lte": "LTE radio", "nr": "NR radio", "gamma_sig": "γ_sig", "peripherals": "peripherals", "m_b": "m_b",
    "gamma_R": "γ_R", "T_amb": "T_amb", "SOH": "SOH"}

# constant-load (surrogate) scenarios: id -> (primitive, label)
SUR_SCEN = {
    "S3": (PRIMITIVES["RD"], "reading"),
    "S4": (PRIMITIVES["SO"], "social, Wi-Fi"),
    "S5a": (PRIMITIVES["VW"], "video, Wi-Fi"),
    "S5b": (PRIMITIVES["VL"], "video, LTE"),
    "S5c": (PRIMITIVES["V5"], "video, 5G NSA"),
    "S6": (PRIMITIVES["NV"], "navigation"),
    "S8a": (PRIMITIVES["VC"], "video call, LTE"),
    "S8b": (PRIMITIVES["VC-W"], "video call, Wi-Fi"),
    "S11a": (replace(PRIMITIVES["SO-L"], RSRP=-110.0), "social, LTE, -110 dBm"),
    "S11c": (replace(PRIMITIVES["VC"], RSRP=-110.0), "video call, LTE, -110 dBm"),
}
# full-L2 scenarios
L2_SCEN = {
    "S7": dict(label="gaming, 25 °C", prim="GM", Tamb=25.0, soh=1.0, T0=None),
    "S9": dict(label="cold hike, −10 °C, SOH 0.8", prim=None, Tamb=-10.0, soh=0.8, T0=20.0),
    "S10a": dict(label="social, −10 °C, SOH 0.8", prim="SO", Tamb=-10.0, soh=0.8, T0=None),
    # S11c is also run in full L2: its surrogate fails the 1 % check at the thermal-shutdown boundary (P > 7 W)
    "S11c": dict(label="video call, LTE, −110 dBm", prim=None, u=SUR_SCEN["S11c"][0], Tamb=25.0, soh=1.0,
                 T0=None),
}
G_FACTORS = dict(T_amb=(-20.0, 10.0), SOH=(0.75, 1.0), m_b=(0.5, 2.0))  # S10a-G explicit factors (uniform)


# ------------------------------------------------------------------------------------------- design
def problem_groups(groups: list[str], extra: dict | None = None) -> tuple[dict, dict]:
    """SALib problem on the unit cube with SALib groups = our factor groups, and the ordered §4 spec."""
    names, glab, spec = [], [], {}
    for g in groups:
        for k in SOBOL_GROUPS[g]:
            names.append(k)
            glab.append(g)
            spec[k] = UQ_SPEC[k]
    for k, (lo, hi) in (extra or {}).items():
        names.append(k)
        glab.append(k)
        spec[k] = ("U", lo, 0.5 * (lo + hi), hi)
    prob = dict(num_vars=len(names), names=names, groups=glab, bounds=[[0.0, 1.0]] * len(names))
    return prob, spec


def design(groups, N, extra=None, seed=SEED):
    from SALib.sample import sobol as ssob
    prob, spec = problem_groups(groups, extra)
    U = ssob.sample(prob, N, calc_second_order=True, scramble=True, seed=seed)
    theta = map_unit(spec, U)
    return prob, theta, U


def factor_list(prob) -> list[str]:
    """Unique group labels in SALib order (= order of first appearance)."""
    seen = []
    for g in prob["groups"]:
        if g not in seen:
            seen.append(g)
    return seen


# ------------------------------------------------------------------------------------------- estimators
def split_saltelli(Y, D):
    """Rows of a second-order Saltelli design: A (N,), B (N,), AB (N, D), BA (N, D)."""
    step = 2 * D + 2
    N = Y.size // step
    Yr = Y.reshape(N, step)
    return Yr[:, 0], Yr[:, -1], Yr[:, 1:D + 1], Yr[:, D + 1:2 * D + 1]


def estimators(A, B, AB):
    """Sobol indices from a Saltelli design (vectorised over a leading bootstrap axis if present).

    S1: Saltelli et al. (2010) estimator  V_i ≈ E[f(B)·(f(A_B^i) − f(A))]  (SALib's first-order estimator);
    ST: Jansen (1999) estimator          V_Ti ≈ E[(f(A) − f(A_B^i))²]/2 (SALib's total-order estimator);
    S1_J: Jansen (1999) first order      V_i ≈ V − E[(f(B) − f(A_B^i))²]/2 (kept as a cross-check; its variance is
    far larger for small S_i, e.g. CI half-width 0.12 for every group of S10a-G at N = 256).
    V = Var(f(A) ∪ f(B)); outputs are centred on the mean of f(A) ∪ f(B) first.
    """
    AB_ = np.concatenate([A, B], axis=-1)
    V = np.var(AB_, axis=-1)[..., None]
    mu = np.mean(AB_, axis=-1)[..., None]  # centring: same expectation, far smaller variance (as SALib does)
    A, B, AB = A - mu, B - mu, AB - mu[..., None]
    S1 = np.mean(B[..., None] * (AB - A[..., None]), axis=-2) / V
    ST = 0.5 * np.mean((A[..., None] - AB) ** 2, axis=-2) / V
    S1J = 1.0 - 0.5 * np.mean((B[..., None] - AB) ** 2, axis=-2) / V
    return S1, ST, S1J


def sobol_indices(Y, prob, n_boot=N_BOOT, seed=SEED, second_order=False) -> dict:
    """S1 (Saltelli 2010) and ST (Jansen) with 95 % percentile-bootstrap CIs; optional SALib S2 (Saltelli 2002)."""
    labels = factor_list(prob)
    D = len(labels)
    Y = np.asarray(Y, dtype=float)
    A, B, AB, _ = split_saltelli(Y, D)
    S1, ST, S1J = estimators(A, B, AB)
    rng = np.random.default_rng(seed)
    N = A.size
    b1, bT = [], []
    for _ in range(0, n_boot, 100):
        idx = rng.integers(0, N, size=(100, N))
        s1, st, _ = estimators(A[idx], B[idx], AB[idx])
        b1.append(s1)
        bT.append(st)
    b1, bT = np.concatenate(b1), np.concatenate(bT)
    lo1, hi1 = np.percentile(b1, [2.5, 97.5], axis=0)
    loT, hiT = np.percentile(bT, [2.5, 97.5], axis=0)
    out = dict(
        S1={g: float(v) for g, v in zip(labels, S1)}, ST={g: float(v) for g, v in zip(labels, ST)},
        S1_conf={g: float(0.5 * (h - l)) for g, l, h in zip(labels, lo1, hi1)},
        ST_conf={g: float(0.5 * (h - l)) for g, l, h in zip(labels, loT, hiT)},
        S1_ci={g: [float(l), float(h)] for g, l, h in zip(labels, lo1, hi1)},
        ST_ci={g: [float(l), float(h)] for g, l, h in zip(labels, loT, hiT)},
        S1_jansen={g: float(v) for g, v in zip(labels, S1J)},
        sum_S1=float(S1.sum()), sum_ST=float(ST.sum()), sum_S1_jansen=float(S1J.sum()),
        interaction_share=float(1.0 - S1.sum()),
        interaction_share_ci=[float(v) for v in np.percentile(1.0 - b1.sum(axis=1), [2.5, 97.5])],
        mean=float(np.mean(np.r_[A, B])), sd=float(np.std(np.r_[A, B])), N=int(N), D=int(D), n_eval=int(Y.size),
        q05_q50_q95=[float(v) for v in np.percentile(Y, [5, 50, 95])])
    out["cv"] = out["sd"] / out["mean"] if out["mean"] else math.nan
    if second_order:
        from SALib.analyze import sobol as asob
        r = asob.analyze(prob, Y, calc_second_order=True, num_resamples=n_boot, conf_level=0.95, seed=seed)
        out["S1_salib"] = {g: float(v) for g, v in zip(labels, r["S1"])}
        out["ST_salib"] = {g: float(v) for g, v in zip(labels, r["ST"])}
        out["_S2_matrix"] = (np.asarray(r["S2"]), np.asarray(r["S2_conf"]), labels)
    return out


# ------------------------------------------------------------------------------------------- surrogate (F2+F4+F5)
def _gl(n=200):
    return np.polynomial.legendre.leggauss(n)


def _node_steady(z, Pns, dem, Tamb_K, p, n_iter=48):
    """Quasi-steady thermal state at SOC z (F5 including the battery heat q_batt, solved by bisection on T_s).

    T_s = T_amb + k_env R_skin [P(T_s) + q_batt],  T = T_amb + (R_th/R_skin)(T_s − T_amb)  (E8),
    q_batt = I²(R0 + R1) − I T k_S ∂U/∂T(z) with the quasi-static current I = P/V_qs (RC equilibrated).
    Arrays broadcast as (m, n) (m SOC nodes per sample) or (n,).  Returns (P, T_K).
    """
    from phonebatt.power import throttle

    def state(Ts):
        P = Pns + throttle(dem, Ts, p)
        T = Tamb_K + p.R_th / p.R_skin * (Ts - Tamb_K)
        I = P / analytic.v_qs(z, P, T, p)
        qb = I * I * (cell.R0(z, T, p) + cell.R1(T, p)) - I * T * p.k_S * cell.dUdT(z, p)
        return P, T, qb

    shape = np.broadcast(np.asarray(z), np.asarray(Pns), np.asarray(dem), np.asarray(Tamb_K)).shape
    lo = np.full(shape, 0.0) + (Tamb_K - 1.0)
    hi = np.full(shape, 0.0) + (Tamb_K + 200.0)
    for _ in range(n_iter):
        mid = 0.5 * (lo + hi)
        P, T, qb = state(mid)
        pos = mid - Tamb_K - p.k_env * p.R_skin * (P + qb) > 0
        hi = np.where(pos, mid, hi)
        lo = np.where(pos, lo, mid)
    P, T, _ = state(0.5 * (lo + hi))
    return P, T


def surrogate(u, theta: dict, Tamb_C=25.0, soh=1.0, chunk=4096, n_gl=200, n_prof=24, with_qbatt=True) -> dict:
    """Vectorised constant-load surrogate (§10.3): F2 at the thermal steady state, F5 power, F4 end point.

    The battery temperature and the (throttled) power are the quasi-steady fixed point of F5 including the battery
    heat q_batt; because q_batt grows as z falls, they are evaluated as profiles T(z), P(z) on ``n_prof`` SOC nodes
    (thermal time constant 7.5 min ≪ TTE) and interpolated onto the Gauss–Legendre nodes:
    TTE = Q_max ∫_{z_end}^{z0} V_qs(z)/P(z) dz,  z_end = max(z_c, z*) with z_c the cut-off root (F2) and z* the F4
    brown-out root (P + ΔP_b, v1 = R1 P/U), both at the end-point temperature.  ``with_qbatt=False`` gives the plain
    spec form (q_batt neglected, T = T_amb + k_env R_th P*, constant P*).
    """
    n = next(iter(theta.values())).shape[-1]
    tte = np.empty(n)
    soff = np.empty(n)
    zend_is_bo = np.empty(n, bool)
    therm = np.empty(n, bool)
    Pst = np.empty(n)
    xg, wg = _gl(n_gl)
    xi_g = 0.5 * (xg + 1.0)  # GL nodes on [0, 1]
    xi_k = np.linspace(0.0, 1.0, n_prof) ** 2  # profile nodes, dense near the end point
    Wint = np.zeros((n_gl, n_prof))  # linear interpolation weights xi_k -> xi_g
    j = np.clip(np.searchsorted(xi_k, xi_g) - 1, 0, n_prof - 2)
    f = (xi_g - xi_k[j]) / (xi_k[j + 1] - xi_k[j])
    Wint[np.arange(n_gl), j] = 1.0 - f
    Wint[np.arange(n_gl), j + 1] = f
    Ta = Tamb_C + K0
    for i0 in range(0, n, chunk):
        sl = slice(i0, min(n, i0 + chunk))
        th = {k: v[..., sl] for k, v in theta.items()}
        p = with_theta(aged(RP63, 1.0 - soh) if soh < 1 else RP63, th)
        Pns, dem, dPb, _ = power_parts(u, p)
        P = np.asarray(steady_power(u, p, Tamb_C), float)  # F5 fixed point with q_batt = 0
        m = P.shape
        Pns, dem = np.broadcast_to(np.asarray(Pns, float), m), np.broadcast_to(np.asarray(dem, float), m)
        dPb = np.broadcast_to(np.asarray(dPb, float), m)
        T_K = Ta + p.k_env * p.R_th * P
        z0 = np.broadcast_to(cell.z_from_s0(1.0, p), m) + 0.0

        def end_point(Pe, Te):
            zc = analytic.cutoff_soc(Pe, Te, p, z0)
            if np.any(dPb > 0):
                zb = np.where(dPb > 0, analytic.brownout_root(Pe, dPb, Te, p)["z_star"], 0.0)
            else:
                zb = np.zeros(m)
            return np.minimum(np.maximum(zc, zb), z0), zb > zc

        zend, isbo = end_point(P, T_K)
        Tsh = np.broadcast_to(np.asarray(p.T_shut, float), m)
        thermal = np.zeros(m, bool)
        hot = np.zeros(m, bool)
        if with_qbatt:
            for _ in range(2):  # end-point temperature with q_batt (2 passes converge to < 1e-4 in z)
                Pe, Te = _node_steady(zend, Pns, dem, Ta, p)
                zend, isbo = end_point(Pe, Te)
            zk = zend[None, :] + (z0 - zend)[None, :] * xi_k[:, None]
            Pk, Tk = _node_steady(zk, Pns[None, :], dem[None, :], Ta, p)
            # thermal shutdown (event g5): the quasi-steady temperature reaches T_shut somewhere on [z_end, z0]
            over = Tk >= Tsh[None, :]
            hot = over[-1]  # already above T_shut at full charge: shutdown during the warm-up transient (F5)
            cross = over.any(axis=0) & ~hot
            if cross.any():  # T(z) rises as z falls: end at the highest node still below/above the limit
                kk = n_prof - 1 - np.argmax(over[::-1], axis=0)  # highest node index with T >= T_shut
                kk = np.minimum(kk, n_prof - 2)
                c_ = np.arange(m[0])
                T1, T2 = Tk[kk, c_], Tk[kk + 1, c_]
                fr = np.clip((T1 - Tsh) / np.where(T1 > T2, T1 - T2, 1.0), 0.0, 1.0)
                xi_th = xi_k[kk] + fr * (xi_k[kk + 1] - xi_k[kk])
                z_th = zend + (z0 - zend) * xi_th
                zend = np.where(cross, z_th, zend)
                thermal = cross
                isbo = isbo & ~cross
                zk = zend[None, :] + (z0 - zend)[None, :] * xi_k[:, None]
                Pk, Tk = _node_steady(zk, Pns[None, :], dem[None, :], Ta, p)
            Pq, Tq = Wint @ Pk, Wint @ Tk
        else:
            Pq, Tq = P[None, :], T_K[None, :]
        half = 0.5 * (z0 - zend)
        zq = zend[None, :] + 2.0 * half[None, :] * xi_g[:, None]
        Vq = analytic.v_qs(zq, Pq, Tq, p)
        tte_c = p.Q_max * half * np.sum(wg[:, None] * Vq / Pq, axis=0)
        ze = cell.z_empty(p)
        soff_c = (zend - ze) / (1.0 - ze)
        if hot.any():  # F5 step response from T0 = T_amb towards T_inf = T(z0) > T_shut
            Tinf = Tk[-1]
            tau = p.k_env * p.R_th * p.C_th
            t_th = tau * np.log(np.maximum(Tinf - Ta, 1e-9) / np.maximum(Tinf - Tsh, 1e-9))
            I0 = Pk[-1] / analytic.v_qs(z0, Pk[-1], Tinf, p)
            z_th = z0 - I0 * t_th / (3600.0 * p.Q_max)
            tte_c = np.where(hot, t_th / 3600.0, tte_c)
            soff_c = np.where(hot, (z_th - ze) / (1.0 - ze), soff_c)
            isbo = isbo & ~hot
        tte[sl] = tte_c
        soff[sl] = soff_c
        zend_is_bo[sl] = isbo
        therm[sl] = thermal | hot
        Pst[sl] = (Wint @ Pk).mean(axis=0) if with_qbatt else P
    return dict(tte_h=tte, s_off=soff, brownout=zend_is_bo, thermal=therm, P_W=Pst)


# ------------------------------------------------------------------------------------------- full L2 jobs
def _l2_job(args):
    """Worker: one full-L2 batch.  args = (tag, kind, theta, extra) -> dict of arrays."""
    tag, kind, theta, extra = args
    t0 = time.time()
    if kind == "const":
        u = extra["u"]
        b = run_batch(Schedule.constant(u), extra.get("Tamb", 25.0), theta, h=H, p=RP63)
    elif kind == "l2":
        sc = L2_SCEN[extra["scen"]]
        p = aged(RP63, 1.0 - sc["soh"]) if sc["soh"] < 1 else RP63
        if extra["scen"] == "S9":
            sch = s9_schedule(sc["Tamb"])
        else:
            sch = Schedule.constant(sc["u"] if "u" in sc else PRIMITIVES[sc["prim"]])
        b = run_batch(sch, sc["Tamb"], theta, h=H, T0_C=sc["T0"], p=p)
    elif kind == "G":
        th = dict(theta)
        Ta = th.pop("T_amb")
        soh = th.pop("SOH")
        th["q_loss"] = 1.0 - soh
        b = run_batch(Schedule.constant(PRIMITIVES["SO"]), Ta, th, h=H, p=RP63)
    else:
        raise ValueError(kind)
    return tag, dict(tte_h=np.asarray(b["tte_h"]), s_off=np.asarray(b["s_off"]), event=np.asarray(b["event"]),
                     T_end_C=np.asarray(b["T_end_C"])), time.time() - t0


def _split_theta(theta, nparts):
    n = next(iter(theta.values())).shape[-1]
    edges = np.linspace(0, n, nparts + 1).astype(int)
    return [{k: v[..., a:b] for k, v in theta.items()} for a, b in zip(edges[:-1], edges[1:])]


# ------------------------------------------------------------------------------------------- main
def run_models(jobs: int = 2) -> dict:
    """Evaluate every design point (surrogate, surrogate check, full L2, interaction study); returns raw outputs."""
    t_all = time.time()
    timing: dict = {}

    # ---- designs
    prob24, th_sur, _ = design(GROUPS, N_SUR)
    prob24_l2, th_l2, _ = design(GROUPS, N_L2)
    groups_G = [g for g in GROUPS if g != "m_b"]
    probG, th_G, _ = design(groups_G, N_G, extra=G_FACTORS)
    # m_b is an explicit factor in S10a-G: rename the uniform factor into the parameter field
    th_G = dict(th_G)

    # ---- surrogate Sobol
    sur_out = {}
    for s, (u, lab) in SUR_SCEN.items():
        t0 = time.time()
        r = surrogate(u, th_sur)
        sur_out[s] = r
        timing[f"surrogate_{s}"] = round(time.time() - t0, 1)
        print(f"surrogate {s}: {time.time() - t0:.1f} s, mean TTE {r['tte_h'].mean():.2f} h", flush=True)

    # ---- L2 jobs (surrogate check + full-L2 Sobol + interaction study), process pool
    rng = np.random.default_rng(SEED + 1)
    n_sur = next(iter(th_sur.values())).shape[-1]
    check_idx = {s: np.sort(rng.choice(n_sur, N_CHECK, replace=False)) for s in SUR_SCEN}
    tasks = []
    for s, (u, lab) in SUR_SCEN.items():
        th = {k: v[..., check_idx[s]] for k, v in th_sur.items()}
        tasks.append((f"check|{s}", "const", th, dict(u=u)))
    nsplit = 2  # split the big batches so the pool balances
    for s in L2_SCEN:
        for i, th in enumerate(_split_theta(th_l2, nsplit)):
            tasks.append((f"l2|{s}|{i}", "l2", th, dict(scen=s)))
    for i, th in enumerate(_split_theta(th_G, nsplit)):
        tasks.append((f"G|{i}", "G", th, {}))
    # largest first
    order = sorted(range(len(tasks)), key=lambda i: (0 if tasks[i][1] in ("l2", "G") else 1))
    tasks = [tasks[i] for i in order]
    outs = {}
    if jobs > 1:
        with get_context("fork").Pool(jobs) as pool:
            for tag, o, dt in pool.imap_unordered(_l2_job, tasks):
                outs[tag] = o
                timing[tag] = round(dt, 1)
                print(f"done {tag}: {dt:.1f} s", flush=True)
    else:
        for t in tasks:
            tag, o, dt = _l2_job(t)
            outs[tag] = o
            timing[tag] = round(dt, 1)
            print(f"done {tag}: {dt:.1f} s", flush=True)

    def merged(prefix):
        parts = [outs[f"{prefix}|{i}"] for i in range(nsplit)]
        return {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}

    l2 = {s: merged(f"l2|{s}") for s in L2_SCEN}
    l2["G"] = merged("G")
    checks_raw = {s: outs[f"check|{s}"] for s in SUR_SCEN}
    return dict(sur_out=sur_out, check_idx=check_idx, checks_raw=checks_raw, l2=l2, timing=timing,
                runtime_s=round(time.time() - t_all, 1))


def analyse(raw: dict) -> dict:
    """Sobol indices, CIs and checks from the raw outputs of :func:`run_models`."""
    res: dict = {}
    prob24, _, _ = design(GROUPS, N_SUR)
    prob24_l2 = prob24  # same factor structure (groups), only N differs
    probG, _, _ = design([g for g in GROUPS if g != "m_b"], N_G, extra=G_FACTORS)
    sur_out, check_idx, l2 = raw["sur_out"], raw["check_idx"], raw["l2"]
    outs = {f"check|{s}": o for s, o in raw["checks_raw"].items()}

    # ---- surrogate scenarios: indices + check
    checks = {}
    for s, (u, lab) in SUR_SCEN.items():
        r = sur_out[s]
        d = sobol_indices(r["tte_h"], prob24)
        d.update(label=lab, method="surrogate F2 at thermal steady state + F4 brown-out end point + F5 power",
                 Tamb_C=25.0, soh=1.0, s0=1.0, brownout_end_share=float(np.mean(r["brownout"])),
                 thermal_end_share=float(np.mean(r["thermal"])),
                 P_mean_W=float(np.mean(r["P_W"])))
        o = outs[f"check|{s}"]
        tl2 = o["tte_h"]
        ts = r["tte_h"][check_idx[s]]
        err = 100.0 * (ts / tl2 - 1.0)
        so_err = 100.0 * (r["s_off"][check_idx[s]] - o["s_off"])
        checks[s] = dict(n=int(N_CHECK), maxabs_err_pct=float(np.max(np.abs(err))), mean_err_pct=float(np.mean(err)),
                         n_over_1pct=int(np.sum(np.abs(err) > CHECK_TOL_PCT)),
                         over_1pct_events_L2=[str(e) for e in o["event"][np.abs(err) > CHECK_TOL_PCT]],
                         over_1pct_P_W=[float(v) for v in r["P_W"][check_idx[s]][np.abs(err) > CHECK_TOL_PCT]],
                         rms_err_pct=float(np.sqrt(np.mean(err ** 2))), q05_q95_err_pct=[float(v) for v in
                                                                                       np.percentile(err, [5, 95])],
                         s_off_maxabs_err_pp=float(np.max(np.abs(so_err))),
                         event_share_L2={e: float(np.mean(o["event"] == e)) for e in np.unique(o["event"])},
                         pass_le_1pct=bool(np.max(np.abs(err)) <= CHECK_TOL_PCT))
        d["surrogate_check"] = checks[s]
        res[s if s not in L2_SCEN else f"{s}_surrogate"] = d

    # ---- full-L2 scenarios
    for s, sc in L2_SCEN.items():
        o = l2[s]
        d = sobol_indices(o["tte_h"], prob24_l2)
        if s in SUR_SCEN:
            c = checks[s]
            d["note"] = (f"full L2 replaces the surrogate for this scenario: surrogate check max |error| "
                         f"{c['maxabs_err_pct']:.2f} % ({c['n_over_1pct']} of {c['n']} points above "
                         f"{CHECK_TOL_PCT:g} %, L2 events {sorted(set(c['over_1pct_events_L2']))}; the TTE jumps at "
                         f"the thermal-shutdown boundary); surrogate indices kept under '{s}_surrogate'")
        d.update(label=sc["label"], method="full L2 (run_batch, rk4_batch, h = 2 s)", Tamb_C=sc["Tamb"],
                 soh=sc["soh"], s0=1.0, T0_C=sc["T0"] if sc["T0"] is not None else sc["Tamb"],
                 event_share={e: float(np.mean(o["event"] == e)) for e in np.unique(o["event"])},
                 n_nonfinite=int(np.sum(~np.isfinite(o["tte_h"]))))
        so = o["s_off"]
        d["s_off"] = sobol_indices(100.0 * so, prob24_l2)
        d["s_off"]["unit"] = "percentage points (D1)"
        res[s] = d

    # ---- interaction study S10a-G
    o = l2["G"]
    dG = sobol_indices(o["tte_h"], probG, second_order=True)
    dGs = sobol_indices(100.0 * o["s_off"], probG, second_order=True)
    trio = list(G_FACTORS)
    for d in (dG, dGs):
        S2m, S2c, labs = d.pop("_S2_matrix")
        idx = {g: i for i, g in enumerate(labs)}
        d["S2"] = {}
        d["S2_conf"] = {}
        for i, a in enumerate(trio):
            for b in trio[i + 1:]:
                ia, ib = sorted((idx[a], idx[b]))
                d["S2"][f"{a}|{b}"] = float(S2m[ia, ib])
                d["S2_conf"][f"{a}|{b}"] = float(S2c[ia, ib])
        d["trio_sum_S1"] = float(sum(d["S1"][g] for g in trio))
        d["trio_sum_S2"] = float(sum(d["S2"].values()))
    dG.update(label="social S4: T_amb ~ U(−20, 10) °C, SOH ~ U(0.75, 1.0), m_b ~ U(0.5, 2.0) explicit",
              method="full L2 (run_batch, rk4_batch, h = 2 s)", explicit_factors={k: list(v) for k, v in
                                                                                   G_FACTORS.items()},
              event_share={e: float(np.mean(o["event"] == e)) for e in np.unique(o["event"])})
    dGs["unit"] = "percentage points (D1)"
    dG["s_off"] = dGs
    res["S10a-G"] = dG

    res["surrogate_check_maxerr_pct"] = float(max(c["maxabs_err_pct"] for c in checks.values()))
    used = {s: c for s, c in checks.items() if s not in L2_SCEN}
    res["surrogate_check_maxerr_pct_used"] = float(max(c["maxabs_err_pct"] for c in used.values()))
    res["surrogate_check_note"] = ("max |error| over all 10 surrogate scenarios (200 points each); '_used' excludes "
                                   "the scenarios whose reported indices come from full L2 instead: "
                                   + ", ".join(s for s in checks if s in L2_SCEN))
    res["surrogate_check"] = checks
    res["_timing_s"] = raw["timing"]
    res["_runtime_s"] = raw["runtime_s"]
    res["_groups"] = {g: list(SOBOL_GROUPS[g]) for g in GROUPS}
    return res


# ------------------------------------------------------------------------------------------- figure
FIG_COLS = ["S3", "S4", "S5a", "S5b", "S5c", "S6", "S8a", "S8b", "S11a", "S11c", "S7", "S9", "S10a", "S10a-G"]
FIG_FIRST_L2 = "S11c"
ROW_MIN = 0.05  # factor rows shown in fig_sobol: max S_Ti over the columns >= ROW_MIN


def make_figure(res: dict):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from style import INK, INK2, MUTED, SEQ_BLUE, W1, apply, panel_label, save

    apply()
    cols = [c for c in FIG_COLS if c in res]
    rows_all = GROUPS + ["T_amb", "SOH"]
    # keep the factor rows that matter somewhere (S_Ti >= 0.02); list the rest in the note
    def st(c, g):
        return res[c]["ST"].get(g, np.nan)

    keep = [g for g in rows_all if np.nanmax([st(c, g) for c in cols]) >= ROW_MIN]
    dropped = [g for g in rows_all if g not in keep]
    M1 = np.array([[res[c]["S1"].get(g, np.nan) for c in cols] for g in keep])
    MT = np.array([[res[c]["ST"].get(g, np.nan) for c in cols] for g in keep])
    nr, nc = M1.shape
    fig, axes = plt.subplots(1, 2, figsize=(W1, 0.15 * nr + 0.95), sharey=True,
                             gridspec_kw=dict(wspace=0.04))
    norm = Normalize(0.0, 1.0)
    for ax, M, title, tag in ((axes[0], M1, "First-order $S_i$", "(a)"), (axes[1], MT, "Total $S_{Ti}$", "(b)")):
        Mp = np.ma.masked_invalid(np.clip(M, 0.0, 1.0))
        cm = LinearSegmentedColormap.from_list("seq_blue_w", ["#f7fafe"] + SEQ_BLUE)  # single hue, 0 = near white
        cm.set_bad("#ecebe7")
        im = ax.imshow(Mp, cmap=cm, norm=norm, aspect="auto", interpolation="nearest")
        for i in range(nr):
            for j in range(nc):
                v = M[i, j]
                if np.isfinite(v) and v >= 0.05:
                    ax.text(j, i, f"{v:.2f}".lstrip("0") if v < 1 else "1.0", ha="center", va="center",
                            fontsize=5.6, color="white" if v > 0.55 else INK)
        ax.set_xticks(range(nc))
        ax.set_xticklabels(cols, rotation=60, ha="right", fontsize=6.6)
        ax.set_yticks(range(nr))
        ax.set_yticklabels([GROUP_LABEL.get(g, g) for g in keep], fontsize=6.6)
        ax.tick_params(length=0)
        ax.grid(False)
        for sp in ax.spines.values():
            sp.set_visible(False)
        # separators: surrogate | full L2 | interaction study
        for x in (cols.index(FIG_FIRST_L2) - 0.5 if FIG_FIRST_L2 in cols else None,
                  cols.index("S10a-G") - 0.5 if "S10a-G" in cols else None):
            if x is not None:
                ax.axvline(x, color=MUTED, lw=0.8)
        ax.set_title(title, fontsize=8.5, loc="left", color=INK)
        panel_label(ax, tag)
    cb = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.01)
    cb.ax.tick_params(labelsize=6.5, length=2)
    cb.outline.set_visible(False)
    cb.set_label("Sobol index (TTE)", fontsize=7, color=INK2)
    note = ("S3–S11a: surrogate (N = 1024); S11c, S7, S9, S10a: full L2 (N = 256); S10a-G: T_amb, SOH explicit. "
            f"Rows with S_Ti < {ROW_MIN:g} in every column omitted: "
            f"{', '.join(GROUP_LABEL.get(g, g) for g in dropped)}. Gray: factor not in that design.")
    fig.text(0.0, -0.02, note, fontsize=5.8, color=MUTED, ha="left", va="top", wrap=True)
    save(fig, "fig_sobol")
    return keep, dropped


# ------------------------------------------------------------------------------------------- entry
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=2, help="worker processes for the L2 batches (default 2)")
    ap.add_argument("--figure-only", action="store_true", help="re-draw fig_sobol from results/sobol.json")
    ap.add_argument("--save-raw", default="", help="pickle the raw model outputs to this path (development)")
    ap.add_argument("--from-raw", default="", help="re-analyse a pickle written by --save-raw (no model runs)")
    a = ap.parse_args(argv)
    if a.figure_only:
        res = json.loads((RESULTS / "sobol.json").read_text())
        make_figure(res)
        return 0
    import pickle
    if a.from_raw:
        with open(a.from_raw, "rb") as fh:
            raw = pickle.load(fh)
    else:
        raw = run_models(jobs=a.jobs)
        if a.save_raw:
            with open(a.save_raw, "wb") as fh:
                pickle.dump(raw, fh)
    res = analyse(raw)
    keep, dropped = make_figure(res)
    out = dict(meta=meta("run_sobol.py", seed=SEED, solver="surrogate F2+F4+F5 (vectorised); full L2: run_batch "
                                                            "(rk4_batch) h = 2 s", h_s=H,
                         N_surrogate=N_SUR, N_L2=N_L2, N_interaction=N_G, n_check=N_CHECK, n_boot=N_BOOT,
                         estimator="S1: Saltelli et al. (2010); ST: Jansen (1999) (SALib's pairing) on SALib's "
                                   "second-order Saltelli design; S1_jansen: Jansen first order (cross-check); "
                                   "S2: SALib (Saltelli 2002); *_conf = half-width of the 95 % percentile-bootstrap "
                                   "CI (1000 resamples of the N base rows)",
                         design="SALib.sample.sobol.sample(scramble=True, seed=2026) on the unit cube, mapped by "
                                "phonebatt.params.map_unit through the FINAL_SPEC §4 distributions",
                         runtime_s=res.pop("_runtime_s"), timing_s=res.pop("_timing_s"),
                         fig_rows_omitted=dropped, output_unit="TTE in h; s_off blocks in percentage points"),
               groups=res.pop("_groups"))
    out.update(res)
    write_json("sobol.json", out)
    print(f"wrote results/sobol.json; surrogate max error {out['surrogate_check_maxerr_pct']:.3f} %")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
