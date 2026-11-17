"""AN11: usage-fluctuation sensitivity (FINAL_SPEC §10.5): how much does *how* the phone is used matter, beyond
the mean power?

Writes results/usage_fluct.json with the §14 keys
* ``rate_scaling``  — §10.5.1: all sojourns × κ, π fixed; F7 predicts SD ∝ √κ (taken from results/ctmc.json when it
  was produced with the same configuration, else recomputed with ``run_ctmc.rate_scaling``).
* ``burstiness``    — §10.5.2: constant 1.2 W vs a 50/50 alternation 0.4 W / 2.0 W, periods 10 s … 1 h, at 25 °C
  fresh and −10 °C SOH 0.8, without and with a 3 W burst envelope; L2 (LSODA) vs L0 and the L1 quasi-static
  limits (period ≪ τ1: R0 sees the alternation, R1 the mean; period ≫ τ1: R0 + R1 see it) and F10 (Jensen).
* ``order_effect``  — §10.5.3: S12 forward vs reversed, at 25 °C fresh and −10 °C SOH 0.8, L2 and L0; plus the
  start-time (phase) scan of both orders.
* ``s0_curvature``  — §10.5.4: TTE(s0) of S4 (social, ΔP_b = 3 W) at 25/0/−10/−20 °C fresh and −10/−20 °C SOH 0.8;
  linear fit on s0 ≥ 0.2, effective empty point vs F4/F2.
* ``heavy_tail``    — §10.5.5 / F9 / X15: lognormal sessions (CV 0.5 … 3): F9 vs L0 Monte Carlo (full chain and the
  exact hub check, also at a 50× energy budget), and L2 Monte Carlo at CV 1.5 coupled to exponential paths.
* ``jitter``        — §10.5.5: per-visit ±30 % jitter of b, APL, λ, d_dl: extended F6/F7 vs paired L0 and L2 MC.
* ``usage_params``  — §10.5.6: λ_on, p (Dirichlet 50·p), d_j, connectivity means with the §4.6 ranges through F6–F8:
  distribution of the Dynkin mean and CLT SD, law of total variance (path vs usage-model parameters), OAT
  tornado, first-order indices, F8-predicted vs exact OAT.

Figure: paper/figures/fig_usage_fluct.pdf.
Run:  cd /home/user/MCM && PYTHONPATH=src:analysis python3 analysis/run_usage_fluct.py      (about 6 min)
"""
from __future__ import annotations

import json
import math
import time

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from basics_common import RESULTS, meta, write_json
from phonebatt import RP63, aged, analytic, cell, run, run_batch
from phonebatt.cell import R0, R1, ocv
from phonebatt.params import K0, USAGE_UQ_SPEC, sample
from phonebatt.power import demand_power, steady_power
from phonebatt.usage import ACTS, CTMC, PRIMITIVES, PROFILES, Schedule, UsageState, s12_items
from stochastic_common import (EU_WH, SESSIONS, ExpandedCTMC, UsageSim, f9_full_sigma2, hub_inputs, jitter_bank,
                               jitter_expanded_inputs, jitter_power_bank, l0_from_paths, state_powers, summ)
from style import (BLUE, CAT, GRID, INK, INK2, MUTED, ORANGE, SEQ_BLUE, W1, apply, panel_label, save)

PROFS = ("typical", "heavy")
CFG = dict(
    burst=dict(P_const=1.2, P_lo=0.4, P_hi=2.0, periods_s=[10, 30, 60, 120, 300, 600, 1800, 3600],
               conditions=[[25.0, 1.0], [-10.0, 0.8]], envelope_W=[0.0, 3.0]),
    order=dict(conditions=[[25.0, 1.0], [-10.0, 0.8]], phase_shifts_h=[0, 4, 8, 12, 16, 20]),
    s0=dict(scenario="S4", grid=[round(0.05 * i, 2) for i in range(1, 21)],
            conditions=[[25.0, 1.0], [0.0, 1.0], [-10.0, 1.0], [-20.0, 1.0], [-10.0, 0.8], [-20.0, 0.8]]),
    heavy_tail=dict(laws=[["exp", 1.0], ["lognormal", 0.5], ["lognormal", 1.0], ["lognormal", 1.5],
                          ["lognormal", 2.0], ["lognormal", 3.0]], n_L0=20000, seed_L0=505,
                    hub_cvs=[1.0, 1.5, 3.0], hub_energy_mult=[1, 50], n_hub=4000, seed_hub=11,
                    n_L2=400, seed_L2=606, cv_L2=1.5, h_s=10.0, t_max_h=200.0),
    jitter=dict(K=256, seed_bank=4242, K_ref=16384, lo=0.7, hi=1.3, n_L0=20000, seed_L0=707, n_L2=400,
                seed_L2=808, h_s=10.0, t_max_h=200.0),
    usage_params=dict(n=4096, seed=2026, oat_beta_q=[0.05, 0.95], n_bins=32),
    Tamb_C=25.0,
)


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


def const_load(P: float, dPb: float = 0.0, name: str = "") -> UsageState:
    """A pure constant-power input: suspend floor + explicit extra load (screen off, no radio, no wake-ups)."""
    return UsageState(name or f"C{P:g}", P_extra=P - float(RP63.P_susp), dPb=dPb)


def soh_params(soh: float):
    return RP63 if soh >= 1.0 else aged(RP63, 1.0 - soh)


def result_dict(r) -> dict:
    return dict(tte_h=r.tte_h, event=r.event, s_off_pct=100 * r.s_off, P_mean_W=r.P_mean_W,
                E_del_Wh=r.P_mean_W * r.tte_h, T_end_C=r.T_end_C)


# ------------------------------------------------------------------------------------------- 1 rate scaling
def rate_scaling_block() -> dict:
    import run_ctmc
    want = {k: run_ctmc.CFG[k] for k in ("kappas", "n_L0_kappa", "seed_L0_kappa", "n_L2_kappa", "seed_L2_kappa",
                                          "cv_LN", "n_LN", "seed_LN", "h_s")}
    src = RESULTS / "ctmc.json"
    rs, source = None, "computed by run_usage_fluct.py (run_ctmc.rate_scaling)"
    if src.exists():
        d = json.load(open(src))
        if d.get("rate_scaling", {}).get("config") == want:
            rs, source = d["rate_scaling"], "results/ctmc.json (AN5, identical configuration)"
    if rs is None:
        rs = {p: run_ctmc.rate_scaling(p) for p in PROFS}
        rs["kappa"] = run_ctmc.CFG["kappas"]
    out = dict(source=source, kappa=rs["kappa"], config=want)
    for prof in PROFS:
        rows = rs[prof]["rows"]
        out[prof] = dict(
            mu_dynkin_h=[r["mu_dynkin_h"] for r in rows], sd_clt_h=[r["sd_clt_h"] for r in rows],
            mc_L0_mean_h=[r["mc_L0"]["mean_h"] for r in rows], mc_L0_sd_h=[r["mc_L0"]["sd_h"] for r in rows],
            mc_L2_mean_h=[r["mc_L2"]["mean_h"] for r in rows], mc_L2_sd_h=[r["mc_L2"]["sd_h"] for r in rows],
            sd_mc_L0_over_clt=[r["sd_mc_L0_over_clt"] for r in rows],
            sd_mc_L2_over_clt=[r["sd_mc_L2_over_clt"] for r in rows],
            lognormal_sd_f9_h=[r["lognormal"]["sd_f9_approx_h"] for r in rows],
            lognormal_mc_L0_sd_h=[r["lognormal"]["mc_L0"]["sd_h"] for r in rows],
            loglog_slope=rs[prof]["loglog_slope"])
    out["proto"] = dict(kappa=[0.25, 1.0, 4.0], mu_dynkin_h=[40.63, 41.25, 43.71], sd_clt_h=[3.61, 7.21, 14.43])
    return out


# ------------------------------------------------------------------------------------------- 2 burstiness
def l1_alternation(p, Tamb_C: float, P_lo: float, P_hi: float, dPb: float, limit: str, s0: float = 1.0,
                   n: int = 6001) -> dict:
    """L1 quasi-static TTE of a 50/50 alternation in the two period limits (exact quadrature, F2 generalised).

    'long'  (τ1 ≪ period ≪ τ_th): each phase sees R_Σ = R0 + R1; cut-off when V(P_hi) = V_cut.
    'long_thermal' (period ≫ τ_th = 7.5 min as well): as 'long', but each phase at its own steady temperature
            T_amb + k_env R_th P_phase (the battery is warmer exactly when the load is high).
    'short' (τ_burst ≪ period ≪ τ1): v1 = R1·Ī (mean current), R0 sees each phase.
    Battery temperature (except 'long_thermal'): the steady self-heated T = T_amb + k_env R_th P̄ (as the §10.3
    surrogate).
    The end point is the first violation of the cut-off or (if dPb > 0) the brown-out (F4) in the high phase.
    TTE = Q_max ∫ dz / Ī [h].
    """
    Pbar = 0.5 * (P_lo + P_hi)
    T = Tamb_C + K0 + float(p.k_env * p.R_th) * Pbar
    z0 = float(cell.z_from_s0(s0, p))
    z = np.linspace(1e-4, z0, n)
    U = ocv(z, T, p)
    r0 = R0(z, T, p)
    r1 = float(R1(T, p))

    def cur(Ue, R, P):
        D = Ue * Ue - 4.0 * R * P
        return np.where(D > 0, 2.0 * P / (Ue + np.sqrt(np.maximum(D, 0.0))), np.nan), D

    if limit == "long":
        I_lo, _ = cur(U, r0 + r1, P_lo)
        I_hi, D = cur(U, r0 + r1, P_hi)
        Ibar = 0.5 * (I_lo + I_hi)
        Ue_hi = U - r1 * I_hi
    elif limit == "long_thermal":
        T_lo = Tamb_C + K0 + float(p.k_env * p.R_th) * P_lo
        T_hi = Tamb_C + K0 + float(p.k_env * p.R_th) * P_hi
        I_lo, _ = cur(ocv(z, T_lo, p), R0(z, T_lo, p) + float(R1(T_lo, p)), P_lo)
        U = ocv(z, T_hi, p)
        r0 = R0(z, T_hi, p)
        r1 = float(R1(T_hi, p))
        I_hi, D = cur(U, r0 + r1, P_hi)
        Ibar = 0.5 * (I_lo + I_hi)
        Ue_hi = U - r1 * I_hi
    elif limit == "short":
        Ibar = Pbar / U
        for _ in range(60):
            Ue = U - r1 * Ibar
            I_lo, _ = cur(Ue, r0, P_lo)
            I_hi, D = cur(Ue, r0, P_hi)
            Ibar = 0.5 * (I_lo + I_hi)
        Ue_hi = U - r1 * Ibar
    else:
        raise ValueError(limit)
    V_hi = Ue_hi - r0 * I_hi
    g = V_hi - float(p.V_cut)
    if dPb > 0:
        Dk = Ue_hi ** 2 - 4.0 * r0 * (P_hi + dPb)
        Vpk = np.where(Dk > 0, 0.5 * (Ue_hi + np.sqrt(np.maximum(Dk, 0.0))), 0.0)
        g = np.minimum(g, Vpk - float(p.V_sys))
    g = np.where(np.isfinite(g), g, -1.0)
    bad = np.nonzero(g <= 0)[0]
    if bad.size == 0:
        z_end = z[0]
    else:
        i = bad[-1]
        if i + 1 >= z.size:
            return dict(tte_h=0.0, z_end=z0, s_eff_pct=100 * float(cell.display_soc(z0, p)))
        z_end = z[i] + (z[i + 1] - z[i]) * (0 - g[i]) / (g[i + 1] - g[i])
    m = z >= z_end
    zz = np.r_[z_end, z[m]]
    inv = np.r_[np.interp(z_end, z, 1.0 / Ibar), 1.0 / Ibar[m]]
    tte = float(p.Q_max) * float(np.sum(0.5 * (inv[1:] + inv[:-1]) * np.diff(zz)))
    ze = float(cell.z_empty(p))
    return dict(tte_h=tte, z_end=float(z_end), s_eff_pct=100 * (float(z_end) - ze) / (1 - ze), T_C=T - K0)


def burstiness() -> dict:
    c = CFG["burst"]
    out = dict(config=c, cases={})
    for Ta, soh in c["conditions"]:
        p = soh_params(soh)
        Eu = float(cell.usable_energy(p))
        a = float(analytic.a_coefficient(p, Ta)["total"])
        cv2 = ((c["P_hi"] - c["P_lo"]) / (c["P_hi"] + c["P_lo"])) ** 2
        for env in c["envelope_W"]:
            key = f"Tamb={Ta:g}|soh={soh:g}|env={env:g}W"
            t0 = time.time()
            lo = const_load(c["P_lo"], env, "ALT-lo")
            hi = const_load(c["P_hi"], env, "ALT-hi")
            rc = run(Schedule.constant(const_load(c["P_const"], env)), Ta, p=p)
            base = result_dict(rc)
            rows = []
            for per in c["periods_s"]:
                row = dict(period_s=per)
                for order, segs in (("lo_first", [(per / 2, lo), (per / 2, hi)]),
                                    ("hi_first", [(per / 2, hi), (per / 2, lo)])):
                    r = result_dict(run(Schedule(segs), Ta, p=p))
                    r["dtte_pct"] = 100 * (r["tte_h"] / base["tte_h"] - 1)
                    r["dE_pct"] = 100 * (r["E_del_Wh"] / base["E_del_Wh"] - 1)
                    row[order] = r
                row["dtte_mean_pct"] = 0.5 * (row["lo_first"]["dtte_pct"] + row["hi_first"]["dtte_pct"])
                row["dE_mean_pct"] = 0.5 * (row["lo_first"]["dE_pct"] + row["hi_first"]["dE_pct"])
                rows.append(row)
            L1 = {lim: l1_alternation(p, Ta, c["P_lo"], c["P_hi"], env, lim) for lim in ("short", "long",
                                                                                         "long_thermal")}
            L1c = l1_alternation(p, Ta, c["P_const"], c["P_const"], env, "long")
            for lim in L1:
                L1[lim]["dtte_pct"] = 100 * (L1[lim]["tte_h"] / L1c["tte_h"] - 1)
            out["cases"][key] = dict(
                Tamb_C=Ta, soh=soh, envelope_W=env, const=base, L0_tte_h=Eu / c["P_const"],
                L2_const_vs_L0_pct=100 * (base["tte_h"] / (Eu / c["P_const"]) - 1), rows=rows,
                L1_const=L1c, L1_short=L1["short"], L1_long=L1["long"], L1_long_thermal=L1["long_thermal"],
                F10_jensen=dict(a_per_W=a, cv_P2=cv2, dtte_pct=-100 * a * c["P_const"] * cv2,
                                note="F10 first-order L1 penalty of the fluctuation, a·P̄·CV_P² (a at T_amb)"))
            log(f"  burstiness {key}: const {base['tte_h']:.3f} h ({base['event']}); 10 s "
                f"{rows[0]['dtte_mean_pct']:+.2f} %, 3600 s {rows[-1]['dtte_mean_pct']:+.2f} %; L1 short "
                f"{L1['short']['dtte_pct']:+.2f} %, long {L1['long']['dtte_pct']:+.2f} %, long+thermal "
                f"{L1['long_thermal']['dtte_pct']:+.2f} %  [{time.time() - t0:.0f} s]")
    return out


# ------------------------------------------------------------------------------------------- 3 order effect
def l0_schedule_tte(items, p, Tamb_C: float, shift_min: float = 0.0) -> float:
    """Exact L0 TTE [h] of a cyclic schedule (quasi-steady F5 powers at T_amb), starting ``shift_min`` into it."""
    Eu = float(cell.usable_energy(p))
    segs = [(m, float(steady_power(PRIMITIVES[n], p, Tamb_C))) for m, n in items]
    segs = _rotate(segs, shift_min)
    E, t = 0.0, 0.0
    while True:
        for m, P in segs:
            e = P * m / 60.0
            if E + e >= Eu:
                return t + (Eu - E) / P
            E += e
            t += m / 60.0


def _rotate(items, shift_min: float):
    """Rotate a cyclic list of (minutes, x) so that it starts ``shift_min`` minutes into the original."""
    if shift_min <= 0:
        return list(items)
    acc = 0.0
    for i, (m, x) in enumerate(items):
        if acc + m > shift_min:
            first = (acc + m - shift_min, x)
            out = [first] + list(items[i + 1:]) + list(items[:i]) + [(shift_min - acc, x)]
            return [(mm, xx) for mm, xx in out if mm > 1e-9]
        acc += m
    return list(items)


def order_effect() -> dict:
    c = CFG["order"]
    fwd = s12_items()
    rev = list(reversed(fwd))
    out = dict(config=c, cases={})
    for Ta, soh in c["conditions"]:
        p = soh_params(soh)
        key = f"Tamb={Ta:g}|soh={soh:g}"
        t0 = time.time()
        case = {}
        for name, items in (("forward", fwd), ("reversed", rev)):
            scan = []
            for sh in c["phase_shifts_h"]:
                it = [(m, n) for m, n in _rotate(items, 60.0 * sh)]
                r = result_dict(run(Schedule.from_minutes(it), Ta, p=p))
                r["L0_tte_h"] = l0_schedule_tte(items, p, Ta, 60.0 * sh)
                r["shift_h"] = sh
                r["L2_vs_L0_pct"] = 100 * (r["tte_h"] / r["L0_tte_h"] - 1)
                scan.append(r)
            case[name] = dict(**{k: v for k, v in scan[0].items()}, phase_scan=scan)
        f, r = case["forward"], case["reversed"]
        case["reversed_minus_forward"] = dict(
            L2_h=r["tte_h"] - f["tte_h"], L2_pct=100 * (r["tte_h"] / f["tte_h"] - 1),
            L0_h=r["L0_tte_h"] - f["L0_tte_h"], L0_pct=100 * (r["L0_tte_h"] / f["L0_tte_h"] - 1),
            electrochemical_pp=100 * (r["tte_h"] / f["tte_h"] - 1) - 100 * (r["L0_tte_h"] / f["L0_tte_h"] - 1),
            s_off_pp=r["s_off_pct"] - f["s_off_pct"])
        for name in ("forward", "reversed"):
            sc = case[name]["phase_scan"]
            case[name]["phase_range_L2_h"] = [min(x["tte_h"] for x in sc), max(x["tte_h"] for x in sc)]
            case[name]["phase_range_L0_h"] = [min(x["L0_tte_h"] for x in sc), max(x["L0_tte_h"] for x in sc)]
            case[name]["L2_vs_L0_pct_range"] = [min(x["L2_vs_L0_pct"] for x in sc), max(x["L2_vs_L0_pct"] for x in sc)]
        out["cases"][key] = case
        d = case["reversed_minus_forward"]
        log(f"  order {key}: fwd {f['tte_h']:.2f} h ({f['event']}, s_off {f['s_off_pct']:.1f} %), rev {r['tte_h']:.2f} h "
            f"({r['event']}, {r['s_off_pct']:.1f} %); L0 {f['L0_tte_h']:.2f}/{r['L0_tte_h']:.2f}  [{time.time() - t0:.0f} s]")
    return out


# ------------------------------------------------------------------------------------------- 4 initial SOC
def s0_curvature() -> dict:
    c = CFG["s0"]
    u = PRIMITIVES["SO"]
    P = float(demand_power(u, RP63))
    out = dict(config=c, P_W=P, dPb_W=float(u.dPb * RP63.m_b), cases={})
    s = np.array(c["grid"])
    for Ta, soh in c["conditions"]:
        p = soh_params(soh)
        key = f"Tamb={Ta:g}|soh={soh:g}"
        t0 = time.time()
        rr = [run(Schedule.constant(u), Ta, s0=float(v), p=p) for v in s]
        tte = np.array([r.tte_h for r in rr])
        ev = [r.event for r in rr]
        soff = np.array([100 * r.s_off for r in rr])
        t1 = tte[-1]
        def fit(lo):
            m = s >= lo - 1e-9
            A = np.vstack([s[m], np.ones(m.sum())]).T
            coef, *_ = np.linalg.lstsq(A, tte[m], rcond=None)
            resid = tte[m] - A @ coef
            ss_tot = np.sum((tte[m] - tte[m].mean()) ** 2)
            # deviation of every grid point s0 >= 0.2 from this line, in % of TTE(1)
            dev = tte - (coef[0] * s + coef[1])
            md = s >= 0.2 - 1e-9
            return dict(range=[lo, 1.0], slope_h=float(coef[0]), intercept_h=float(coef[1]),
                        s0_at_zero_pct=float(-100 * coef[1] / coef[0]), r2=float(1 - np.sum(resid ** 2) / ss_tot),
                        max_abs_resid_pct_of_tte1=float(100 * np.max(np.abs(resid)) / t1),
                        max_abs_dev_s0ge020_pct_of_tte1=float(100 * np.max(np.abs(dev[md])) / t1))
        i06 = int(np.argmin(np.abs(s - 0.6)))
        i02 = int(np.argmin(np.abs(s - 0.2)))
        # F4 / F2 effective empty point at the steady battery temperature
        T_K = Ta + K0 + float(p.k_env * p.R_th) * P
        zs = float(analytic.brownout_root(P, u.dPb * p.m_b, T_K, p)["z_star"])
        zc = float(analytic.cutoff_soc(P, T_K, p))
        ze = float(cell.z_empty(p))
        s_eff = 100 * (max(zs, zc) - ze) / (1 - ze)
        out["cases"][key] = dict(
            Tamb_C=Ta, soh=soh, s0=list(s), tte_h=list(tte), event=ev, s_off_pct=list(soff),
            ratio_50=float(tte[s == 0.5][0] / t1), ratio_20=float(tte[i02] / t1),
            prop_dev_pp=list(100 * (tte / t1 - s)),
            lin_fit_020=fit(0.2), lin_fit_050=fit(0.5),
            n_immediate=int(np.sum(tte <= 1e-9)),
            chord_curvature=float((t1 + tte[i02] - 2 * tte[i06]) / t1),
            F4_s_off_pct=100 * (zs - ze) / (1 - ze), F2_cutoff_s_pct=100 * (zc - ze) / (1 - ze), F_eff_pct=s_eff,
            T_steady_C=T_K - K0, s_off_L2_at_s0_1_pct=float(soff[-1]))
        log(f"  s0 {key}: TTE(1) {t1:.2f} h, ratio50 {out['cases'][key]['ratio_50']:.3f}, ratio20 "
            f"{out['cases'][key]['ratio_20']:.3f}, intercept {out['cases'][key]['lin_fit_050']['s0_at_zero_pct']:.1f} %, "
            f"F4/F2 {s_eff:.1f} %, L2 s_off {soff[-1]:.1f} %  [{time.time() - t0:.0f} s]")
    return out


# ------------------------------------------------------------------------------------------- 5 heavy tails
def heavy_tail() -> dict:
    c = CFG["heavy_tail"]
    out = dict(config=c, full_chain={}, hub_exact={}, L2={})
    for prof in PROFS:
        ct = CTMC(prof)
        P = state_powers(ct)
        m = analytic.ctmc_moments(ct.generator(), P, EU_WH)
        rows = []
        for law, cv in c["laws"]:
            x = UsageSim(ct, law=law, cv=cv).run(c["n_L0"], seed=c["seed_L0"], P=P)[0]
            s = summ(x)
            f9 = f9_full_sigma2(ct, P, cv ** 2)
            rows.append(dict(law=law, cv=cv, mc_L0=dict(mean_h=s["mean_h"], sd_h=s["sd_h"], se_sd_h=s["se_sd_h"],
                                                       se_mean_h=s["se_mean_h"], cv=s["cv"], skew=s["skew"],
                                                       q05_q50_q95_h=s["q05_q50_q95_h"], p_lt_24h=s["p_lt_24h"]),
                             sd_f9_approx_h=f9["sd_clt_h"], sigma2_f9_approx_W2h=f9["sigma2"],
                             sd_mc_over_f9=s["sd_h"] / f9["sd_clt_h"]))
        e = rows[0]["mc_L0"]
        for r in rows:
            r["sd_ratio_to_exp_mc"] = r["mc_L0"]["sd_h"] / e["sd_h"]
            r["sd_ratio_to_exp_f9"] = r["sd_f9_approx_h"] / rows[0]["sd_f9_approx_h"]
            r["mean_shift_vs_exp_h"] = r["mc_L0"]["mean_h"] - e["mean_h"]
        out["full_chain"][prof] = dict(rows=rows, mu_dynkin_exp_h=m["mu_dynkin"], sd_clt_exp_h=m["sd_clt"])
        log(f"  heavy tail {prof}: " + ", ".join(f"{r['law'][:3]}{r['cv']:g}: sd {r['mc_L0']['sd_h']:.2f} "
                                                 f"(F9 {r['sd_f9_approx_h']:.2f})" for r in rows))
        # exact hub check (Wi-Fi only, F9 is exact asymptotically): physical E_u and a 50x budget
        lam, pj, dj, P0, Pj = hub_inputs(ct, P)
        Pw = np.repeat(np.r_[P0, Pj], 2)
        Pbar = (P0 / lam + np.sum(pj * dj * Pj)) / (1 / lam + np.sum(pj * dj))
        hub = []
        for cv in c["hub_cvs"]:
            s2 = analytic.semi_markov_sigma2(lam, pj, dj, P0, Pj, cv ** 2, 1.0)
            for mult in c["hub_energy_mult"]:
                Eu = EU_WH * mult
                x = UsageSim(ct, law="exp" if cv == 1.0 else "lognormal", cv=cv, connectivity=False).run(
                    c["n_hub"], seed=c["seed_hub"], P=Pw, Eu=Eu, t_max_h=1e7)[0]
                s = summ(x)
                sd_pred = math.sqrt(s2 * Eu / Pbar ** 3)
                hub.append(dict(cv=cv, energy_mult=mult, sigma2_f9_W2h=s2, sd_f9_h=sd_pred, mc_sd_h=s["sd_h"],
                                mc_se_sd_h=s["se_sd_h"], ratio=s["sd_h"] / sd_pred, mc_mean_h=s["mean_h"],
                                fluid_h=Eu / Pbar))
        out["hub_exact"][prof] = hub
        # L2 Monte Carlo: exponential vs lognormal CV 1.5 on coupled paths (common random numbers)
        res = {}
        for law, cv in (("exp", 1.0), ("lognormal", c["cv_L2"])):
            t0 = time.time()
            pa = UsageSim(ct, law=law, cv=cv).run(c["n_L2"], seed=c["seed_L2"], store=True, t_max_h=c["t_max_h"])
            b = run_batch(ct, CFG["Tamb_C"], None, paths=pa, h=c["h_s"], t_max_h=c["t_max_h"])
            l0 = l0_from_paths(pa, P)
            res[law] = dict(tte=b["tte_h"], l0=l0, ev=b["event"])
            log(f"    L2 {prof} {law}: mean {np.nanmean(b['tte_h']):.2f} sd {np.nanstd(b['tte_h']):.2f}  "
                f"[{time.time() - t0:.0f} s]")
        e2, l2 = res["exp"], res["lognormal"]
        se = summ(e2["tte"])
        sl = summ(l2["tte"])
        d = l2["tte"] - e2["tte"]
        out["L2"][prof] = dict(
            exp=dict(mean_h=se["mean_h"], sd_h=se["sd_h"], se_sd_h=se["se_sd_h"], q05_q50_q95_h=se["q05_q50_q95_h"],
                     p_lt_24h=se["p_lt_24h"], brownout_share=float(np.mean(e2["ev"] == "brownout")),
                     L0_same_paths_mean_h=float(np.nanmean(e2["l0"])), L0_same_paths_sd_h=float(np.nanstd(e2["l0"], ddof=1))),
            lognormal=dict(cv=c["cv_L2"], mean_h=sl["mean_h"], sd_h=sl["sd_h"], se_sd_h=sl["se_sd_h"],
                           q05_q50_q95_h=sl["q05_q50_q95_h"], p_lt_24h=sl["p_lt_24h"],
                           brownout_share=float(np.mean(l2["ev"] == "brownout")),
                           L0_same_paths_mean_h=float(np.nanmean(l2["l0"])),
                           L0_same_paths_sd_h=float(np.nanstd(l2["l0"], ddof=1))),
            sd_ratio_L2=sl["sd_h"] / se["sd_h"],
            sd_ratio_L0_same_paths=float(np.nanstd(l2["l0"], ddof=1) / np.nanstd(e2["l0"], ddof=1)),
            paired_mean_diff_h=float(np.nanmean(d)), paired_mean_diff_se_h=float(np.nanstd(d, ddof=1) / math.sqrt(d.size)),
            L2_over_L0_exp=float(np.nanmean(e2["tte"]) / np.nanmean(e2["l0"])),
            L2_over_L0_lognormal=float(np.nanmean(l2["tte"]) / np.nanmean(l2["l0"])))
    return out


# ------------------------------------------------------------------------------------------- 6 jitter
def jitter() -> dict:
    c = CFG["jitter"]
    mult = jitter_bank(c["K"], c["lo"], c["hi"], c["seed_bank"])
    mref = jitter_bank(c["K_ref"], c["lo"], c["hi"], c["seed_bank"] + 1)
    out = dict(config=c, profiles={})
    for prof in PROFS:
        ct = CTMC(prof)
        P = state_powers(ct)
        L = ct.generator()
        base = analytic.ctmc_moments(L, P, EU_WH)
        B = jitter_power_bank(ct, mult)  # (16, K)
        Bref = jitter_power_bank(ct, mref)
        mx, vx = B.mean(1), B.var(1)
        mj = analytic.ctmc_moments(L, mx, EU_WH)
        D = np.repeat(np.r_[1.0 / ct.lam_on, [ct.d_min[a] / 60.0 for a in SESSIONS]], 2)  # activity sojourn [h]
        s2_add = float(2.0 * np.sum(mj["pi"] * D * vx))
        s2 = mj["sigma2"] + s2_add
        sd = math.sqrt(s2 * EU_WH / mj["Pbar"] ** 3)
        analytic_ = dict(Pbar_W=mj["Pbar"], Pbar_base_W=base["Pbar"], mu_dynkin_h=mj["mu_dynkin"],
                         mu_dynkin_base_h=base["mu_dynkin"], sigma2_mean_part_W2h=mj["sigma2"],
                         sigma2_visit_part_W2h=s2_add, sd_clt_h=sd, sd_clt_base_h=base["sd_clt"],
                         sd_ratio=sd / base["sd_clt"], mean_ratio=mj["mu_dynkin"] / base["mu_dynkin"],
                         note="F7 with per-state mean powers plus the per-visit term 2 Σ_x π_x D_a(x) Var_visit(P_x) "
                              "(connectivity assumed fixed within a visit)")
        state_cv = {f"{a}|{cc}": float(math.sqrt(v) / m) for (a, cc), m, v in zip(ct.states, mx, vx)}
        bank_check = dict(mean_rel_diff_max=float(np.max(np.abs(mx / Bref.mean(1) - 1))),
                          sd_rel_diff_max=float(np.max(np.abs(np.sqrt(vx) / np.sqrt(Bref.var(1)) - 1))),
                          K=c["K"], K_ref=c["K_ref"])
        # L0 Monte Carlo, paired: identical activity paths with and without jitter (separate jitter stream)
        x0 = UsageSim(ct).run(c["n_L0"], seed=c["seed_L0"], P=P)[0]
        xj = UsageSim(ct, jitter_K=c["K"]).run(c["n_L0"], seed=c["seed_L0"], P=B.ravel())[0]
        s0, sj = summ(x0), summ(xj)
        dl0 = xj - x0
        # L2 Monte Carlo, paired
        t0 = time.time()
        pc = UsageSim(ct).run(c["n_L2"], seed=c["seed_L2"], store=True, t_max_h=c["t_max_h"])
        pj = UsageSim(ct, jitter_K=c["K"]).run(c["n_L2"], seed=c["seed_L2"], store=True, t_max_h=c["t_max_h"])
        same = bool(np.array_equal(pc.states, pj.states // c["K"]) and np.array_equal(pc.times, pj.times))
        bc = run_batch(ct, CFG["Tamb_C"], None, paths=pc, h=c["h_s"], t_max_h=c["t_max_h"])
        ej = ExpandedCTMC(ct, jitter_expanded_inputs(ct, mult))
        bj = run_batch(ej, CFG["Tamb_C"], None, paths=pj, h=c["h_s"], t_max_h=c["t_max_h"])
        sc2, sj2 = summ(bc["tte_h"]), summ(bj["tte_h"])
        d2 = bj["tte_h"] - bc["tte_h"]
        out["profiles"][prof] = dict(
            analytic=analytic_, state_power_cv=state_cv, bank_check=bank_check,
            mc_L0=dict(base=dict(mean_h=s0["mean_h"], sd_h=s0["sd_h"], se_sd_h=s0["se_sd_h"]),
                       jitter=dict(mean_h=sj["mean_h"], sd_h=sj["sd_h"], se_sd_h=sj["se_sd_h"],
                                   q05_q50_q95_h=sj["q05_q50_q95_h"]),
                       sd_ratio=sj["sd_h"] / s0["sd_h"], paired_mean_diff_h=float(dl0.mean()),
                       paired_mean_diff_se_h=float(dl0.std(ddof=1) / math.sqrt(dl0.size)),
                       sd_jitter_vs_analytic=sj["sd_h"] / sd),
            mc_L2=dict(base=dict(mean_h=sc2["mean_h"], sd_h=sc2["sd_h"], se_sd_h=sc2["se_sd_h"],
                                 brownout_share=float(np.mean(bc["event"] == "brownout"))),
                       jitter=dict(mean_h=sj2["mean_h"], sd_h=sj2["sd_h"], se_sd_h=sj2["se_sd_h"],
                                   brownout_share=float(np.mean(bj["event"] == "brownout")),
                                   q05_q50_q95_h=sj2["q05_q50_q95_h"]),
                       sd_ratio=sj2["sd_h"] / sc2["sd_h"], paired_mean_diff_h=float(np.nanmean(d2)),
                       paired_mean_diff_se_h=float(np.nanstd(d2, ddof=1) / math.sqrt(d2.size)),
                       paired_sd_of_diff_h=float(np.nanstd(d2, ddof=1)), same_activity_paths=same))
        log(f"  jitter {prof}: analytic sd x{analytic_['sd_ratio']:.3f} mean x{analytic_['mean_ratio']:.4f}; L0 sd "
            f"x{sj['sd_h'] / s0['sd_h']:.3f} dmean {dl0.mean():+.3f} h; L2 sd x{sj2['sd_h'] / sc2['sd_h']:.3f} dmean "
            f"{np.nanmean(d2):+.3f} h  [{time.time() - t0:.0f} s]")
    return out


# ------------------------------------------------------------------------------------------- 7 usage parameters
MULTS = ["lam_on_mul"] + [f"d_mul_{a}" for a in SESSIONS] + ["conn_mul_wifi", "conn_mul_cell"]


def _ctmc_from(prof: str, lam_mul=1.0, d_mul=None, conn=(1.0, 1.0), p=None) -> CTMC:
    lam0, p0, d0 = PROFILES[prof]
    d_mul = d_mul or {}
    return CTMC(prof, lam_on=lam0 * lam_mul, p=p if p is not None else p0,
                d_min={a: d0[a] * d_mul.get(a, 1.0) for a in SESSIONS}, mean_w=2.5 * conn[0], mean_c=2.5 * conn[1])


def _moments(ct: CTMC, P) -> dict:
    m = analytic.ctmc_moments(ct.generator(), P, EU_WH)
    pia = m["pi"].reshape(-1, 2).sum(1)
    return dict(Pbar=m["Pbar"], mu_fluid=m["mu_fluid"], mu_dynkin=m["mu_dynkin"], sd=m["sd_clt"], lam_IG=m["lam_IG"],
                screen_on=24 * (1 - pia[0] - pia[ACTS.index("audio")]))


def first_order_bins(x, y, nb: int) -> float:
    """First-order index Var(E[y|x])/Var(y) by quantile binning (bias-corrected for the within-bin noise)."""
    edges = np.quantile(x, np.linspace(0, 1, nb + 1))
    idx = np.clip(np.searchsorted(edges, x, side="right") - 1, 0, nb - 1)
    vy = y.var()
    num = 0.0
    corr = 0.0
    for b in range(nb):
        yy = y[idx == b]
        if yy.size > 1:
            num += yy.size * (yy.mean() - y.mean()) ** 2
            corr += yy.var(ddof=1)
    return float(max(num / y.size - corr / y.size, 0.0) / vy)


def usage_params() -> dict:
    c = CFG["usage_params"]
    out = dict(config=c, profiles={})
    for prof in PROFS:
        t0 = time.time()
        P = state_powers(CTMC(prof))
        lam0, p0, d0 = PROFILES[prof]
        spec = {k: USAGE_UQ_SPEC[k] for k in MULTS}
        spec[f"p_{prof}"] = USAGE_UQ_SPEC[f"p_{prof}"]
        S = sample(spec, c["n"], seed=c["seed"])
        pm = S[f"p_{prof}"]  # (7, n)
        rows = []
        for i in range(c["n"]):
            ct = _ctmc_from(prof, S["lam_on_mul"][i], {a: S[f"d_mul_{a}"][i] for a in SESSIONS},
                            (S["conn_mul_wifi"][i], S["conn_mul_cell"][i]), dict(zip(SESSIONS, pm[:, i])))
            rows.append(_moments(ct, P))
        arr = {k: np.array([r[k] for r in rows]) for k in rows[0]}
        base = _moments(CTMC(prof), P)
        mu, sd = arr["mu_dynkin"], arr["sd"]
        var_path = float(np.mean(sd ** 2))
        var_par = float(np.var(mu))
        # mixture TTE law (IG per usage sample) vs the nominal IG
        tg = np.linspace(1.0, 120.0, 2381)
        cdf = np.mean(stats.invgauss.cdf(tg[None, :], (mu / arr["lam_IG"])[:, None], scale=arr["lam_IG"][:, None]),
                      axis=0)
        mixq = [float(np.interp(q, cdf, tg)) for q in (0.05, 0.5, 0.95)]
        q = lambda x: [float(v) for v in np.percentile(x, [5, 50, 95])]
        # first-order indices on mu_dynkin and sd
        X = {k: S[k] for k in MULTS}
        X.update({f"p_{a}": pm[j] for j, a in enumerate(SESSIONS)})
        S1_mu = {k: first_order_bins(v, mu, c["n_bins"]) for k, v in X.items()}
        S1_sd = {k: first_order_bins(v, sd, c["n_bins"]) for k, v in X.items()}
        # OAT tornado (exact F6/F7 on the full chain) and the F8 closed-form prediction (hub, Wi-Fi only)
        lamh, pj, dj, P0, Pj = hub_inputs(CTMC(prof), P)
        hub = analytic.hub_closed_form(lamh, pj, dj, P0, Pj, SESSIONS)
        oat = {}
        for k in MULTS:
            dist, lo, mode, hi = USAGE_UQ_SPEC[k]
            res = {}
            for tag, v in (("lo", lo), ("hi", hi)):
                if k == "lam_on_mul":
                    ct = _ctmc_from(prof, lam_mul=v)
                    el = hub["elast_lam_on"]
                elif k.startswith("d_mul_"):
                    a = k[6:]
                    ct = _ctmc_from(prof, d_mul={a: v})
                    el = hub["elast_d"][a]
                else:
                    ct = _ctmc_from(prof, conn=(v, 1.0) if k.endswith("wifi") else (1.0, v))
                    el = math.nan
                mm = _moments(ct, P)
                res[tag] = dict(value=v, mu_dynkin_pct=100 * (mm["mu_dynkin"] / base["mu_dynkin"] - 1),
                                mu_fluid_pct=100 * (mm["mu_fluid"] / base["mu_fluid"] - 1),
                                sd_pct=100 * (mm["sd"] / base["sd"] - 1),
                                F8_mu_fluid_pct=100 * (math.exp(el * math.log(v)) - 1) if np.isfinite(el) else math.nan)
            oat[k] = res
        conc = float(USAGE_UQ_SPEC[f"p_{prof}"][1])
        for j, a in enumerate(SESSIONS):
            res = {}
            for tag, qq in zip(("lo", "hi"), c["oat_beta_q"]):
                v = float(stats.beta.ppf(qq, conc * p0[a], conc * (1 - p0[a])))
                pnew = {b: (v if b == a else p0[b] * (1 - v) / (1 - p0[a])) for b in SESSIONS}
                mm = _moments(_ctmc_from(prof, p=pnew), P)
                res[tag] = dict(value=v, mu_dynkin_pct=100 * (mm["mu_dynkin"] / base["mu_dynkin"] - 1),
                                mu_fluid_pct=100 * (mm["mu_fluid"] / base["mu_fluid"] - 1),
                                sd_pct=100 * (mm["sd"] / base["sd"] - 1))
            oat[f"p_{a}"] = res
        out["profiles"][prof] = dict(
            base=dict(mu_dynkin_h=base["mu_dynkin"], sd_clt_h=base["sd"], Pbar_W=base["Pbar"],
                      screen_on_h_per_day=base["screen_on"]),
            mu_dynkin_q05_q50_q95_h=q(mu), sd_clt_q05_q50_q95_h=q(sd), Pbar_q05_q50_q95_W=q(arr["Pbar"]),
            screen_on_q05_q50_q95_h=q(arr["screen_on"]), mu_fluid_q05_q50_q95_h=q(arr["mu_fluid"]),
            mu_dynkin_cv=float(mu.std() / mu.mean()), p_mu_lt_24h=float(np.mean(mu < 24.0)),
            total_variance=dict(path_var_h2=var_path, param_var_h2=var_par, path_share=var_path / (var_path + var_par),
                                param_share=var_par / (var_path + var_par),
                                total_sd_h=math.sqrt(var_path + var_par)),
            mixture_q05_q50_q95_h=mixq, mixture_p_lt_24h=float(np.interp(24.0, tg, cdf)),
            nominal_ig_q05_q50_q95_h=[float(v) for v in analytic.ig_quantiles(base["mu_dynkin"], base["lam_IG"])],
            S1_mu_dynkin=S1_mu, S1_sd_clt=S1_sd, S1_mu_sum=float(sum(S1_mu.values())), oat=oat)
        log(f"  usage params {prof}: mu {q(mu)}, path share {var_path / (var_path + var_par):.3f}; mixture "
            f"{[round(v, 2) for v in mixq]}  [{time.time() - t0:.0f} s]")
    return out


# ------------------------------------------------------------------------------------------- figure
COND_COL = {"Tamb=25|soh=1": BLUE, "Tamb=-10|soh=0.8": ORANGE}
COND_LAB = {"Tamb=25|soh=1": "25 °C, new", "Tamb=-10|soh=0.8": "−10 °C, SOH 0.8"}
PROF_COL = dict(typical=BLUE, heavy=ORANGE)


def fig_usage_fluct(res: dict):
    """2 × 3 panels: (a) burstiness, (b) order effect, (c) initial SOC, (d) heavy-tailed sessions,
    (e) sources of TTE spread (Monte Carlo SD with closed-form markers), (f) usage-model parameter tornado."""
    from matplotlib.lines import Line2D
    apply()
    fig, axs = plt.subplots(2, 3, figsize=(W1, 4.35), gridspec_kw=dict(wspace=0.45, hspace=0.6))
    # (a) burstiness at fixed mean power
    ax = axs[0, 0]
    bu = res["burstiness"]["cases"]
    for key, cs in bu.items():
        ck = key.rsplit("|env", 1)[0]
        env = cs["envelope_W"]
        if ck == "Tamb=25|soh=1" and env > 0:
            continue  # identical to the no-envelope curve at 25 °C (the brown-out never binds)
        per = np.array([r["period_s"] for r in cs["rows"]])
        mid = np.array([r["dtte_mean_pct"] for r in cs["rows"]])
        lo = np.array([min(r["lo_first"]["dtte_pct"], r["hi_first"]["dtte_pct"]) for r in cs["rows"]])
        hi = np.array([max(r["lo_first"]["dtte_pct"], r["hi_first"]["dtte_pct"]) for r in cs["rows"]])
        col = COND_COL[ck]
        ls = "-" if env == 0 else "--"
        ax.plot(per, mid, ls=ls, color=col, lw=1.3, marker="o", ms=2.6)
        ax.vlines(per, lo, hi, color=col, lw=0.7, alpha=0.7)
        lt = cs.get("L1_long_thermal", cs["L1_long"])
        ax.plot([6.5], [cs["L1_short"]["dtte_pct"]], marker=">", color=col, ms=4.2, mfc="white")
        ax.plot([5500.0], [lt["dtte_pct"]], marker="<", color=col, ms=4.2, mfc="white")
    ax.text(12, 0.15, "25 °C, new (± bursts)", fontsize=6.3, color=INK2, va="bottom")
    ax.text(12, -3.05, "−10 °C, SOH 0.8", fontsize=6.3, color=INK2, va="top")
    ax.text(12, -6.75, "same, +3 W bursts", fontsize=6.3, color=INK2, va="bottom")
    ax.set_xscale("log")
    ax.set_xlim(5, 7000)
    ax.set_ylim(-10.2, 1.8)
    ax.set_xticks([10, 60, 600, 3600])
    ax.set_xticklabels(["10 s", "1 min", "10 min", "1 h"])
    ax.minorticks_off()
    ax.axhline(0, color=MUTED, lw=0.7)
    ax.set_xlabel("Period of the 0.4 / 2.0 W alternation")
    ax.set_ylabel("ΔTTE vs constant 1.2 W (%)")
    ax.legend(handles=[Line2D([], [], color=INK2, marker=">", mfc="white", ls="none", ms=4,
                              label="L1 quasi-static limits")],
              loc="lower left", fontsize=6.0, handlelength=1.0, borderaxespad=0.1)
    panel_label(ax, "(a)")

    # (b) order effect: S12 forward vs reversed
    ax = axs[0, 1]
    oe = res["order_effect"]["cases"]
    xs = []
    for i, (key, cs) in enumerate(oe.items()):
        col = COND_COL[key]
        for j, name in enumerate(("forward", "reversed")):
            x = i * 2.6 + j
            c = cs[name]
            ax.bar(x, c["tte_h"], width=0.8, color=col, alpha=1.0 if j == 0 else 0.5, edgecolor="none")
            r2 = c["phase_range_L2_h"]
            ax.vlines(x, r2[0], r2[1], color=INK, lw=0.9)
            ax.plot([x - 0.32, x + 0.32], [c["L0_tte_h"]] * 2, color=INK, lw=1.6)
            xs.append((x, "fwd" if j == 0 else "rev"))
        ax.text(i * 2.6 + 0.5, -0.2, COND_LAB[key], transform=ax.get_xaxis_transform(), ha="center", va="top",
                fontsize=6.6, color=INK2)
    ax.set_xticks([x for x, _ in xs])
    ax.set_xticklabels([t for _, t in xs], fontsize=6.8)
    ax.set_ylabel("S12 time to empty (h)")
    ax.legend(handles=[Line2D([], [], color=INK, lw=1.6, label="L0 (energy only)"),
                       Line2D([], [], color=INK, lw=0.9, label="range over start time")],
              loc="upper right", fontsize=6.2, handlelength=1.3, borderaxespad=0.1)
    ax.set_ylim(0, 60)
    ax.grid(axis="x", visible=False)
    panel_label(ax, "(b)")

    # (c) initial SOC: deviation from proportionality
    ax = axs[0, 2]
    s0c = res["s0_curvature"]["cases"]
    shade = {25.0: SEQ_BLUE[3], 0.0: SEQ_BLUE[6], -10.0: SEQ_BLUE[9], -20.0: SEQ_BLUE[12]}
    ax.axhline(0, color=MUTED, lw=0.7)
    for key, cs in s0c.items():
        s = np.array(cs["s0"])
        dev = 100 * (np.array(cs["tte_h"]) / cs["tte_h"][-1] - s)
        ls = "-" if cs["soh"] >= 1 else "--"
        lab = f"{cs['Tamb_C']:g} °C" + ("" if cs["soh"] >= 1 else ", SOH 0.8")
        ax.plot(100 * s, dev, color=shade[cs["Tamb_C"]], ls=ls, lw=1.3, label=lab)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Initial displayed SOC s0 (%)")
    ax.set_ylabel("TTE(s0)/TTE(100 %) − s0 (pp)")
    ax.set_ylim(-36, 2)
    ax.legend(fontsize=5.8, loc="lower right", handlelength=1.8, borderaxespad=0.1, labelspacing=0.15, ncol=2,
              columnspacing=0.8)
    panel_label(ax, "(c)")

    # (d) heavy-tailed sessions: SD vs session CV
    ax = axs[1, 0]
    ht = res["heavy_tail"]
    for prof in PROFS:
        rows = ht["full_chain"][prof]["rows"]
        col = PROF_COL[prof]
        ln = [r for r in rows if r["law"] == "lognormal"]
        cv = np.array([r["cv"] for r in ln])
        ax.plot(cv, [r["sd_f9_approx_h"] for r in ln], color=col, lw=1.0, ls="--")
        ax.errorbar(cv, [r["mc_L0"]["sd_h"] for r in ln], yerr=[2 * r["mc_L0"]["se_sd_h"] for r in ln], fmt="o",
                    mfc="white", color=col, ms=3.6, lw=0.8, mew=1.0, capsize=0)
        e = rows[0]
        ax.plot([1.0], [e["mc_L0"]["sd_h"]], "o", color=col, ms=3.8)
        l2 = ht["L2"][prof]
        ax.plot([1.06], [l2["exp"]["sd_h"]], "s", color=col, ms=3.2)
        ax.plot([l2["lognormal"]["cv"] * 1.06], [l2["lognormal"]["sd_h"]], "s", mfc="white", color=col, ms=3.2,
                mew=1.0)
    ax.text(3.0, 16.4, "typical", color=INK2, fontsize=6.3, ha="right", va="bottom")
    ax.text(3.0, 11.9, "heavy", color=INK2, fontsize=6.3, ha="right", va="bottom")
    ax.set_ylim(3.5, 18.2)
    ax.set_xscale("log")
    ax.set_xticks([0.5, 1, 1.5, 2, 3])
    ax.set_xticklabels(["0.5", "1", "1.5", "2", "3"])
    ax.minorticks_off()
    ax.set_xlabel("CV of session durations")
    ax.set_ylabel("SD of time to empty (h)")
    ax.legend(handles=[Line2D([], [], color=INK2, lw=1.0, ls="--", label="F9 (CLT)"),
                       Line2D([], [], color=INK2, marker="o", mfc="white", ls="none", ms=3.6, label="MC L0, lognormal"),
                       Line2D([], [], color=INK2, marker="o", ls="none", ms=3.6, label="MC L0, exponential"),
                       Line2D([], [], color=INK2, marker="s", ls="none", ms=3.2, label="MC L2 (same fill)")],
              loc="upper left", fontsize=6.0, handlelength=1.4, borderaxespad=0.1, labelspacing=0.2)
    panel_label(ax, "(d)")

    # (e) sources of TTE spread: Monte Carlo SD (bars) and the closed form (ticks)
    ax = axs[1, 1]
    rs = res["rate_scaling"]
    k4 = rs["kappa"].index(4.0)
    cats = ["Markov", "jitter ±30 %", "lognormal 1.5", "sojourns ×4", "usage params"]
    w = 0.38
    for kk, prof in enumerate(PROFS):
        jt = res["jitter"]["profiles"][prof]
        r15 = [r for r in ht["full_chain"][prof]["rows"] if r["law"] == "lognormal" and abs(r["cv"] - 1.5) < 1e-9][0]
        up = res["usage_params"]["profiles"][prof]
        mc = [jt["mc_L0"]["base"]["sd_h"], jt["mc_L0"]["jitter"]["sd_h"], r15["mc_L0"]["sd_h"],
              rs[prof]["mc_L0_sd_h"][k4], math.nan]
        cf = [up["base"]["sd_clt_h"], jt["analytic"]["sd_clt_h"], r15["sd_f9_approx_h"], rs[prof]["sd_clt_h"][k4],
              up["total_variance"]["total_sd_h"]]
        x = np.arange(len(cats)) + (kk - 0.5) * w
        ax.bar(x, [v if np.isfinite(v) else 0 for v in mc], width=w, color=PROF_COL[prof], edgecolor="none",
               label=prof)
        ax.bar(x[-1:], cf[-1:], width=w, color=PROF_COL[prof], edgecolor="none", alpha=0.45)
        ax.plot(np.c_[x - w / 2, x + w / 2].T, np.c_[cf, cf].T, color=INK, lw=1.3)
    ax.set_xticks(np.arange(len(cats)))
    ax.set_xticklabels(cats, fontsize=6.0, rotation=32, ha="right", rotation_mode="anchor")
    ax.set_ylabel("SD of time to empty (h)")
    ax.grid(axis="x", visible=False)
    ax.legend(handles=[Line2D([], [], color=BLUE, lw=4, label="typical (MC L0)"),
                       Line2D([], [], color=ORANGE, lw=4, label="heavy (MC L0)"),
                       Line2D([], [], color=INK, lw=1.3, label="closed form")], loc="upper left", fontsize=6.0,
              handlelength=1.2, borderaxespad=0.1, labelspacing=0.2)
    ax.set_ylim(0, 19)
    panel_label(ax, "(e)")

    # (f) usage-model parameter tornado (typical, exact F6 Dynkin mean)
    ax = axs[1, 2]
    up = res["usage_params"]["profiles"]["typical"]
    oat = up["oat"]
    rng_ = {k: (v["lo"]["mu_dynkin_pct"], v["hi"]["mu_dynkin_pct"]) for k, v in oat.items()}
    top = sorted(rng_, key=lambda k: -abs(rng_[k][1] - rng_[k][0]))[:8][::-1]

    def nice(k):
        if k == "lam_on_mul":
            return r"$\lambda_\mathrm{on}$"
        if k.startswith("d_mul_"):
            return r"$d_\mathrm{" + k[6:] + "}$"
        if k.startswith("p_"):
            return r"$p_\mathrm{" + k[2:] + "}$"
        return k.replace("conn_mul_", "sojourn ")

    for i, k in enumerate(top):
        lo, hi = rng_[k]
        ax.barh(i, lo, height=0.62, color=SEQ_BLUE[3], edgecolor="none")
        ax.barh(i, hi, height=0.62, color=BLUE, edgecolor="none")
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([nice(k) for k in top], fontsize=7)
    ax.axvline(0, color=MUTED, lw=0.7)
    ax.set_xlabel("Δ Dynkin mean, typical day (%)")
    ax.grid(axis="y", visible=False)
    ax.legend(handles=[Line2D([], [], color=SEQ_BLUE[3], lw=4, label="low end"),
                       Line2D([], [], color=BLUE, lw=4, label="high end")], loc="lower right", fontsize=6.0,
              handlelength=1.0, borderaxespad=0.1, labelspacing=0.2)
    panel_label(ax, "(f)")
    return save(fig, "fig_usage_fluct")


def main():
    import sys
    if "--figs-only" in sys.argv:  # redraw from results/usage_fluct.json
        fig_usage_fluct(json.load(open(RESULTS / "usage_fluct.json")))
        return
    t_start = time.time()
    res = {}
    for name, fn in (("rate_scaling", rate_scaling_block), ("burstiness", burstiness), ("order_effect", order_effect),
                     ("s0_curvature", s0_curvature), ("heavy_tail", heavy_tail), ("jitter", jitter),
                     ("usage_params", usage_params)):
        t0 = time.time()
        log(name)
        res[name] = fn()
        res[name]["runtime_s"] = time.time() - t0
    out = dict(meta=meta("run_usage_fluct.py",
                         seed=dict(heavy_tail_L0=CFG["heavy_tail"]["seed_L0"], heavy_tail_L2=CFG["heavy_tail"]["seed_L2"],
                                   jitter_L0=CFG["jitter"]["seed_L0"], jitter_L2=CFG["jitter"]["seed_L2"],
                                   jitter_bank=CFG["jitter"]["seed_bank"], usage_params=CFG["usage_params"]["seed"]),
                         solver="LSODA (deterministic schedules, model.run); rk4_batch h = 10 s (CTMC L2); exact "
                                "event-driven (L0)", h_s=10.0, Tamb_C=CFG["Tamb_C"]),
               config=CFG, **res)
    out["runtime_s"] = time.time() - t_start
    p = write_json("usage_fluct.json", out)
    log(f"wrote {p} ({out['runtime_s']:.0f} s)")
    fig_usage_fluct(json.load(open(p)))


if __name__ == "__main__":
    main()
