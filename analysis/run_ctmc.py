"""AN5: stochastic usage days S13 (typical) and S14 (heavy): closed forms F6–F10 vs Monte Carlo (FINAL_SPEC §5.3, §7).

Writes
* results/ctmc.json — per profile (``typical``, ``heavy``): ``Pbar_W, sigma2_W2h, mu_fluid_h, mu_dynkin_h, sd_clt_h,
  ig_q05_q50_q95_h``, ``mc_L0{mean_h, sd_h, q05_q50_q95_h, p_lt_24h, ...}`` (20 000 exact L0 paths, seed 1),
  ``mc_L2{...}`` (2000 CTMC paths, rk4_batch h = 10 s, seed 2026), ``mc_L2_proto400`` (the [PROTO] configuration:
  400 paths, seed 7), ``energy_share{}``, ``component_share{}``, ``variance_share{}`` (hub, Wi-Fi only),
  ``elasticities{}`` (F8 closed form and finite differences on the full 16-state chain), ``checks``, ``proto``;
  plus ``rate_scaling`` (κ ∈ {0.25, 0.5, 1, 2, 4}: F7 √κ law vs L0/L2 Monte Carlo, exponential vs lognormal
  sessions), ``raster_paths`` and ``checks`` (two-state formula, generator scaling, step convergence at κ = 0.25).
* paper/figures/fig_ctmc_tte.pdf — (a) activity raster and (b) displayed-SOC paths of three typical days;
  (c)/(d) L2 TTE histograms with the IG law (F7) and the Dynkin mean (F6).
* paper/figures/fig_cv_scaling.pdf — SD, CV and MC/analytic ratio vs the rate scaling κ.

Run:  cd /home/user/MCM && PYTHONPATH=src:analysis python3 analysis/run_ctmc.py      (about 6 min on one core)
"""
from __future__ import annotations

import math
import time

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy import stats

from basics_common import RESULTS as ROOT_RESULTS, meta, write_json
from phonebatt import RP63, analytic, cell, run_batch
from phonebatt.power import power_parts, steady_power
from phonebatt.usage import ACTS, CTMC, PROFILES
from stochastic_common import (EU_WH, SESSIONS, UsageSim, f9_full_sigma2, hub_inputs, ig_summary, l0_from_paths,
                               scaled_ctmc, state_powers, subset_paths, summ, boot_ci)
from style import (BLUE, CAT, GRID, INK, INK2, MUTED, ORANGE, W1, apply, panel_label, save)

CFG = dict(n_L0=20000, seed_L0=1, n_L2=2000, seed_L2=2026, n_L2_proto=400, seed_L2_proto=7, h_s=10.0,
           t_max_h=120.0, Tamb_C=25.0, kappas=[0.25, 0.5, 1.0, 2.0, 4.0], n_L0_kappa=20000, seed_L0_kappa=1,
           n_L2_kappa=400, seed_L2_kappa=7, cv_LN=1.5, n_LN=20000, seed_LN=303, n_dynkin_check=20000,
           seed_dynkin_check=101, step_check=dict(kappa=0.25, n=200, seed=17, h=[10.0, 5.0]))
PROFS = ("typical", "heavy")
SID = dict(typical="S13", heavy="S14")

# FINAL_SPEC §7 "[PROTO] results for F6–F10" (regression targets, not validation)
PROTO = {
    "typical": dict(Pbar_W=0.4400, screen_on_h_per_day=6.5, mu_fluid_h=40.43, mu_dynkin_h=41.25, mc_L0_mean_h=41.34,
                    sigma2_W2h=0.2492, sd_clt_h=7.21, mc_L0_sd_h=7.04, cv_tte=0.175, tau_c_h=0.234, cv_P=1.659,
                    ig_q=[30.25, 40.59, 54.49], mc_L0_q=[29.58, 41.46, 52.69], ig_p_lt_24h=0.002,
                    mc_L0_p_lt_24h=0.009, mc_L2_mean_h=39.53, mc_L2_sd_h=7.05, mc_L2_q=[28.36, 39.38, 51.34],
                    ratio_L2_L0=0.956, s_off_mean_pct=0.70, brownout_share=0.07, jensen=0.987,
                    energy_share=dict(social=0.279, game=0.282, nav=0.106, video=0.100, audio=0.078, vcall=0.074,
                                      idle=0.058, camera=0.024),
                    component_share={"SoC": 0.374, "display": 0.221, "radio": 0.182, "awake": 0.084,
                                     "suspend": 0.045, "audio/BT": 0.043, "camera": 0.027, "GNSS": 0.019},
                    variance_share=dict(game=0.807, nav=0.102, social=0.032, vcall=0.032, audio=0.021),
                    hub_sigma2_W2h=0.236495,
                    elasticities=dict(d_game=-0.257, d_social=-0.169, d_audio=0.065, lam_on=-0.531),
                    rate_scaling=dict(kappa=[0.25, 1.0, 4.0], mu_dynkin_h=[40.63, 41.25, 43.71],
                                      sd_clt_h=[3.61, 7.21, 14.43])),
    "heavy": dict(Pbar_W=0.8224, screen_on_h_per_day=11.5, mu_fluid_h=21.63, mu_dynkin_h=22.50, mc_L0_mean_h=22.52,
                  sigma2_W2h=0.8388, sd_clt_h=5.18, mc_L0_sd_h=5.11, cv_tte=0.230, tau_c_h=0.362, cv_P=1.310,
                  ig_q=[14.72, 21.85, 32.50], mc_L0_q=[14.36, 22.36, 31.16], ig_p_lt_24h=0.650,
                  mc_L0_p_lt_24h=0.623, mc_L2_mean_h=21.22, mc_L2_sd_h=4.77, mc_L2_q=[14.13, 20.87, 29.94],
                  ratio_L2_L0=0.942, s_off_mean_pct=0.91, brownout_share=0.11,
                  energy_share=dict(game=0.545, social=0.157, video=0.119, nav=0.071, vcall=0.059),
                  component_share={"SoC": 0.502, "display": 0.202, "radio": 0.158},
                  variance_share=dict(game=0.94), hub_sigma2_W2h=0.820905),
}
# activity colours (fixed categorical order following the entity; idle = no fill)
ACT_COL = {a: CAT[i] for i, a in enumerate(SESSIONS)}
COMP_GROUPS = (("SoC", "soc"), ("display", "disp"), ("radio", "net"), ("awake", "awake"), ("suspend", "susp"),
               ("wake-ups", "wake"), ("audio/BT", "aud"), ("camera", "cam"), ("GNSS", "gnss"), ("sensors", "sens"))


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


# ------------------------------------------------------------------------------------------- closed forms
def closed_forms(c: CTMC, P: np.ndarray) -> dict:
    L = c.generator()
    m = analytic.ctmc_moments(L, P, EU_WH)
    pia = m["pi"].reshape(-1, 2).sum(1)
    out = dict(Pbar_W=m["Pbar"], sigma2_W2h=m["sigma2"], mu_fluid_h=m["mu_fluid"], mu_dynkin_h=m["mu_dynkin"],
               sd_clt_h=m["sd_clt"], cv_tte=m["cv_tte"], cv_P=m["cv_P"], var_P_W2=m["var_P"], tau_c_h=m["tau_c"],
               lam_IG_h=m["lam_IG"], Eu_Wh=EU_WH,
               screen_on_h_per_day=24.0 * (1.0 - pia[0] - pia[ACTS.index("audio")]),
               pi_activity={a: float(v) for a, v in zip(ACTS, pia)},
               dynkin_terms=dict(h_x0_Wh=float(m["h"][0]),
                                 E_h_dep_Wh=float(np.sum(m["pi"] * P * m["h"]) / m["Pbar"])))
    ig = ig_summary(m["mu_dynkin"], m["lam_IG"])
    out["ig_q05_q50_q95_h"] = ig["q05_q50_q95_h"]
    out["ig_p_lt_24h"] = ig["p_lt_24h"]
    out["ig_p_lt_16h"] = float(analytic.ig_cdf(16.0, m["mu_dynkin"], m["lam_IG"]))
    out["ig_skew"] = ig["skew"]
    z = stats.norm.ppf([0.05, 0.5, 0.95])
    out["normal_q05_q50_q95_h"] = [float(m["mu_dynkin"] + v * m["sd_clt"]) for v in z]
    return out, m


def energy_shares(c: CTMC, P: np.ndarray, pi: np.ndarray) -> tuple[dict, dict]:
    Pbar = float(pi @ P)
    by_act = {a: float(pi.reshape(-1, 2)[i] @ P.reshape(-1, 2)[i]) / Pbar for i, a in enumerate(ACTS)}
    comp = {g: 0.0 for g, _ in COMP_GROUPS}
    for x, u in enumerate(c.primitives()):
        Pns, dem, _, br = power_parts(u, RP63)
        br = {k: float(v) for k, v in br.items()}
        br["soc"] = float(P[x]) - float(Pns)  # throttled SoC power at the F5 fixed point
        for g, k in COMP_GROUPS:
            comp[g] += pi[x] * br[k] / Pbar
    return by_act, comp


def elasticities(prof: str, P: np.ndarray, hub: dict, delta: float = 0.01) -> dict:
    """F8 closed-form elasticities (hub, Wi-Fi only) and central differences on the full 16-state chain."""
    lam0, p0, d0 = PROFILES[prof]

    def metrics(**kw):
        m = analytic.ctmc_moments(CTMC(prof, **kw).generator(), P, EU_WH)
        return np.array([m["mu_fluid"], m["mu_dynkin"], m["sd_clt"]])

    fd = {}
    up, dn = math.log(1 + delta), math.log(1 - delta)
    for a in SESSIONS:
        hi = metrics(d_min={**d0, a: d0[a] * (1 + delta)})
        lo = metrics(d_min={**d0, a: d0[a] * (1 - delta)})
        fd[f"d_{a}"] = (np.log(hi) - np.log(lo)) / (up - dn)
    for key, kw in (("lam_on", lambda f: dict(lam_on=lam0 * f)), ("mean_wifi", lambda f: dict(mean_w=2.5 * f)),
                    ("mean_cell", lambda f: dict(mean_c=2.5 * f))):
        fd[key] = (np.log(metrics(**kw(1 + delta))) - np.log(metrics(**kw(1 - delta)))) / (up - dn)
    closed = {f"d_{a}": float(v) for a, v in hub["elast_d"].items()}
    closed["lam_on"] = float(hub["elast_lam_on"])
    return dict(closed_form_hub_mu_fluid=closed,
                fd_full_mu_fluid={k: float(v[0]) for k, v in fd.items()},
                fd_full_mu_dynkin={k: float(v[1]) for k, v in fd.items()},
                fd_full_sd_clt={k: float(v[2]) for k, v in fd.items()})


# ------------------------------------------------------------------------------------------- Monte Carlo
def l2_summary(b: dict, P: np.ndarray, L0_mean: float, with_paired: bool = True) -> dict:
    tte, ev, so = b["tte_h"], b["event"], b["s_off"]
    s = summ(tte)
    names, cnt = np.unique(ev, return_counts=True)
    s["event_share"] = {str(k): float(v) / ev.size for k, v in zip(names, cnt)}
    s["brownout_share"] = float(np.mean(ev == "brownout"))
    s["s_off_mean_pct"] = float(100 * np.nanmean(so))
    s["s_off_q05_q50_q95_pct"] = [float(v) for v in 100 * np.nanpercentile(so, [5, 50, 95])]
    s["s_off_mean_brownout_pct"] = float(100 * np.nanmean(so[ev == "brownout"])) if (ev == "brownout").any() \
        else math.nan
    s["T_end_mean_C"] = float(np.nanmean(b["T_end_C"]))
    s["ratio_L2_L0"] = s["mean_h"] / L0_mean
    if with_paired:
        l0 = l0_from_paths(b["paths"], P)
        ok = np.isfinite(l0) & np.isfinite(tte)
        r = tte[ok] / l0[ok]
        pair = np.stack([tte[ok], l0[ok]], axis=1)
        rng = np.random.default_rng(0)
        bs = []
        for _ in range(1000):
            i = rng.integers(0, ok.sum(), ok.sum())
            bs.append(pair[i, 0].mean() / pair[i, 1].mean())
        s["paired_L0_same_paths"] = dict(
            n=int(ok.sum()), L0_mean_h=float(l0[ok].mean()), L0_sd_h=float(l0[ok].std(ddof=1)),
            ratio_of_means=float(tte[ok].mean() / l0[ok].mean()),
            ratio_of_means_ci95=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
            mean_ratio=float(r.mean()), ratio_q05_q50_q95=[float(v) for v in np.percentile(r, [5, 50, 95])],
            corr_L2_L0=float(np.corrcoef(tte[ok], l0[ok])[0, 1]),
            ratio_cutoff=float(np.mean(r[ev[ok] == "cutoff"])) if (ev[ok] == "cutoff").any() else math.nan,
            ratio_brownout=float(np.mean(r[ev[ok] == "brownout"])) if (ev[ok] == "brownout").any() else math.nan)
    return s


def mc_checks(c: CTMC, P: np.ndarray, m: dict, x_l0: np.ndarray, cf: dict) -> dict:
    """Closed forms vs the L0 Monte Carlo: Dynkin mean, CLT SD, IG law, the depletion-state approximation."""
    s = summ(x_l0)
    out = dict(
        dynkin_minus_mc_h=cf["mu_dynkin_h"] - s["mean_h"], dynkin_z=(cf["mu_dynkin_h"] - s["mean_h"]) / s["se_mean_h"],
        fluid_minus_mc_h=cf["mu_fluid_h"] - s["mean_h"],
        sd_mc_over_clt=s["sd_h"] / cf["sd_clt_h"], sd_mc_over_clt_se=s["se_sd_h"] / cf["sd_clt_h"],
        ig_minus_mc_q_h=[a - b for a, b in zip(cf["ig_q05_q50_q95_h"], s["q05_q50_q95_h"])],
        normal_minus_mc_q_h=[a - b for a, b in zip(cf["normal_q05_q50_q95_h"], s["q05_q50_q95_h"])],
        ig_minus_mc_p_lt_24h=cf["ig_p_lt_24h"] - s["p_lt_24h"],
        skew_mc=s["skew"], skew_ig=cf["ig_skew"],
        ks_D_ig=float(stats.kstest(x_l0, lambda t: analytic.ig_cdf(t, cf["mu_dynkin_h"], cf["lam_IG_h"])).statistic),
        ks_D_normal=float(stats.kstest(x_l0, lambda t: stats.norm.cdf(t, cf["mu_dynkin_h"], cf["sd_clt_h"])).statistic),
    )
    # Dynkin identity (exact in L0): E[TTE] = (E_u - h(x0) + E[h(X_TTE)])/P̄.  Check it pathwise with the
    # event-driven simulator, which also returns the depletion state, and test the approximation
    # E[h(X_TTE)] ≈ Σ π_x P_x h_x / P̄ (inspection paradox).
    tte, xdep = UsageSim(c).run(CFG["n_dynkin_check"], seed=CFG["seed_dynkin_check"], P=P)
    hv, Pbar = m["h"], m["Pbar"]
    Y = tte - (EU_WH - hv[0] + hv[xdep]) / Pbar
    emp = np.bincount(xdep, minlength=P.size) / xdep.size
    pred = m["pi"] * P / Pbar
    out["dynkin_identity"] = dict(
        n=int(tte.size), mean_residual_h=float(Y.mean()), se_residual_h=float(Y.std(ddof=1) / math.sqrt(Y.size)),
        E_h_dep_empirical_Wh=float(np.sum(emp * hv)), E_h_dep_formula_Wh=float(np.sum(pred * hv)),
        mu_dynkin_with_empirical_dep_h=float((EU_WH - hv[0] + np.sum(emp * hv)) / Pbar),
        dep_state_tv_distance=float(0.5 * np.abs(emp - pred).sum()),
        dep_activity_empirical={a: float(v) for a, v in zip(ACTS, emp.reshape(-1, 2).sum(1))},
        dep_activity_formula={a: float(v) for a, v in zip(ACTS, pred.reshape(-1, 2).sum(1))},
        sim_mean_h=float(tte.mean()), sim_sd_h=float(tte.std(ddof=1)))
    return out


def rate_scaling(prof: str) -> dict:
    """§10.5 item 1 / F7: all sojourns × κ (π fixed).  Exponential: F7 √κ law vs L0 and L2 Monte Carlo;
    lognormal sessions (CV 1.5, X15): F9-increment approximation vs L0 Monte Carlo."""
    base = CTMC(prof)
    P = state_powers(base)
    L1 = base.generator()
    rows = []
    for kap in CFG["kappas"]:
        t0 = time.time()
        c = scaled_ctmc(prof, kap)
        L = c.generator()
        assert np.allclose(L, L1 / kap, rtol=1e-12, atol=1e-12)
        m = analytic.ctmc_moments(L, P, EU_WH)
        x0 = analytic.mc_L0(L, P, EU_WH, n=CFG["n_L0_kappa"], seed=CFG["seed_L0_kappa"])
        s0 = summ(x0)
        tmax = 120.0 if kap <= 1 else 240.0
        b = run_batch(c, CFG["Tamb_C"], None, n=CFG["n_L2_kappa"], h=CFG["h_s"], seed=CFG["seed_L2_kappa"],
                      t_max_h=tmax)
        s2 = summ(b["tte_h"])
        xl = UsageSim(c, law="lognormal", cv=CFG["cv_LN"]).run(CFG["n_LN"], seed=CFG["seed_LN"], P=P)[0]
        sl = summ(xl)
        f9 = f9_full_sigma2(c, P, CFG["cv_LN"] ** 2)
        rows.append(dict(
            kappa=kap, mu_dynkin_h=m["mu_dynkin"], mu_fluid_h=m["mu_fluid"], sd_clt_h=m["sd_clt"],
            cv_clt=m["sd_clt"] / m["mu_dynkin"],
            mc_L0=dict(mean_h=s0["mean_h"], sd_h=s0["sd_h"], se_mean_h=s0["se_mean_h"], se_sd_h=s0["se_sd_h"],
                       cv=s0["cv"], q05_q50_q95_h=s0["q05_q50_q95_h"]),
            mc_L2=dict(mean_h=s2["mean_h"], sd_h=s2["sd_h"], se_mean_h=s2["se_mean_h"], se_sd_h=s2["se_sd_h"],
                       cv=s2["cv"], n_censored=s2["n_censored"], brownout_share=float(np.mean(b["event"] == "brownout"))),
            lognormal=dict(cv_session=CFG["cv_LN"], sd_f9_approx_h=f9["sd_clt_h"], sigma2_f9_approx_W2h=f9["sigma2"],
                           mc_L0=dict(mean_h=sl["mean_h"], sd_h=sl["sd_h"], se_mean_h=sl["se_mean_h"],
                                      se_sd_h=sl["se_sd_h"], cv=sl["cv"], q05_q50_q95_h=sl["q05_q50_q95_h"])),
            sd_mc_L0_over_clt=s0["sd_h"] / m["sd_clt"], sd_mc_L2_over_clt=s2["sd_h"] / m["sd_clt"],
            sd_mc_LN_over_f9=sl["sd_h"] / f9["sd_clt_h"], mean_mc_L0_minus_dynkin_h=s0["mean_h"] - m["mu_dynkin"],
            mean_mc_L2_over_mc_L0=s2["mean_h"] / s0["mean_h"]))
        log(f"  {prof} kappa={kap}: dynkin {m['mu_dynkin']:.2f} clt {m['sd_clt']:.2f} | L0 {s0['mean_h']:.2f}/"
            f"{s0['sd_h']:.2f} | L2 {s2['mean_h']:.2f}/{s2['sd_h']:.2f} | LN {sl['mean_h']:.2f}/{sl['sd_h']:.2f} "
            f"(F9~{f9['sd_clt_h']:.2f})  [{time.time() - t0:.0f} s]")
    k = np.array(CFG["kappas"])
    sd_mc = np.array([r["mc_L0"]["sd_h"] for r in rows])
    sd_l2 = np.array([r["mc_L2"]["sd_h"] for r in rows])
    sd_ln = np.array([r["lognormal"]["mc_L0"]["sd_h"] for r in rows])
    fit = lambda y: float(np.polyfit(np.log(k), np.log(y), 1)[0])
    return dict(rows=rows, loglog_slope=dict(clt=0.5, mc_L0=fit(sd_mc), mc_L2=fit(sd_l2), mc_LN=fit(sd_ln)),
                sd_ratio_to_kappa1=dict(clt=[math.sqrt(v) for v in k],
                                        mc_L0=list(sd_mc / sd_mc[list(k).index(1.0)]),
                                        mc_L2=list(sd_l2 / sd_l2[list(k).index(1.0)])))


# ------------------------------------------------------------------------------------------- one profile
def analyse(prof: str) -> tuple[dict, dict]:
    log(f"profile {prof}")
    c = CTMC(prof)
    P = state_powers(c)
    cf, m = closed_forms(c, P)
    res = dict(scenario=SID[prof], state_powers_W={f"{a}|{cc}": float(v) for (a, cc), v in zip(c.states, P)},
               state_demand_W={f"{a}|{cc}": float(v) for (a, cc), v in zip(c.states, c.state_powers(RP63, 25.0,
                                                                                                   capped=False))},
               **cf)
    # L0 Monte Carlo (exact, 20 000 paths, seed 1 as in §15.3)
    t0 = time.time()
    x_l0 = analytic.mc_L0(c.generator(), P, EU_WH, n=CFG["n_L0"], seed=CFG["seed_L0"])
    res["mc_L0"] = dict(**summ(x_l0), seed=CFG["seed_L0"])
    res["mc_L0_checks"] = mc_checks(c, P, m, x_l0, cf)
    log(f"  L0 MC done ({time.time() - t0:.1f} s): mean {res['mc_L0']['mean_h']:.3f} sd {res['mc_L0']['sd_h']:.3f}")
    # L2 Monte Carlo: [PROTO] configuration and the larger main sample
    t0 = time.time()
    b4 = run_batch(c, CFG["Tamb_C"], None, n=CFG["n_L2_proto"], h=CFG["h_s"], seed=CFG["seed_L2_proto"],
                   t_max_h=CFG["t_max_h"])
    res["mc_L2_proto400"] = dict(**l2_summary(b4, P, res["mc_L0"]["mean_h"]), seed=CFG["seed_L2_proto"])
    log(f"  L2 400 ({time.time() - t0:.0f} s): mean {res['mc_L2_proto400']['mean_h']:.3f}")
    t0 = time.time()
    b = run_batch(c, CFG["Tamb_C"], None, n=CFG["n_L2"], h=CFG["h_s"], seed=CFG["seed_L2"], t_max_h=CFG["t_max_h"])
    res["mc_L2"] = dict(**l2_summary(b, P, res["mc_L0"]["mean_h"]), seed=CFG["seed_L2"])
    log(f"  L2 {CFG['n_L2']} ({time.time() - t0:.0f} s): mean {res['mc_L2']['mean_h']:.3f} sd {res['mc_L2']['sd_h']:.3f}"
        f" brownout {res['mc_L2']['brownout_share']:.3f}")
    a = analytic.a_coefficient(RP63, CFG["Tamb_C"])["total"]
    res["jensen"] = dict(a_per_W=float(a), factor=float(1 - a * cf["Pbar_W"] * (1 + cf["cv_P"] ** 2)),
                         mc_ratio_L2_L0=res["mc_L2"]["ratio_L2_L0"],
                         paired_ratio=res["mc_L2"]["paired_L0_same_paths"]["ratio_of_means"])
    # shares
    by_act, comp = energy_shares(c, P, m["pi"])
    res["energy_share"] = by_act
    res["component_share"] = comp
    lam, pj, dj, P0, Pj = hub_inputs(c, P)
    hub = analytic.hub_closed_form(lam, pj, dj, P0, Pj, SESSIONS)
    mm = analytic.ctmc_moments(c.activity_generator(), np.r_[P0, Pj], EU_WH)
    res["variance_share"] = {k: float(v) for k, v in hub["share"].items()}
    res["hub"] = dict(sigma2_closed_W2h=hub["sigma2"], sigma2_devmatrix_W2h=mm["sigma2"],
                      abs_diff=abs(hub["sigma2"] - mm["sigma2"]), Pbar_W=hub["Pbar"], EL_h=hub["EL"],
                      pi_idle=hub["pi_idle"], note="activity chain alone, Wi-Fi state powers (F8)")
    # variance attribution with π_c-averaged powers (connectivity folded into the session powers)
    hubA = analytic.hub_closed_form(*hub_inputs(c, P, average_connectivity=True), SESSIONS)
    res["variance_share_conn_avg"] = {k: float(v) for k, v in hubA["share"].items()}
    res["elasticities"] = elasticities(prof, P, hub)
    res["proto"] = PROTO[prof]
    return res, dict(b=b, P=P, x_l0=x_l0, m=m)


def proto_deviation(res: dict) -> dict:
    """Signed deviations from the §7 [PROTO] numbers (ours − PROTO)."""
    out = {}
    for prof in PROFS:
        r, pr = res[prof], PROTO[prof]
        d = {}
        for k in ("Pbar_W", "screen_on_h_per_day", "mu_fluid_h", "mu_dynkin_h", "sigma2_W2h", "sd_clt_h", "cv_tte",
                  "tau_c_h", "cv_P"):
            d[k] = r[k] - pr[k]
        d["ig_q_h"] = [a - b for a, b in zip(r["ig_q05_q50_q95_h"], pr["ig_q"])]
        d["ig_p_lt_24h"] = r["ig_p_lt_24h"] - pr["ig_p_lt_24h"]
        d["mc_L0_mean_h"] = r["mc_L0"]["mean_h"] - pr["mc_L0_mean_h"]
        d["mc_L0_sd_h"] = r["mc_L0"]["sd_h"] - pr["mc_L0_sd_h"]
        d["mc_L0_q_h"] = [a - b for a, b in zip(r["mc_L0"]["q05_q50_q95_h"], pr["mc_L0_q"])]
        d["mc_L0_p_lt_24h"] = r["mc_L0"]["p_lt_24h"] - pr["mc_L0_p_lt_24h"]
        for tag in ("mc_L2_proto400", "mc_L2"):
            q = r[tag]
            d[tag] = dict(mean_pct=100 * (q["mean_h"] / pr["mc_L2_mean_h"] - 1), sd_h=q["sd_h"] - pr["mc_L2_sd_h"],
                          q_h=[a - b for a, b in zip(q["q05_q50_q95_h"], pr["mc_L2_q"])],
                          ratio_L2_L0=q["ratio_L2_L0"] - pr["ratio_L2_L0"],
                          s_off_mean_pp=q["s_off_mean_pct"] - pr["s_off_mean_pct"],
                          brownout_share=q["brownout_share"] - pr["brownout_share"])
        d["energy_share"] = {k: r["energy_share"][k] - v for k, v in pr["energy_share"].items()}
        d["component_share"] = {k: r["component_share"][k] - v for k, v in pr["component_share"].items()}
        d["variance_share"] = {k: r["variance_share"][k] - v for k, v in pr["variance_share"].items()}
        d["hub_sigma2_W2h"] = r["hub"]["sigma2_closed_W2h"] - pr["hub_sigma2_W2h"]
        if "elasticities" in pr:
            el = r["elasticities"]["closed_form_hub_mu_fluid"]
            d["elasticities"] = {k: el[k] - v for k, v in pr["elasticities"].items()}
        if "jensen" in pr:
            d["jensen"] = r["jensen"]["factor"] - pr["jensen"]
        out[prof] = d
    return out


# ------------------------------------------------------------------------------------------- figures
def pick_raster_paths(b: dict, qs=(0.05, 0.5, 0.95)) -> list[int]:
    tte = b["tte_h"]
    out = []
    for q in qs:
        target = np.nanquantile(tte, q)
        out.append(int(np.nanargmin(np.abs(tte - target))))
    return out


def raster_data(prof: str, b: dict) -> dict:
    rows = pick_raster_paths(b)
    sub = subset_paths(b["paths"], rows)
    rb = run_batch(CTMC(prof), CFG["Tamb_C"], None, paths=sub, h=CFG["h_s"], t_max_h=CFG["t_max_h"], record=True)
    ze = float(cell.z_empty(RP63))
    tr = rb["traj"]
    s = (tr["x"][:, 0, :] - ze) / (1 - ze)
    return dict(rows=rows, tte_batch_h=[float(b["tte_h"][r]) for r in rows], tte_rerun_h=list(map(float, rb["tte_h"])),
                event=[str(e) for e in rb["event"]], s_off_pct=[float(100 * v) for v in rb["s_off"]],
                t_h=tr["t"] / 3600.0, soc=s, paths=sub)


def fig_ctmc_tte(res: dict, ras: dict, mc: dict):
    apply()
    fig = plt.figure(figsize=(W1, 3.05))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.42, 1.0], height_ratios=[0.62, 1.0], hspace=0.16, wspace=0.24,
                          left=0.07, right=0.99, top=0.86, bottom=0.13)
    ax_r = fig.add_subplot(gs[0, 0])
    ax_s = fig.add_subplot(gs[1, 0], sharex=ax_r)
    gs2 = gs[:, 1].subgridspec(2, 1, hspace=0.12)
    ax_t = fig.add_subplot(gs2[0])
    ax_h = fig.add_subplot(gs2[1], sharex=ax_t)

    # (a) activity raster of three typical paths (5 %, 50 %, 95 % TTE quantiles)
    pa = ras["paths"]
    labels = ["5 % path", "median", "95 % path"]
    tte = ras["tte_rerun_h"]
    for i in range(3):
        y = 2 - i
        t = pa.times[i] / 3600.0
        st = pa.states[i]
        end = tte[i]
        for j in range(t.size - 1):
            t0, t1 = t[j], min(t[j + 1], end)
            if not np.isfinite(t0) or t0 >= end:
                break
            act = ACTS[st[j] // 2]
            if act != "idle":
                ax_r.broken_barh([(t0, t1 - t0)], (y - 0.32, 0.64), facecolors=ACT_COL[act], edgecolor="none",
                                 linewidth=0)
            if st[j] % 2 == 1:  # cellular connectivity: thin gray strip under the row
                ax_r.broken_barh([(t0, t1 - t0)], (y - 0.47, 0.1), facecolors=MUTED, edgecolor="none", linewidth=0)
        ax_r.plot([end], [y], marker="|", ms=8, mew=1.4, color=INK)
    ax_r.set_yticks([2, 1, 0])
    ax_r.set_yticklabels(labels)
    ax_r.set_ylim(-0.6, 2.6)
    ax_r.grid(False)
    ax_r.spines["left"].set_visible(False)
    ax_r.tick_params(axis="y", length=0)
    plt.setp(ax_r.get_xticklabels(), visible=False)
    handles = [Patch(facecolor=ACT_COL[a], label=a) for a in SESSIONS] + \
              [Patch(facecolor=MUTED, label="on LTE")]
    ax_r.legend(handles=handles, ncol=8, loc="lower left", bbox_to_anchor=(0.0, 1.0), fontsize=6.5,
                handlelength=0.9, handletextpad=0.3, columnspacing=0.6, borderaxespad=0.2)
    panel_label(ax_r, "(a)")

    # (b) displayed SOC (D1) along the same paths, L2 model
    cols = [INK, INK2, MUTED]
    t_h = ras["t_h"]
    for i in range(3):
        keep = t_h <= tte[i]
        ax_s.plot(np.r_[t_h[keep], tte[i]], 100 * np.r_[ras["soc"][keep, i], ras["s_off_pct"][i] / 100],
                  color=cols[i], lw=1.3)
        ax_s.annotate(f"{tte[i]:.1f} h", (tte[i], 100 * ras["s_off_pct"][i] / 100), xytext=(2, 4),
                      textcoords="offset points", fontsize=6.5, color=INK2, ha="left", va="bottom")
    ax_s.set_ylim(-3, 103)
    ax_s.set_xlim(0, max(tte) * 1.09)
    ax_s.set_ylabel("Displayed SOC (%)")
    ax_s.set_xlabel("Time since full charge (h)  —  typical profile S13, 25 °C")
    panel_label(ax_s, "(b)")

    # (c)/(d) TTE histograms (L2 MC) with the L0 MC outline, IG law and Dynkin mean
    bins = np.arange(0, 80.5, 1.0)
    tgrid = np.linspace(1, 80, 600)
    for ax, prof, tag in ((ax_t, "typical", "(c)"), (ax_h, "heavy", "(d)")):
        r = res[prof]
        x2 = mc[prof]["b"]["tte_h"]
        ax.hist(x2[np.isfinite(x2)], bins=bins, density=True, color="#cfcec6", edgecolor="white", linewidth=0.3,
                label=f"L2 Monte Carlo ({CFG['n_L2']} paths)")
        ax.hist(mc[prof]["x_l0"], bins=bins, density=True, histtype="step", color=INK2, linewidth=0.8,
                label=f"L0 Monte Carlo ({CFG['n_L0'] // 1000}k paths)")
        ax.plot(tgrid, analytic.ig_pdf(tgrid, r["mu_dynkin_h"], r["lam_IG_h"]), color=BLUE, lw=1.5,
                label="Inverse Gaussian (F7)")
        ymax = ax.get_ylim()[1]
        ax.axvline(r["mu_dynkin_h"], color=BLUE, lw=1.0, ls="--", label="Dynkin mean (F6)")
        ax.axvline(r["mc_L2"]["mean_h"], color=INK, lw=0.9, ls=":", label="L2 mean")
        ax.set_ylim(0, ymax * 1.18)
        name = "Typical day S13" if prof == "typical" else "Heavy day S14"
        left = prof == "typical"  # put the label where the distribution is not
        ax.text(0.02 if left else 0.98, 0.95, f"{name}\nDynkin {r['mu_dynkin_h']:.1f} h\nL2 mean {r['mc_L2']['mean_h']:.1f} h",
                transform=ax.transAxes, ha="left" if left else "right", va="top", fontsize=6.8, color=INK)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.grid(axis="y", visible=False)
        panel_label(ax, tag)
    plt.setp(ax_t.get_xticklabels(), visible=False)
    ax_t.set_xlim(5, 72)
    ax_h.set_xlabel("Time to empty (h)")
    ax_t.legend(loc="lower left", bbox_to_anchor=(-0.02, 1.0), ncol=2, fontsize=6.3, handlelength=1.4,
                columnspacing=0.8, borderaxespad=0.2)
    return save(fig, "fig_ctmc_tte")


def fig_cv_scaling(rs: dict):
    apply()
    fig, axs = plt.subplots(1, 3, figsize=(W1, 2.35), gridspec_kw=dict(wspace=0.36))
    kk = np.geomspace(0.2, 5.0, 50)
    pc = dict(typical=BLUE, heavy=ORANGE)
    for prof in PROFS:
        rows = rs[prof]["rows"]
        k = np.array([r["kappa"] for r in rows])
        col = pc[prof]
        i1 = list(k).index(1.0)
        sd1 = rows[i1]["sd_clt_h"]
        sdl1 = rows[i1]["lognormal"]["sd_f9_approx_h"]
        mu = np.array([r["mu_dynkin_h"] for r in rows])
        # (a) SD
        ax = axs[0]
        ax.plot(kk, sd1 * np.sqrt(kk), color=col, lw=1.3)
        ax.plot(kk, sdl1 * np.sqrt(kk), color=col, lw=1.0, ls="--")
        ax.errorbar(k, [r["mc_L0"]["sd_h"] for r in rows], yerr=[2 * r["mc_L0"]["se_sd_h"] for r in rows], fmt="o",
                    color=col, ms=3.8, lw=0.8, capsize=0)
        ax.errorbar(k * 1.08, [r["mc_L2"]["sd_h"] for r in rows], yerr=[2 * r["mc_L2"]["se_sd_h"] for r in rows],
                    fmt="s", color=col, ms=3.2, lw=0.8, capsize=0)
        ax.errorbar(k / 1.08, [r["lognormal"]["mc_L0"]["sd_h"] for r in rows],
                    yerr=[2 * r["lognormal"]["mc_L0"]["se_sd_h"] for r in rows], fmt="o", mfc="white", color=col,
                    ms=3.8, lw=0.8, capsize=0, mew=1.0)
        # (b) CV
        ax = axs[1]
        ax.plot(k, np.array([r["sd_clt_h"] for r in rows]) / mu, color=col, lw=1.3)
        ax.plot(k, np.array([r["lognormal"]["sd_f9_approx_h"] for r in rows]) / mu, color=col, lw=1.0, ls="--")
        ax.plot(k, [r["mc_L0"]["cv"] for r in rows], "o", color=col, ms=3.8)
        ax.plot(k * 1.08, [r["mc_L2"]["cv"] for r in rows], "s", color=col, ms=3.2)
        ax.plot(k / 1.08, [r["lognormal"]["mc_L0"]["cv"] for r in rows], "o", mfc="white", color=col, ms=3.8,
                mew=1.0)
        # (c) MC / analytic SD
        ax = axs[2]
        ax.errorbar(k, [r["sd_mc_L0_over_clt"] for r in rows],
                    yerr=[2 * r["mc_L0"]["se_sd_h"] / r["sd_clt_h"] for r in rows], fmt="o-", color=col, ms=3.8,
                    lw=0.9, capsize=0)
        ax.errorbar(k * 1.08, [r["sd_mc_L2_over_clt"] for r in rows],
                    yerr=[2 * r["mc_L2"]["se_sd_h"] / r["sd_clt_h"] for r in rows], fmt="s", color=col, ms=3.2,
                    lw=0.8, capsize=0)
        ax.errorbar(k / 1.08, [r["sd_mc_LN_over_f9"] for r in rows],
                    yerr=[2 * r["lognormal"]["mc_L0"]["se_sd_h"] / r["lognormal"]["sd_f9_approx_h"] for r in rows],
                    fmt="o--", mfc="white", color=col, ms=3.8, lw=0.9, capsize=0, mew=1.0)
    for ax in axs:
        ax.set_xscale("log")
        ax.set_xticks([0.25, 0.5, 1, 2, 4])
        ax.set_xticklabels(["0.25", "0.5", "1", "2", "4"])
        ax.minorticks_off()
        ax.set_xlabel("Sojourn scaling κ")
    axs[0].set_yscale("log")
    axs[0].minorticks_off()
    axs[0].set_yticks([2, 4, 8, 16])
    axs[0].set_yticklabels(["2", "4", "8", "16"])
    axs[0].set_ylabel("SD of time to empty (h)")
    axs[1].set_ylabel("Coefficient of variation")
    axs[2].set_ylabel("Monte Carlo SD / closed form")
    axs[2].axhline(1.0, color=MUTED, lw=0.8)
    axs[0].text(0.3, 0.9, "slope ½", transform=axs[0].transAxes, fontsize=6.8, color=INK2)
    panel_label(axs[0], "(a)")
    panel_label(axs[1], "(b)")
    panel_label(axs[2], "(c)")
    h1 = [Line2D([], [], color=BLUE, lw=1.3, label="Typical S13"), Line2D([], [], color=ORANGE, lw=1.3,
                                                                          label="Heavy S14")]
    h2 = [Line2D([], [], color=INK2, lw=1.3, label="F7 CLT, exp. sessions"),
          Line2D([], [], color=INK2, lw=1.0, ls="--", label="F9, lognormal CV 1.5"),
          Line2D([], [], color=INK2, marker="o", ls="none", ms=3.8, label="MC L0, exp."),
          Line2D([], [], color=INK2, marker="s", ls="none", ms=3.2, label="MC L2, exp."),
          Line2D([], [], color=INK2, marker="o", mfc="white", ls="none", ms=3.8, label="MC L0, lognormal")]
    fig.legend(handles=h1 + h2, loc="lower center", ncol=7, fontsize=6.6, bbox_to_anchor=(0.5, 0.99),
               handlelength=1.5, columnspacing=0.9, handletextpad=0.4)
    return save(fig, "fig_cv_scaling")


# ------------------------------------------------------------------------------------------- table
def write_table(res: dict):
    """paper/tables/tab_ctmc.tex: closed forms (F6, F7) vs Monte Carlo L0 and L2 for S13 and S14."""
    from basics_common import write_tex
    T, H = res["typical"], res["heavy"]
    q = lambda v: "/".join(f"{x:.1f}" for x in v)

    def row(label, f):
        return f"{label} & {f(T)} & {f(H)} \\\\"

    lines = [
        r"\begin{tabular}{lrr}", r"\toprule", r"Quantity & Typical (S13) & Heavy (S14) \\", r"\midrule",
        row(r"Mean power $\bar P$ (W); screen-on (h/day)", lambda r: f"{r['Pbar_W']:.3f}; {r['screen_on_h_per_day']:.1f}"),
        row(r"$\mu_{\mathrm{fluid}}$ / $\mu_{\mathrm{Dynkin}}$ (h)", lambda r: f"{r['mu_fluid_h']:.2f} / {r['mu_dynkin_h']:.2f}"),
        row(r"MC L0 mean $\pm$ SE (h)", lambda r: f"{r['mc_L0']['mean_h']:.2f} $\\pm$ {r['mc_L0']['se_mean_h']:.2f}"),
        row(r"SD: CLT / MC L0 / MC L2 (h)", lambda r: f"{r['sd_clt_h']:.2f} / {r['mc_L0']['sd_h']:.2f} / {r['mc_L2']['sd_h']:.2f}"),
        row(r"5/50/95\,\%: IG (F7)", lambda r: q(r["ig_q05_q50_q95_h"])),
        row(r"5/50/95\,\%: MC L0", lambda r: q(r["mc_L0"]["q05_q50_q95_h"])),
        row(r"5/50/95\,\%: MC L2", lambda r: q(r["mc_L2"]["q05_q50_q95_h"])),
        row(r"$\Pr(\mathrm{TTE}<24$\,h): IG / L0 / L2", lambda r: f"{r['ig_p_lt_24h']:.3f} / {r['mc_L0']['p_lt_24h']:.3f} / {r['mc_L2']['p_lt_24h']:.3f}"),
        row(r"L2/L0 mean ratio (same paths)", lambda r: f"{r['mc_L2']['paired_L0_same_paths']['ratio_of_means']:.3f}"),
        row(r"Brown-out share; mean $s_{\mathrm{off}}$ (\%)", lambda r: f"{100 * r['mc_L2']['brownout_share']:.1f}\\,\\%; {r['mc_L2']['s_off_mean_pct']:.2f}"),
        row(r"Game share of energy / of $\sigma^2$ (hub)", lambda r: f"{r['energy_share']['game']:.2f} / {r['variance_share']['game']:.2f}"),
        r"\bottomrule", r"\end{tabular}", ""]
    return write_tex("tab_ctmc.tex", "\n".join(lines))


# ------------------------------------------------------------------------------------------- main
def _figs_from_cache(path: str):
    """Developer shortcut: redraw the figures from a pickle written by a full run (env PB_FIG_CACHE)."""
    import json
    import pickle
    with open(path, "rb") as fh:
        res, ras, mc = pickle.load(fh)
    rs = json.load(open(ROOT_RESULTS / "ctmc.json"))["rate_scaling"]
    write_table(res)
    fig_ctmc_tte(res, ras, mc)
    fig_cv_scaling(rs)


def main():
    import os
    cache = os.environ.get("PB_FIG_CACHE")
    if cache and os.path.exists(cache) and "--figs-only" in os.sys.argv:
        _figs_from_cache(cache)
        return
    t_start = time.time()
    res, mc = {}, {}
    for prof in PROFS:
        res[prof], mc[prof] = analyse(prof)
    # checks independent of the profile
    a, b_, P1, P0 = 3.0, 0.7, 2.4, 0.05
    m2 = analytic.ctmc_moments(np.array([[-a, a], [b_, -b_]]), np.array([P0, P1]), 10.0)
    checks = dict(two_state_sigma2=dict(devmatrix=m2["sigma2"], closed=analytic.two_state_sigma2(a, b_, P1, P0)))
    sc = CFG["step_check"]
    c = scaled_ctmc("typical", sc["kappa"])
    t0 = time.time()
    r10 = run_batch(c, CFG["Tamb_C"], None, n=sc["n"], h=sc["h"][0], seed=sc["seed"], t_max_h=CFG["t_max_h"])
    r5 = run_batch(c, CFG["Tamb_C"], None, h=sc["h"][1], paths=r10["paths"], t_max_h=CFG["t_max_h"])
    checks["step_convergence_kappa025"] = dict(
        n=sc["n"], mean_h10=float(np.mean(r10["tte_h"])), mean_h5=float(np.mean(r5["tte_h"])),
        rel_diff_pct=float(100 * (np.mean(r10["tte_h"]) / np.mean(r5["tte_h"]) - 1)),
        max_abs_path_diff_h=float(np.max(np.abs(r10["tte_h"] - r5["tte_h"]))),
        same_event_share=float(np.mean(r10["event"] == r5["event"])), tolerance_pct=0.5,
        note="T1b-style check (same continuous-time paths) at the fastest switching used here")
    log(f"step check kappa=0.25: {checks['step_convergence_kappa025']['rel_diff_pct']:+.3f} % "
        f"({time.time() - t0:.0f} s)")
    # rate scaling
    rs = {}
    for prof in PROFS:
        log(f"rate scaling {prof}")
        rs[prof] = rate_scaling(prof)
    rs["kappa"] = CFG["kappas"]
    rs["config"] = {k: CFG[k] for k in ("kappas", "n_L0_kappa", "seed_L0_kappa", "n_L2_kappa", "seed_L2_kappa", "cv_LN",
                                         "n_LN", "seed_LN", "h_s")}
    rs["proto"] = PROTO["typical"]["rate_scaling"]
    # raster paths and figures
    ras = raster_data("typical", mc["typical"]["b"])
    checks["raster_rerun_identical"] = float(np.max(np.abs(np.array(ras["tte_rerun_h"]) - ras["tte_batch_h"])))
    out = dict(meta=meta("run_ctmc.py", seed=dict(L0=CFG["seed_L0"], L2=CFG["seed_L2"], L2_proto=CFG["seed_L2_proto"]),
                         solver="rk4_batch (L2, event-synced, model.run_batch); exact event-driven (L0)",
                         h_s=CFG["h_s"], Tamb_C=CFG["Tamb_C"], profiles=dict(typical="S13", heavy="S14")),
               config=CFG, **res, rate_scaling=rs, checks=checks,
               raster_paths=dict(profile="typical", rows=ras["rows"], tte_h=ras["tte_rerun_h"], event=ras["event"],
                                 s_off_pct=ras["s_off_pct"]),
               proto_deviation=proto_deviation(res))
    out["runtime_s"] = time.time() - t_start
    p = write_json("ctmc.json", out)
    log(f"wrote {p}")
    if cache:
        import pickle
        slim = {k: dict(b=dict(tte_h=v["b"]["tte_h"]), x_l0=v["x_l0"]) for k, v in mc.items()}
        with open(cache, "wb") as fh:
            pickle.dump((res, ras, slim), fh)
    write_table(res)
    fig_ctmc_tte(res, ras, mc)
    fig_cv_scaling(rs)
    log(f"figures written; total {time.time() - t_start:.0f} s")


if __name__ == "__main__":
    main()
