"""AN13: aging scenarios A1-A7 over three years and their effect on time-to-empty (FINAL_SPEC §2.12, §11, §14).

Writes
* results/aging.json
  - ``calibration.{k_c, k_H, k_L, q_cal_1000, mean_g, soh_check_1000, proto, variants}`` (protocol C0, §11.1;
    variants N_80 = 500 / 800 FCE with the A1 years to SOH 0.8)
  - ``g_stress`` (E16 at z = 0.1 ... 1.0 vs [PROTO])
  - ``"{A}".{label, t_yr[], soh[], r_ratio[], q_cal[], q_H[], q_L[], years_to_80, soh_1_2_3yr, loss_split_2yr,
    tte_S13_h[], tte_S12_h[], tte_S4_0C_h[], tte_S4_m10C_h[], s_off_0C_pct[], s_off_m10C_pct[], event_S4_0C[],
    event_S4_m10C[], winter_onset_yr, proto, dev}`` for A1, A2, A3, A4, A5, A6a, A6b, A6c, A7 (monthly grid,
    t = 0 ... 3 yr).  tte_S13_h is the Dynkin mean of the typical CTMC (F6, L0: capacity channel only); the S12 and
    S4 values are full L2 LSODA runs on the aged pack (q_loss, r_growth = γ_R q_loss), S4 cold-soaked (T0 = T_amb).
  - ``soh_grid`` (SOH 1.00 ... 0.60): S13 Dynkin mean, S12 L2 at 25 °C, S4 L2 at 25/0/-10/-20 °C (TTE, s_off,
    event) and the -10 °C / 25 °C channel split (capacity only: r_growth = 0; resistance only: q_loss = 0)
  - ``winter_threshold``: SOH at which the S4 shutdown SOC reaches 5 % displayed (illustrative threshold)
  - ``F11_check``: L2 TTE ratio vs the closed form (1 - q)(1 - a γ_R q P) at 25 °C (§7 F11)
  - ``beta_L_band``: A1 and A6a/b/c with β_L at the ends of U(0.5, 2.64) (recalibrated), as §11.3 requires
  - ``uq``: joint Sobol' sample (N = 128, seed 2026) of the aging factors of §4.5 (k_c, E_cal, g_hi, κ_HV, w_L, E_H,
    E_L, β_L, γ_R) and ΔT_ch, each sample recalibrated with protocol C0; SOH quantiles, years to 0.8 and the
    paired SOH differences against A1 at two years
* paper/figures/fig_aging.pdf (+ png preview): (a) SOH vs years per strategy, (b) TTE relative to new vs SOH,
  (c) shutdown SOC of the S4 social load vs SOH in the cold
* paper/tables/tab_aging.tex

Run:  cd /home/user/MCM && PYTHONPATH=src:analysis python3 analysis/run_aging.py [--workers 2] [--n-uq 128]
      (--figure-only redraws fig_aging and tab_aging from results/aging.json; full run ~5 min on 2 workers)
"""
from __future__ import annotations

import argparse
import math
import time
from dataclasses import replace
from multiprocessing import Pool

import matplotlib.pyplot as plt
import numpy as np

from basics_common import meta, write_json, write_tex
from phonebatt import RP63, CTMC, aged, run
from phonebatt import aging as G
from phonebatt import analytic as A
from phonebatt.params import UQ_SPEC, sample
from phonebatt.power import demand_power
from phonebatt.usage import PRIMITIVES, Schedule, s12_schedule
from style import CAT, INK, INK2, MUTED, SEQ_BLUE, W1, apply, panel_label, save

SEED = 2026
YEARS = 3
MONTHS = np.arange(0, 12 * YEARS + 1)
DAYS = np.round(MONTHS * 365 / 12).astype(int)  # state after this many days (0 = new)
SOH_GRID = np.round(np.arange(1.00, 0.595, -0.01), 2)
T_S4 = (25.0, 0.0, -10.0, -20.0)
WINTER_SOFF_PCT = 5.0  # illustrative "noticeable winter shutdown" threshold (displayed SOC at shutdown)

ORDER = ("A1", "A2", "A3", "A4", "A5", "A6a", "A6b", "A6c", "A7")
LABELS = {
    "A1": "charge to 100 %, held overnight",
    "A2": "80 % charge limit",
    "A3": "optimised charging (hold 80 %, finish 1 h before unplug)",
    "A4": "hot climate, T_amb = 35 °C",
    "A5": "2 h/day at 45 °C (car dashboard), otherwise A1",
    "A6a": "fast charge: CC 1.5C to 80 %",
    "A6b": "charging in a 0 °C room at 1C (not JEITA-compliant)",
    "A6c": "charging in a 0 °C room at 0.5C (JEITA rule)",
    "A7": "1 FCE/day (heavy user, ΔSOC 0.95)",
}
SHORT = {"A1": "A1 100 % nightly", "A2": "A2 80 % limit", "A3": "A3 optimised", "A4": "A4 35 °C climate",
         "A5": "A5 car, 45 °C", "A6a": "A6a fast 1.5C", "A6b": "A6b 0 °C, 1C",
         "A6c": "A6c 0 °C, 0.5C", "A7": "A7 1 FCE/day"}
# §11.3 [PROTO]: SOH at 1/2/3 years, years to SOH 0.8 (None = beyond 3 years)
PROTO = {"A1": ((0.883, 0.834, 0.795), 2.86), "A2": ((0.897, 0.853, 0.818), None),
         "A3": ((0.889, 0.841, 0.804), None), "A4": ((0.796, 0.709, 0.641), 0.96),
         "A5": ((0.841, 0.774, 0.721), 1.58), "A6a": ((0.858, 0.796, 0.747), 1.93),
         "A6b": ((0.849, 0.783, 0.730), 1.72), "A6c": ((0.896, 0.853, 0.819), None),
         "A7": ((0.872, 0.820, 0.779), 2.46)}
PROTO_CAL = dict(k_c=1.413e-3, k_H=2.697e-3, k_L=1.288e-3, q_cal_1000=0.0782, mean_g=1.749)
PROTO_G = (0.702, 1.026, 1.724, 1.749, 2.043, 2.423)  # g(0.1, 0.5, 0.7, 0.8, 0.9, 1.0)
UQ_FACTORS = ("k_c", "E_cal", "g_hi", "k_HV", "w_L", "E_H", "E_L", "beta_L", "gamma_R", "dT_ch")


# ------------------------------------------------------------------------------------------- helpers
def calibrated(p=RP63, **kw):
    c = G.calibrate_details(p, **kw)
    return replace(p, k_H=c["k_H"], k_L=c["k_L"]), c


def years_series(sid, p):
    return G.run_years(years=YEARS, p=p, **G.AGING_SCENARIOS[sid])


def at_days(r, key):
    """Value of a daily series after DAYS days (index 0 = new pack)."""
    new = {"soh": 1.0, "q_cal": 0.0, "q_H": 0.0, "q_L": 0.0, "r_ratio": 1.0}[key]
    full = np.concatenate([[new], np.asarray(r[key], float)])
    return full[DAYS]


def _fast_job(job):
    """One L2 LSODA run on an aged pack.  job = (kind, Tamb_C, q_loss, r_growth)."""
    kind, Ta, q, rg = job
    p = aged(RP63, q, rg)
    sch = s12_schedule() if kind == "S12" else Schedule.constant("SO")
    r = run(sch, Ta, p=p)
    return dict(tte_h=r.tte_h, event=r.event, s_off_pct=100.0 * r.s_off)


def _uq_job(args):
    """Recalibrate (protocol C0) one aging-factor sample and run all aging scenarios for three years."""
    i, theta = args
    p = replace(RP63, **theta)
    p, c = calibrated(p)
    out = {}
    for sid in ORDER:
        r = years_series(sid, p)
        out[sid] = dict(soh=[float(r["soh"][365 * y - 1]) for y in (1, 2, 3)], years_to_80=float(r["years_to_80"]),
                        r_ratio_3yr=float(r["r_ratio"][-1]))
    return i, dict(k_H=c["k_H"], k_L=c["k_L"], res=out)


def qtl(x, qs=(5, 50, 95)):
    x = np.asarray(x, float)
    return [float(v) for v in np.percentile(x, qs)]


def first_crossing(x, y, level):
    """First x where y (sampled on x) reaches ``level`` from below, by linear interpolation (None if never)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    idx = np.nonzero(y >= level)[0]
    if idx.size == 0:
        return None
    k = idx[0]
    if k == 0:
        return float(x[0])
    return float(x[k - 1] + (level - y[k - 1]) * (x[k] - x[k - 1]) / (y[k] - y[k - 1]))


def time_to_soh(t_yr, soh, level):
    """First time (yr) at which the daily SOH series falls to ``level`` (None if not within the horizon)."""
    if level is None:
        return None
    if level >= 1.0 - 1e-12:  # threshold already exceeded by a new pack
        return 0.0
    return first_crossing(t_yr, -np.asarray(soh), -level)


# ------------------------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--n-uq", type=int, default=128)
    ap.add_argument("--figure-only", action="store_true", help="redraw the figure and table from results/aging.json")
    a = ap.parse_args(argv)
    t_start = time.time()
    if a.figure_only:
        import json
        from basics_common import RESULTS
        o = json.loads((RESULTS / "aging.json").read_text())
        pcal = replace(RP63, k_H=o["calibration"]["k_H"], k_L=o["calibration"]["k_L"])
        make_table(o)
        make_figure(o, {sid: years_series(sid, pcal) for sid in ORDER})
        return o
    out: dict = {}

    # ---- calibration (C0, §11.1)
    pcal, c = calibrated(RP63)
    st = G.AgingState()
    prof = G.day_profile("full", 25.0, 1.0, pcal, cal_mode=True)
    for _ in range(1000):
        st = G.advance_day(st, prof, pcal)
    cal = dict(k_c=float(pcal.k_c), k_H=c["k_H"], k_L=c["k_L"], q_cal_1000=c["q_cal_N"], mean_g=c["mean_g"],
               soh_check_1000=st.soh, proto=PROTO_CAL,
               dev_rel={k: c[k2] / PROTO_CAL[k] - 1 for k, k2 in (("k_H", "k_H"), ("k_L", "k_L"),
                                                                     ("q_cal_1000", "q_cal_N"),
                                                                     ("mean_g", "mean_g"))},
               note="mean_g is the RMS of g(z) over C0 (q_cal/(k_c sqrt(1000 d)))")
    variants = {}
    for n80 in (500, 800):
        pv, cv = calibrated(RP63, N_80=n80)
        r = years_series("A1", pv)
        variants[f"N80_{n80}"] = dict(k_H=cv["k_H"], k_L=cv["k_L"], q_cal_N=cv["q_cal_N"],
                                      A1_soh_1_2_3yr=[float(r["soh"][365 * y - 1]) for y in (1, 2, 3)],
                                      A1_years_to_80=float(r["years_to_80"]))
    cal["variants"] = variants
    out["calibration"] = cal
    zg = np.array([0.1, 0.5, 0.7, 0.8, 0.9, 1.0])
    out["g_stress"] = dict(z=zg, g=G.g_stress(zg, RP63), proto=PROTO_G)
    print(f"calibration: k_H={c['k_H']:.4e} k_L={c['k_L']:.4e} q_cal(1000 d)={c['q_cal_N']:.4f} "
          f"SOH(1000)={st.soh:.4f}")

    # ---- scenarios A1-A7 (quasi-static daily driver, §11.2)
    series = {sid: years_series(sid, pcal) for sid in ORDER}
    ct = CTMC("typical")

    # ---- fast-model jobs: per-scenario monthly ages + SOH grid
    jobs = set()
    q_month = {sid: 1.0 - at_days(series[sid], "soh") for sid in ORDER}
    for sid in ORDER:
        for q in q_month[sid]:
            q = round(float(q), 10)
            jobs |= {("S12", 25.0, q, None), ("S4", 0.0, q, None), ("S4", -10.0, q, None)}
    for s in SOH_GRID:
        q = round(1.0 - float(s), 10)
        jobs.add(("S12", 25.0, q, None))
        for T in T_S4:
            jobs.add(("S4", T, q, None))
        for T in (25.0, -10.0):
            jobs.add(("S4", T, q, 0.0))  # capacity channel only
            jobs.add(("S4", T, 0.0, round(float(RP63.gamma_R) * q, 10)))  # resistance channel only
    jobs = sorted(jobs, key=lambda j: (j[0], j[1], j[2], -1 if j[3] is None else j[3]))
    t0 = time.time()
    with Pool(a.workers) as pool:
        res = pool.map(_fast_job, jobs, chunksize=8)
    fast = dict(zip(jobs, res))
    print(f"{len(jobs)} L2 runs in {time.time() - t0:.1f} s")

    def F(kind, T, q, rg=None):
        return fast[(kind, T, round(float(q), 10), None if rg is None else round(float(rg), 10))]

    def dynkin(q):
        return float(A.ctmc_mean_tte(ct, aged(RP63, float(q))))

    # SOH grid (figure panels b, c)
    grid = dict(soh=SOH_GRID, tte_S13_h=[dynkin(1 - s) for s in SOH_GRID],
                tte_S12_h=[F("S12", 25.0, 1 - s)["tte_h"] for s in SOH_GRID], S4={})
    for T in T_S4:
        rows = [F("S4", T, 1 - s) for s in SOH_GRID]
        grid["S4"][f"{T:g}"] = dict(tte_h=[r["tte_h"] for r in rows], s_off_pct=[r["s_off_pct"] for r in rows],
                                    event=[r["event"] for r in rows])
    grid["channels"] = {}
    for T in (25.0, -10.0):
        both = np.array(grid["S4"][f"{T:g}"]["tte_h"])
        qonly = np.array([F("S4", T, 1 - s, 0.0)["tte_h"] for s in SOH_GRID])
        ronly = np.array([F("S4", T, 0.0, float(RP63.gamma_R) * (1 - s))["tte_h"] for s in SOH_GRID])
        grid["channels"][f"S4|{T:g}"] = dict(
            both_rel=both / both[0], capacity_only_rel=qonly / qonly[0], resistance_only_rel=ronly / ronly[0],
            s_off_capacity_only_pct=[F("S4", T, 1 - s, 0.0)["s_off_pct"] for s in SOH_GRID],
            s_off_resistance_only_pct=[F("S4", T, 0.0, float(RP63.gamma_R) * (1 - s))["s_off_pct"]
                                       for s in SOH_GRID])
    grid["note"] = ("aged(RP63, q_loss = 1 - SOH) with r_growth = gamma_R q_loss; S4 = constant social feed (SO) "
                    "cold-soaked (T0 = T_amb); S12 from 07:00 at 25 °C; S13 = Dynkin mean of the typical CTMC (L0)")
    out["soh_grid"] = grid

    # winter-shutdown threshold (illustrative): SOH at which s_off(S4) reaches WINTER_SOFF_PCT
    thr = {}
    for T in (0.0, -10.0, -20.0):
        soff = grid["S4"][f"{T:g}"]["s_off_pct"]
        thr[f"{T:g}"] = first_crossing(-SOH_GRID, soff, WINTER_SOFF_PCT)
        thr[f"{T:g}"] = None if thr[f"{T:g}"] is None else -thr[f"{T:g}"]
    out["winter_threshold"] = dict(s_off_level_pct=WINTER_SOFF_PCT, soh_at_level=thr,
                                   note="illustrative threshold: SOH at which the S4 shutdown happens at >= 5 % "
                                        "displayed (linear interpolation on the SOH grid)")

    # F11 closed form at 25 °C
    acoef = float(A.a_coefficient(RP63)["total"])
    P_S4 = float(demand_power(PRIMITIVES["SO"], RP63))
    l2 = np.array(grid["S4"]["25"]["tte_h"]) / grid["S4"]["25"]["tte_h"][0]
    qg = 1 - SOH_GRID
    cf = (1 - qg) * (1 - acoef * float(RP63.gamma_R) * qg * P_S4)
    out["F11_check"] = dict(a_per_W=acoef, P_W=P_S4, q=qg, l2_ratio=l2, closed_form=cf,
                            max_abs_diff=float(np.max(np.abs(l2 - cf))),
                            at_soh_0p8=dict(l2=float(l2[20]), closed_form=float(cf[20])))

    # per-scenario outputs
    for sid in ORDER:
        r = series[sid]
        qm = q_month[sid]
        soh_y = [float(r["soh"][365 * y - 1]) for y in (1, 2, 3)]
        y80 = float(r["years_to_80"])
        e = dict(label=LABELS[sid], options=G.AGING_SCENARIOS[sid], t_yr=MONTHS / 12.0,
                 soh=at_days(r, "soh"), r_ratio=at_days(r, "r_ratio"), q_cal=at_days(r, "q_cal"),
                 q_H=at_days(r, "q_H"), q_L=at_days(r, "q_L"), years_to_80=y80 if math.isfinite(y80) else None,
                 years_to_80_note="None = SOH stays above 0.80 for the 3-year horizon",
                 soh_1_2_3yr=soh_y,
                 loss_split_2yr=dict(q_cal=float(r["q_cal"][729]), q_H=float(r["q_H"][729]),
                                     q_L=float(r["q_L"][729])))
        e["tte_S13_h"] = [dynkin(q) for q in qm]
        e["tte_S12_h"] = [F("S12", 25.0, q)["tte_h"] for q in qm]
        for T, tag in ((0.0, "0C"), (-10.0, "m10C")):
            rows = [F("S4", T, q) for q in qm]
            e[f"tte_S4_{tag}_h"] = [x["tte_h"] for x in rows]
            e[f"s_off_{tag}_pct"] = [x["s_off_pct"] for x in rows]
            e[f"event_S4_{tag}"] = [x["event"] for x in rows]
        e["winter_onset_yr"] = {T: time_to_soh(r["t_yr"], r["soh"], thr[T]) for T in ("0", "-10", "-20")}
        ps, py = PROTO[sid]
        e["proto"] = dict(soh_1_2_3yr=ps, years_to_80=py)
        e["dev"] = dict(soh_1_2_3yr=[s - q for s, q in zip(soh_y, ps)],
                        years_to_80=None if py is None or not math.isfinite(y80) else y80 - py,
                        years_to_80_both_beyond_3yr=(py is None and not math.isfinite(y80)))
        out[sid] = e
        print(f"{sid:4s} SOH 1/2/3 yr = {soh_y[0]:.3f}/{soh_y[1]:.3f}/{soh_y[2]:.3f} (proto {ps}); "
              f"yrs to 0.8 {y80:.2f} (proto {py}); s_off -10C at 3 yr {e['s_off_m10C_pct'][-1]:.1f} %")

    # A2-A7 vs A1 at two years (pp of SOH)
    out["diff_vs_A1_2yr_pp"] = {sid: 100 * (out[sid]["soh_1_2_3yr"][1] - out["A1"]["soh_1_2_3yr"][1])
                                for sid in ORDER if sid != "A1"}
    out["proto_diff_2yr_pp"] = dict(A2=1.9, A3=0.7, A4=-12.5, A6c_vs_A6b=7.0)
    out["diff_A6c_vs_A6b_2yr_pp"] = 100 * (out["A6c"]["soh_1_2_3yr"][1] - out["A6b"]["soh_1_2_3yr"][1])

    # ---- β_L band for the charging-stress scenarios (recalibrated at each end)
    band = {}
    for bl in (0.5, 2.64):
        pb, cb = calibrated(replace(RP63, beta_L=bl))
        for sid in ("A1", "A6a", "A6b", "A6c"):
            r = years_series(sid, pb)
            y80 = float(r["years_to_80"])
            band.setdefault(sid, {})[f"beta_L={bl:g}"] = dict(
                soh_1_2_3yr=[float(r["soh"][365 * y - 1]) for y in (1, 2, 3)],
                years_to_80=y80 if math.isfinite(y80) else None, k_H=cb["k_H"], k_L=cb["k_L"])
    out["beta_L_band"] = band

    # ---- joint UQ of the aging factors (each sample recalibrated with C0)
    th = sample(list(UQ_FACTORS), a.n_uq, method="sobol", seed=SEED)
    tasks = [(i, {k: float(th[k][i]) for k in UQ_FACTORS}) for i in range(a.n_uq)]
    t0 = time.time()
    with Pool(a.workers) as pool:
        uq_res = dict(pool.map(_uq_job, tasks, chunksize=4))
    print(f"aging UQ: {a.n_uq} samples x {len(ORDER)} scenarios in {time.time() - t0:.1f} s")
    uq = dict(n=a.n_uq, method="scrambled Sobol', seed 2026", factors={k: UQ_SPEC[k] for k in UQ_FACTORS})
    uq["k_H_q05_q50_q95"] = qtl([uq_res[i]["k_H"] for i in range(a.n_uq)])
    uq["k_L_q05_q50_q95"] = qtl([uq_res[i]["k_L"] for i in range(a.n_uq)])
    for sid in ORDER:
        S = np.array([uq_res[i]["res"][sid]["soh"] for i in range(a.n_uq)])
        Y = np.array([uq_res[i]["res"][sid]["years_to_80"] for i in range(a.n_uq)])
        Rr = np.array([uq_res[i]["res"][sid]["r_ratio_3yr"] for i in range(a.n_uq)])
        Yf = np.where(np.isfinite(Y), Y, 1e9)  # beyond the horizon: ranked last, reported as None
        yq = np.percentile(Yf, (5, 50, 95), method="lower")
        yq = np.where(yq > YEARS + 1e-9, np.inf, yq)
        uq[sid] = dict(soh_1yr_q05_q50_q95=qtl(S[:, 0]), soh_2yr_q05_q50_q95=qtl(S[:, 1]),
                       soh_3yr_q05_q50_q95=qtl(S[:, 2]),
                       years_to_80_q05_q50_q95=[float(v) if math.isfinite(v) else None for v in yq],
                       p_reach_80_within_3yr=float(np.mean(np.isfinite(Y))),
                       r_ratio_3yr_q05_q50_q95=qtl(Rr))
    S2 = {sid: np.array([uq_res[i]["res"][sid]["soh"][1] for i in range(a.n_uq)]) for sid in ORDER}
    uq["diff_vs_A1_2yr_pp"] = {sid: dict(q05_q50_q95=qtl(100 * (S2[sid] - S2["A1"])),
                                         p_positive=float(np.mean(S2[sid] > S2["A1"])))
                               for sid in ORDER if sid != "A1"}
    d = 100 * (S2["A6c"] - S2["A6b"])
    uq["diff_A6c_vs_A6b_2yr_pp"] = dict(q05_q50_q95=qtl(d), p_positive=float(np.mean(d > 0)))
    uq["note"] = ("years_to_80 quantiles: None = beyond the 3-year horizon. Calibration fixes SOH(1000 FCE, C0) = 0.8 "
                  "in every sample, so the spread reflects conditions that differ from C0 (temperature, SOC, rate).")
    out["uq"] = uq

    # ---- write JSON
    runtime = time.time() - t_start
    out = dict(meta=meta("run_aging.py", seed=SEED, solver="quasi-static daily driver (aging.run_years, exact for "
                                                           "piecewise-constant segments); LSODA (model.run) for "
                                                           "L2 TTE", h_s=None, runtime_s=round(runtime, 1),
                         workers=a.workers),
               **out)
    write_json("aging.json", out)
    make_table(out)
    make_figure(out, series)
    print(f"done in {time.time() - t_start:.1f} s")
    return out


# ------------------------------------------------------------------------------------------- table
def make_table(o):
    def yr(x):
        return "$>3$" if x is None else f"{x:.2f}"

    lines = [r"\begin{tabular}{@{}llcccc@{}}", r"\toprule",
             r"ID & Habit & SOH at 1 / 2 / 3 yr & Yr to 0.8 & SOH 2 yr, 90\,\% & $s_\mathrm{off}$, 2 yr \\",
             r"\midrule"]
    tex = {"A1": r"100\,\%, held overnight", "A2": r"80\,\% charge limit", "A3": r"optimised charging",
           "A4": r"35\,$^\circ$C climate", "A5": r"car, 2\,h/day at 45\,$^\circ$C", "A6a": r"fast charge, 1.5C",
           "A6b": r"charged at 0\,$^\circ$C, 1C", "A6c": r"charged at 0\,$^\circ$C, 0.5C",
           "A7": r"1 FCE/day"}
    for sid in ORDER:
        e = o[sid]
        s = e["soh_1_2_3yr"]
        u = o["uq"][sid]["soh_2yr_q05_q50_q95"]
        lines.append(f"{sid} & {tex[sid]} & {s[0]:.3f} / {s[1]:.3f} / {s[2]:.3f} & {yr(e['years_to_80'])} & "
                     f"{u[0]:.3f}--{u[2]:.3f} & {e['s_off_m10C_pct'][24]:.1f}\\,\\% \\\\")
    lines += [r"\bottomrule", r"\end{tabular}", ""]
    write_tex("tab_aging.tex", "\n".join(lines))


# ------------------------------------------------------------------------------------------- figure
def _repel(y, gap):
    """Spread label positions so neighbours are at least ``gap`` apart (keeps the mean of each cluster)."""
    y = np.asarray(y, float)
    order = np.argsort(y)
    ys = y[order].copy()
    for _ in range(200):
        moved = False
        for i in range(1, ys.size):
            if ys[i] - ys[i - 1] < gap:
                mid = 0.5 * (ys[i] + ys[i - 1])
                ys[i - 1], ys[i] = mid - gap / 2, mid + gap / 2
                moved = True
        if not moved:
            break
    out = np.empty_like(ys)
    out[order] = ys
    return out


def _deg(T):
    return f"{T:g}".replace("-", "\u2212") + " °C"


def make_figure(o, series):
    apply()
    fig = plt.figure(figsize=(W1, 2.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.6, 1, 1], wspace=0.36)
    ax1, ax2, ax3 = (fig.add_subplot(gs[0, i]) for i in range(3))

    # (a) SOH vs years, direct end labels inside the axes
    colors = {"A1": INK}
    colors.update({sid: CAT[i] for i, sid in enumerate(ORDER[1:])})
    ends = []
    for sid in ORDER:
        r = series[sid]
        t = np.concatenate([[0.0], r["t_yr"]])
        s = np.concatenate([[1.0], r["soh"]])
        ax1.plot(t, s, color=colors[sid], lw=1.6 if sid == "A1" else 1.25, zorder=3 if sid == "A1" else 2)
        ends.append(s[-1])
    ax1.plot([0, 3], [0.8, 0.8], color=MUTED, lw=0.8, ls=(0, (3, 2)), zorder=1)
    ax1.text(0.04, 0.802, "SOH 0.8", color=INK2, fontsize=6.5, va="bottom")
    ly = _repel(ends, 0.0205)
    for sid, y_end, y_lab in zip(ORDER, ends, ly):
        ax1.annotate(SHORT[sid], xy=(3.0, y_end), xytext=(3.18, y_lab), textcoords="data", fontsize=6.2,
                     color=INK, va="center", ha="left",
                     arrowprops=dict(arrowstyle="-", color=colors[sid], lw=0.7, shrinkA=0, shrinkB=0))
    ax1.set_xlim(0, 4.75)
    ax1.spines["bottom"].set_bounds(0, 3)
    ax1.set_ylim(0.62, 1.005)
    ax1.set_xticks([0, 1, 2, 3])
    ax1.grid(axis="x", visible=False)
    ax1.set_xlabel("Years of use", x=3 / 4.75 / 2)
    ax1.set_ylabel("State of health")
    ax1.set_title("SOH by charging habit", loc="left")
    panel_label(ax1, "(a)")

    # (b) TTE relative to new vs SOH
    g = o["soh_grid"]
    soh = np.asarray(g["soh"])
    tcol = {25.0: SEQ_BLUE[4], 0.0: SEQ_BLUE[7], -10.0: SEQ_BLUE[9], -20.0: SEQ_BLUE[12]}
    d = np.asarray(g["tte_S13_h"])
    ax2.plot(soh, 100 * d / d[0], color=INK, lw=1.1, ls=(0, (4, 2)), label="S13 mean (L0)")
    for T in T_S4:
        tte = np.asarray(g["S4"][f"{T:g}"]["tte_h"])
        ax2.plot(soh, 100 * tte / tte[0], color=tcol[T], lw=1.6, label=f"S4 {_deg(T)}")
    ax2.set_xlim(1.0, 0.6)
    ax2.set_xticks([1.0, 0.9, 0.8, 0.7, 0.6])
    ax2.set_ylim(30, 101)
    ax2.set_xlabel("State of health")
    ax2.set_ylabel("TTE vs new pack (%)")
    ax2.set_title("Runtime vs health", loc="left")
    ax2.legend(loc="lower left", handlelength=1.6, borderaxespad=0.1, labelspacing=0.25, fontsize=6.5)
    panel_label(ax2, "(b)")

    # (c) s_off vs SOH (cold); A1 ages marked on the -10 °C curve
    for T in (0.0, -10.0, -20.0):
        ax3.plot(soh, g["S4"][f"{T:g}"]["s_off_pct"], color=tcol[T], lw=1.6, label=f"S4 {_deg(T)}")
    ax3.plot([1.0, 0.6], [WINTER_SOFF_PCT] * 2, color=MUTED, lw=0.8, ls=(0, (3, 2)))
    ax3.text(0.605, WINTER_SOFF_PCT - 0.7, "5 %", color=INK2, fontsize=6.3, va="top", ha="right")
    soff10 = np.asarray(g["S4"]["-10"]["s_off_pct"])
    for k, s in enumerate(o["A1"]["soh_1_2_3yr"], 1):
        y = float(np.interp(-s, -soh, soff10))
        ax3.plot([s], [y], marker="o", ms=3.2, color=INK, zorder=5)
        ax3.annotate(("A1, " if k == 1 else "") + f"{k} y", xy=(s, y), xytext=(2 if k == 1 else 0, 4.5),
                     textcoords="offset points", fontsize=6.2, color=INK2, ha="right" if k == 1 else "center",
                     va="bottom")
    ax3.set_xlim(1.0, 0.6)
    ax3.set_xticks([1.0, 0.9, 0.8, 0.7, 0.6])
    top = max(max(g["S4"]["-20"]["s_off_pct"]), 10)
    ax3.set_ylim(0, top * 1.06)
    ax3.set_xlabel("State of health")
    ax3.set_ylabel("Displayed SOC at shutdown (%)")
    ax3.set_title("Stranded charge", loc="left")
    ax3.legend(loc="upper left", handlelength=1.6, borderaxespad=0.1, labelspacing=0.25, fontsize=6.5)
    panel_label(ax3, "(c)")
    save(fig, "fig_aging")


if __name__ == "__main__":
    main()
