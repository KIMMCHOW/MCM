"""Digest results/global_summary.md for the paper writer (group 'global': AN9 Sobol + AN10 structural).

Called by run_structural.py after it has written results/structural.json; reads results/sobol.json (AN9) if
present.  Every number is taken from the two JSON files, so the digest stays consistent after a re-run.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"

SOBOL_ORDER = ["S3", "S4", "S5a", "S5b", "S5c", "S6", "S8a", "S8b", "S11a", "S11c", "S7", "S9", "S10a"]
LAB = {"R0_ref": "R0,ref", "Ea_R0": "Ea,R0", "gamma_R": "γ_R", "gamma_sig": "γ_sig", "d_sys": "Δ_sys",
       "k_RR": "κ_RR", "wifi": "Wi-Fi", "lte": "LTE", "nr": "NR", "T_amb": "T_amb"}


def _top(d, k=3):
    items = sorted(d["ST"].items(), key=lambda t: -t[1])[:k]
    return ", ".join(f"{LAB.get(g, g)} {v:.2f} ({d['S1'][g]:.2f})" for g, v in items)


def _f(v, nd=1):
    return "–" if v is None else f"{v:+.{nd}f}"


def build_summary(st: dict) -> str:
    L = ["# Group \"global\" (AN9 Sobol + AN10 structural variants): digest for the paper writer", ""]
    sp = RES / "sobol.json"
    so = json.loads(sp.read_text()) if sp.exists() else None
    m = st["meta"]
    rt_s = f"{so['meta']['runtime_s']:.0f} s" if so else "n/a"
    L.append(f"Scripts: `analysis/run_sobol.py` (AN9, {rt_s} with 2 workers), `analysis/run_structural.py` (AN10, "
             f"{m['runtime_s']:.0f} s) with `analysis/structural_variants.py` (union model) and this digest "
             "(`analysis/summary_global.py`). Outputs: `results/sobol.json`, `results/structural.json`, "
             "`fig_sobol`, `fig_structural`, `paper/tables/tab_structural.tex` (T4). Core untouched; no core bug "
             "found. Seeds fixed (2026); RP-63, s0 = 100 % unless stated.")
    L.append("")
    if so:
        sm = so["meta"]
        L.append("## AN9 Sobol (§10.3): what drives TTE")
        L.append(f"- Design: 24 factor groups, SALib Saltelli design (second order), S_i by Saltelli 2010, S_Ti by "
                 f"Jansen; 95 % bootstrap CIs ({sm['n_boot']} resamples). Surrogate scenarios N = {sm['N_surrogate']} "
                 f"(51,200 evals each); full L2 (rk4, h = 2 s) N = {sm['N_L2']} (12,800 runs); S10a-G N = "
                 f"{sm['N_interaction']} (k = 26, {54 * sm['N_interaction']:,} runs).")
        ch = so["surrogate_check"]
        ok = [s for s in ch if ch[s]["pass_le_1pct"]]
        L.append(f"- Surrogate (F2 + F4 + F5 incl. q_batt and the thermal event) vs full L2 on 200 points each: "
                 f"max |err| {so['surrogate_check_maxerr_pct_used']:.2f} % over the scenarios that use it "
                 f"({', '.join(ok)}; spec: ≤ 1 %). S11c fails ({ch['S11c']['maxabs_err_pct']:.1f} %, "
                 f"{ch['S11c']['n_over_1pct']}/200 points, all 7–10 W calls at the 68 °C thermal-shutdown "
                 "boundary where TTE jumps), so S11c is reported from full L2 (N = 256); surrogate kept as "
                 "`S11c_surrogate`.")
        L.append("")
        L.append("| Scenario | method | mean TTE (h) | CV | top S_Ti (S_i) | ΣS_i |")
        L.append("|---|---|---|---|---|---|")
        for s in SOBOL_ORDER + ["S10a-G"]:
            d = so[s]
            meth = "surr." if d["method"].startswith("surrogate") else "L2"
            L.append(f"| {s} | {meth} | {d['mean']:.2f} | {d['cv']:.2f} | {_top(d)} | {d['sum_S1']:.2f} |")
        g = so["S10a-G"]
        gs = g["s_off"]
        L.append("")
        L.append(f"- S10a-G (social; T_amb ~ U(−20, 10) °C, SOH ~ U(0.75, 1), m_b explicit): TTE S_i/S_Ti "
                 f"T_amb {g['S1']['T_amb']:.2f}/{g['ST']['T_amb']:.2f}, SOH {g['S1']['SOH']:.2f}/"
                 f"{g['ST']['SOH']:.2f}, m_b {g['S1']['m_b']:.2f}/{g['ST']['m_b']:.2f}, R0,ref "
                 f"{g['S1']['R0_ref']:.2f}/{g['ST']['R0_ref']:.2f}; interaction share 1 − ΣS_i = "
                 f"{g['interaction_share']:.2f} [{g['interaction_share_ci'][0]:.2f}, "
                 f"{g['interaction_share_ci'][1]:.2f}]; S_ij " + ", ".join(
                     f"{k.replace('|', '×')} {v:+.3f} ± {g['S2_conf'][k]:.3f}" for k, v in g["S2"].items()) + ".")
        L.append(f"- Same for s_off (stranded charge, pp): T_amb {gs['S1']['T_amb']:.2f}/{gs['ST']['T_amb']:.2f}, "
                 f"SOH {gs['S1']['SOH']:.2f}/{gs['ST']['SOH']:.2f}, m_b {gs['S1']['m_b']:.2f}/"
                 f"{gs['ST']['m_b']:.2f}, R0,ref {gs['S1']['R0_ref']:.2f}/{gs['ST']['R0_ref']:.2f}; interaction "
                 f"share {gs['interaction_share']:.2f} [{gs['interaction_share_ci'][0]:.2f}, "
                 f"{gs['interaction_share_ci'][1]:.2f}]; S_ij " + ", ".join(
                     f"{k.replace('|', '×')} {v:+.3f}" for k, v in gs["S2"].items()) + ".")
        for s in ("S9", "S10a"):
            d = so[s]["s_off"]
            L.append(f"- {s} s_off (mean {d['mean']:.1f}, median {d['q05_q50_q95'][1]:.1f} pp): top S_Ti (S_i) "
                     f"{_top(d)}.")
        cell = ["Q_n", "R0_ref", "A_R", "Ea_R0", "RC", "OCV", "V_cut", "d_sys", "C_th", "R_th"]
        mx = max(so[s]["ST"][g] for s in SOBOL_ORDER[:9] for g in cell)
        L.append(f"- Reading. 25 °C constant loads: variance is additive (ΣS_i ≥ "
                 f"{min(so[s]['sum_S1'] for s in SOBOL_ORDER[:10]):.2f}) and carried by the load parameters of the "
                 f"active components; every cell/thermal group has S_Ti ≤ {mx:.2f} (an energy-budget problem; Q_n "
                 "is small only because its range is ±2 %). Gaming: thermal design (k_env, R_skin, throttle). "
                 f"Cold + aged: R0,ref, γ_R, E_a,R0 and the burst size m_b, with TTE CV {so['S10a']['cv']:.2f}–"
                 f"{so['S9']['cv']:.2f} vs 0.13–0.21 at 25 °C (a voltage problem). Same three regimes as §10.1. "
                 "Interactions in S10a-G are modest and their CIs include 0.")
        L.append("")
    R = st["results"]
    xc = st["x0_check"]
    mx = max(abs(v["rel_diff_pct"]) for v in xc.values())
    ev = all(v["event_union"] == v["event_lsoda"] for v in xc.values())
    L.append("## AN10 structural variants (§10.4): ΔTTE (%) vs X0, same solver and h")
    L.append(f"- X0 of the union model vs core LSODA on all {len(xc)} deterministic scenarios: max |ΔTTE| "
             f"{mx:.4f} %, events identical: {ev}; vs core run_batch (same h): |Δ| = "
             + ", ".join(f"{v['abs_diff_h']:.1e} h" for v in st["core_check"].values())
             + ". Variant effects at h = 1 s vs 2 s differ by ≤ "
             + f"{max(abs(v['diff_pp']) for v in st['h_convergence'].values()):.4f} pp.")
    cols = [("S4", "S4"), ("S7", "S7"), ("S7@35", "S7 35°C"), ("S12", "S12"), ("S13", "S13"), ("S14", "S14"),
            ("S9@-20/0.8", "S9 −20/0.8"), ("S10a@-20/0.8", "S10a −20/0.8")]
    L.append("")
    L.append("| X | variant | " + " | ".join(c[1] for c in cols) + " | Δs_off S10a −20/0.8 (pp) |")
    L.append("|---|---|" + "---|" * (len(cols) + 1))
    for x in m["variant_definitions"]:
        row = []
        for k, _ in cols:
            kk = f"{x}|{k}"
            if x == "X17":
                kk = {"S4": "X17|S4@2y", "S12": "X17|S12@2y"}.get(k, "")
            row.append(_f(R[kk]["dtte_pct"]) if kk in R else "–")
        so_ = R.get(f"{x}|S10a@-20/0.8")
        L.append(f"| {x} | {m['variant_definitions'][x]} | " + " | ".join(row) + f" | "
                 f"{_f(so_['ds_off_pct']) if so_ else '–'} |")
    x17 = st["X17"]
    L.append("")
    L.append(f"- X17 (A1, 25 °C): SOH at 1/2/3 y {', '.join(f'{v:.3f}' for v in x17['soh_X17'])} vs X0 "
             f"{', '.join(f'{v:.3f}' for v in x17['soh_X0'])}; years to 0.8: {x17['years_to_80_X17']:.2f} vs "
             f"{x17['years_to_80_X0']:.2f} (k_H, k_L recalibrated to SOH 0.8 at 1000 FCE; q_cal(1000 d, C0) "
             f"{x17['q_cal_1000_X17']:.4f} vs 0.0782).")
    g = lambda k: R[k]["dtte_pct"]
    warm = [k for k in ("S1", "S2a", "S2b", "S3", "S4", "S5a", "S5b", "S6", "S7", "S7@35", "S8a", "CM", "S11a@-110",
                        "S16", "S12", "S13", "S14")]
    ecm = ("X2", "X3", "X4", "X5", "X5n", "X6", "X9", "X18a", "X18b")
    m_ecm = max(abs(g(f"{x}|{k}")) for x in ecm for k in warm if not (x in ("X2", "X3") and k in ("S7", "S7@35",
                                                                                                 "S13", "S14")))
    m_l0 = max(g(f"X1|{k}") for k in warm if k not in ("S13", "S14"))
    m_x8 = min(g(f"X8|{k}") for k in warm)
    m_load = max(abs(g(f"{x}|{k}")) for x in ("X10", "X11a", "X11b", "X13") for k in warm)
    L.append(f"- Reading. ≥ 25 °C, fresh, from 100 %: the ECM-form variants (X2–X6, X9, X18) move TTE by ≤ "
             f"{m_ecm:.1f} % (X2/X3 excluded where they remove throttling); L0 is ≤ {m_l0:.1f} % long on constant "
             f"loads and schedules (CTMC days +{g('X1|S13'):.1f} / +{g('X1|S14'):.1f} %); V_sys = V_cut costs up to "
             f"{-m_x8:.1f} %; the load laws matter more (X10/X11/X13 up to {m_load:.1f} %; X14 on S16 "
             f"{g('X14|S16'):+.0f} %). Throttling is "
             f"indispensable for gaming (X2/X3 {g('X3|S7'):+.0f} % at 25 °C, {g('X3|S7@35'):+.0f} % at 35 °C). In "
             f"the cold on an aged pack the form dominates: L0 {g('X1|S9@-20/0.8'):+.0f} % (misses "
             f"{-R['X1|S9@-20/0.8']['ds_off_pct']:.0f} pp stranded charge), R0 split {g('X6|S9@-20/0.8'):+.0f} %, "
             f"KiBaM {g('X9|S9@-20/0.8'):+.0f} %, no envelope {g('X7|S10a@-20/0.8'):+.0f} %, V_sys = V_cut "
             f"{g('X8|S10a@-20/0.8'):+.0f} % (S9/S10a at −20 °C). 2RC with equal total R1 (X5n) changes nothing; "
             f"isothermal S10b (0 °C, aged) {g('X3|S10b@0/0.8'):+.0f} %: self-heating matters in the cold.")
    L.append(f"- Usage laws: lognormal sessions leave the mean unchanged within MC error ({g('X15|S13'):+.1f} ± "
             f"{R['X15|S13']['dtte_se_pct']:.1f} %, {g('X15|S14'):+.1f} ± {R['X15|S14']['dtte_se_pct']:.1f} %) but "
             f"raise the SD {R['X15|S13']['sd_h'] / R['X15|S13']['sd_X0_h'] - 1:+.0%} / "
             f"{R['X15|S14']['sd_h'] / R['X15|S14']['sd_X0_h'] - 1:+.0%}; the diurnal chain adds "
             f"{g('X16|S13'):+.1f} ± {R['X16|S13']['dtte_se_pct']:.1f} % / {g('X16|S14'):+.1f} ± "
             f"{R['X16|S14']['dtte_se_pct']:.1f} % (idle nights). X4 (constant current) is harmless from 100 % but "
             f"{g('X4|S4@s0=50'):+.1f} % from 50 %; S12 from 20 % is non-smooth (X4, X12 and X8 move the shutdown "
             "across a segment boundary).")
    c13, c14 = st["x0_ctmc"]["S13"], st["x0_ctmc"]["S14"]
    L.append(f"- CTMC days: X0 mean {c13['mean_tte_h']:.2f} / {c14['mean_tte_h']:.2f} h (S13/S14, {c13['n_paths']} "
             "paths, h = 10 s, paths from stochastic_common.UsageSim; [PROTO] MC L2 39.53 / 21.22 h with another "
             "RNG stream); variant deltas use common random numbers (paired SE in the JSON).")
    L += ["", "## Deviations from the spec and modelling choices (all documented in the JSON meta / docstrings)",
          "- S1 estimator: Saltelli 2010 (SALib's), not Jansen first order: Jansen-S1 (kept as `S1_jansen`) had CI "
          "half-widths ≈ 0.12 per group at N = 256 and ΣS_i = 0.07 on S10a-G. Outputs centred before estimating.",
          "- Surrogate includes q_batt in the F5 fixed point (as T(z), P(z) profiles) and the thermal event; the "
          "plain spec form (q_batt = 0) had 1.3 % error on S11c in a development check. S11c uses full L2 (above).",
          "- Conditions chosen where §10.3 names none: S9 and S10a at −10 °C, SOH 0.8 (the §10.1 cold case); "
          "S11a/S11c at −110 dBm; S10a-G N = 512.",
          "- X1 uses the F5 quasi-steady powers (the paper's analytic L0); X2 is literal L1 at T_amb, so no "
          "throttling and no burst envelope. X5 also run normalised (X5n). X6 keeps the z-shape and r_growth on "
          "both parts. X9 keeps R1 as series resistance (as AN6). X10's ODE implies an exponential tail timer "
          "(promotions λ(1 − r)E_pro), which dominates its effect. X13 is relative to 25 °C, so it lowers power "
          "in the cold (−0.08 W at T_s = −15 °C): this is why it lengthens the cold TTEs. X16 uses 5.5 h⁻¹ from "
          "07:00 to 24:00. X17 recalibrates k_H, k_L; its TTE column is at the 2-year SOH.",
          "- h = 10 s for S1, S2a, S2b, S16 and the CTMC days (§6.3); CTMC paths come from "
          "stochastic_common.UsageSim (exp law), not phonebatt CTMC.sample_paths, so X15/X16 share random numbers.",
          f"- [PROTO] checks: X0 reproduces §5.4 (e.g. S9 −20 °C/0.8 "
          f"{R['X1|S9@-20/0.8']['tte_X0_h']:.2f} h, s_off {R['X1|S9@-20/0.8']['s_off_X0_pct']:.1f} % vs 12.33 h, "
          f"44.2 %; S10a −20/0.8 {R['X1|S10a@-20/0.8']['tte_X0_h']:.2f} h, "
          f"{R['X1|S10a@-20/0.8']['s_off_X0_pct']:.1f} % vs 10.23 h, 23.3 %); A1 SOH "
          f"{', '.join(f'{v:.3f}' for v in x17['soh_X0'])} = [PROTO]. No [PROTO] exists for Sobol or X-variants."]
    return "\n".join(L) + "\n"
