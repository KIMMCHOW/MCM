"""AN3: deterministic scenarios S1-S16 of FINAL_SPEC §5.2 (L2 model, LSODA reference solver).

Writes
* results/scenarios.json  - key "{S}|Tamb={T}|soh={SOH}|s0={s0}" -> {tte_h, event, s_off_pct, soe_off_pct,
  t_disp0_h, P_mean_W, T_end_C, Ts_end_C, rho_perf, ...}; plus "S15" (initial-SOC curves and curvature
  metrics), "proto_comparison" (every §5.4 [PROTO] number next to ours) and "meta".
* paper/figures/fig_trajectories.pdf, fig_tte_s0.pdf, fig_throttle.pdf
* paper/tables/tab_baseline.tex

Run:  cd /home/user/MCM && PYTHONPATH=src:analysis python3 analysis/run_scenarios.py
Deterministic runs use model.run (LSODA, rtol 1e-7, §6.2).  The stochastic days S13/S14 use model.run_batch
(rk4_batch, h = 10 s, 400 CTMC paths, seed 7; §6.3) and are summarised by their path statistics.
"""
from __future__ import annotations

import json
import math
import subprocess
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from phonebatt import CTMC, RP63, run_batch, scenario
from phonebatt.cell import z_empty
from phonebatt.power import demand_power, rho_perf
from style import (BLUE, CAT, GRID, INK, INK2, MUTED, ORANGE, RED_ARM, SEQ_BLUE, W1, apply, panel_label,
                   save)

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
TAB = ROOT / "paper" / "tables"
SPEC = "FINAL_SPEC 2026-11-17"
SEED = 7
N_PATHS = 400
H_MC = 10.0

# ------------------------------------------------------------------------------------------- [PROTO] targets
# FINAL_SPEC §5.4 (TTE h, event, s_off %, T_end °C, Ts_end °C); key -> target.  All at 25 °C, fresh, s0 = 1
# unless the key says otherwise.  P is the §5.4 "Mean P" column (demand before throttling; for S7 it is the
# throttled end power, which final_proto.py prints as P_end).
PROTO_MAIN = {
    "S1": (442.7, "cutoff", 0.02, 25.1, 25.2, 0.040), "S2a": (378.3, "cutoff", 0.02, 25.1, 25.3, 0.047),
    "S2b": (173.8, "cutoff", 0.04, 25.3, 25.6, 0.102), "S3": (38.88, "cutoff", 0.18, 26.4, 27.7, 0.456),
    "S4": (16.64, "cutoff", 0.40, 28.4, 31.2, 1.060), "S5a": (32.32, "cutoff", 0.21, 26.7, 28.2, 0.548),
    "S5b": (25.49, "cutoff", 0.27, 27.2, 29.0, 0.694), "S5c": (21.10, "cutoff", 0.32, 27.7, 29.9, 0.837),
    "S6": (11.83, "cutoff", 0.56, 29.8, 33.7, 1.486), "S8a": (6.93, "cutoff", 0.93, 33.1, 39.9, 2.518),
    "S8b": (10.83, "cutoff", 0.61, 30.2, 34.5, 1.622), "CM": (10.40, "brownout", 0.68, 30.4, 34.9, 1.687),
    "S16": (79.20, "cutoff", 0.09, 25.7, 26.3, 0.224),
}
PROTO_S7 = {15: (4.40, "brownout", 1.89, 27.6, 38.2, 3.89), 25: (5.68, "cutoff", 1.07, 34.5, 42.4, 2.94),
            35: (10.29, "cutoff", 0.57, 40.3, 44.7, 1.65)}
PROTO_S7_NOTHROTTLE_25 = 4.43
PROTO_S12 = {1.0: 34.92, 0.5: 11.79, 0.2: 8.70}
PROTO_S12_SOFF = 0.93
PROTO_S12_PDEM = 0.528
PROTO_S9 = {(25, 1.0): (26.99, 1.7), (-10, 1.0): (25.33, 7.2), (-20, 1.0): (22.67, 16.9),
            (25, 0.8): (21.32, 2.7), (-10, 0.8): (17.00, 22.2), (-20, 0.8): (12.33, 44.2)}
PROTO_S9_COLDSOAK = 17.00  # -10 °C, SOH 0.8, T0 = -10 °C
PROTO_S10A = {  # (s_off %, TTE h); brownout except 25 °C / SOH 1.0
    (25, 1.0): (0.4, 16.64), (25, 0.9): (0.5, 14.94), (25, 0.8): (0.9, 13.22),
    (10, 1.0): (0.7, 16.57), (10, 0.9): (1.3, 14.80), (10, 0.8): (2.1, 13.05),
    (0, 1.0): (1.4, 16.43), (0, 0.9): (2.4, 14.62), (0, 0.8): (3.8, 12.80),
    (-10, 1.0): (2.8, 16.18), (-10, 0.9): (4.9, 14.22), (-10, 0.8): (8.5, 12.15),
    (-20, 1.0): (6.4, 15.54), (-20, 0.9): (13.1, 12.99), (-20, 0.8): (23.3, 10.23)}
PROTO_S10B = {(25, 0.8): (29.19, "cutoff", None), (0, 1.0): (34.50, "brownout", 1.9),
              (0, 0.8): (24.78, "brownout", 4.7)}
PROTO_S11A = {-85: 14.46, -95: 13.72, -105: 13.05, -110: 12.73, -120: 12.15}
PROTO_S11C = {-85: 6.93, -95: 6.17, -105: 5.58, -110: 5.33, -120: 4.87}
PROTO_S6_COLD = 11.43  # navigation at -10 °C, fresh, cold-soaked
PROTO_MC = {"S13": dict(mean=39.53, sd=7.05, q=(28.36, 39.38, 51.34), s_off=0.70, brown=0.07),
            "S14": dict(mean=21.22, sd=4.77, q=(14.13, 20.87, 29.94), s_off=0.91, brown=0.11)}

CONST = ("S1", "S2a", "S2b", "S3", "S4", "S5a", "S5b", "S5c", "S6", "S7", "S8a", "S8b", "S16", "CM")
S0_SWEEP = (1.0, 0.8, 0.5, 0.2, 0.1)  # S15
S0_FINE = tuple(np.round(np.arange(0.05, 1.0001, 0.05), 2))
TGRID = (25, 10, 0, -10, -20)
SOHS = (1.0, 0.9, 0.8)


def key(S, T, soh, s0):
    return f"{S}|Tamb={T:g}|soh={soh:.1f}|s0={s0:.2f}"


def clean(x):
    """JSON-safe copy: NaN/inf -> None, NumPy scalars -> Python."""
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, np.ndarray):
        return clean(x.tolist())
    if isinstance(x, (np.floating, float)):
        x = float(x)
        return x if math.isfinite(x) else None
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def git_meta():
    def sh(*a):
        try:
            return subprocess.run(a, cwd=ROOT, capture_output=True, text=True, timeout=10).stdout.strip()
        except Exception:  # noqa: BLE001
            return ""
    return dict(git_sha=sh("git", "rev-parse", "--short", "HEAD"),
                core_uncommitted=bool(sh("git", "status", "--porcelain", "src/phonebatt")))


def record(sc, r, label) -> dict:
    """§14 AN3 fields of one Result (plus demand power and end values)."""
    pdem = sc.schedule.mean(lambda u: float(demand_power(u, sc.p)))
    return dict(label=label, Tamb_C=sc.Tamb_C, soh=round(1.0 - float(sc.p.q_loss), 3), s0=sc.s0,
                T0_C=sc.T0_C if sc.T0_C is not None else sc.Tamb_C,
                tte_h=r.tte_h, event=r.event, s_off_pct=100 * r.s_off, soe_off_pct=100 * r.soe_off,
                t_disp0_h=r.t_disp0_h, P_mean_W=r.P_mean_W, P_demand_W=pdem, P_end_W=r.P_end_W,
                T_end_C=r.T_end_C, Ts_end_C=r.Ts_end_C, V_end=r.V_end, rho_perf=r.rho_perf)


def main():
    t_start = time.time()
    apply()
    out: dict = {}
    trajs: dict = {}

    def do(label, sid, T=25.0, soh=1.0, s0=1.0, keep=False, T0_C=None, **opts):
        sc = scenario(sid, Tamb_C=float(T), s0=float(s0), soh=float(soh), T0_C=T0_C, **opts)
        r = sc.run(record=keep)
        k = key(label, T, soh, s0)
        out[k] = record(sc, r, label)
        if keep:
            trajs[k] = (r.traj, sc)
        return out[k]

    # ---- constant scenarios at 25 °C fresh, s0 in {1, 0.5, 0.2} (§5.2) and the S15 sweep for S4/S5a/S7
    for sid in CONST:
        for s0 in sorted(set((1.0, 0.5, 0.2) + (S0_SWEEP if sid in ("S4", "S5a", "S7") else ())), reverse=True):
            do(sid, sid, s0=s0, keep=(s0 == 1.0 and sid in ("S1", "S5a", "S7")))
    # ---- S7 at 15 / 35 °C and the no-throttling variant (throttling and burst envelope off, as in PROTO)
    for T in (15, 25, 35):
        if T != 25:
            do("S7", "S7", T=T, keep=True)
        do("S7-nothrottle", "S7", T=T, throttling=False, keep=True)
    # ---- navigation in the cold, fresh and cold-soaked (§5.4 "Other")
    do("S6", "S6", T=-10)
    # ---- S9 cold hike: T0 = 20 °C; full T x SOH grid (§5.4 reports {25, -10, -20} x {1.0, 0.8})
    for T in TGRID:
        for soh in SOHS:
            do("S9", "S9", T=T, soh=soh, keep=(T, soh) in ((25, 1.0), (-20, 0.8)))
    for T in (-10, -20):  # cold-soaked start
        do("S9-T0=Tamb", "S9", T=T, soh=0.8, T0_C=float(T))
    # ---- S10a (S4 at T x SOH) and S10b (S12 at 25 / 0 °C x SOH)
    for T in TGRID:
        for soh in SOHS:
            do("S10a", "S4", T=T, soh=soh)
    for T in (25, 0):
        for soh in SOHS:
            do("S10b", "S12", T=T, soh=soh)
    # ---- S11 weak signal: S11a SO-L and S11c VC over RSRP; S11b SO over Wi-Fi RSSI
    for rs in (-85, -95, -105, -110, -120):
        do(f"S11a(RSRP={rs})", "S11", variant="a", RSRP=rs)
        do(f"S11c(RSRP={rs})", "S11", variant="c", RSRP=rs)
    for rssi in (-85, -90):
        do(f"S11b(RSSI={rssi})", "S11", variant="b", RSSI=rssi)
    # ---- S12 scheduled day from the S15 initial SOCs
    for s0 in S0_SWEEP:
        do("S12", "S12", s0=s0)

    # ---- S15 fine curves: 25 °C fresh S4, S5a, S7, S12; cold/aged S4 (-10 / -20 °C, SOH 0.8) and S9 (-20, 0.8)
    curves = {}
    for label, sid, T, soh in (("S4", "S4", 25, 1.0), ("S5a", "S5a", 25, 1.0), ("S7", "S7", 25, 1.0),
                               ("S12", "S12", 25, 1.0), ("S4", "S4", -10, 0.8), ("S4", "S4", -20, 0.8),
                               ("S9", "S9", -20, 0.8)):
        tt, so = [], []
        for s0 in S0_FINE:
            k = key(label, T, soh, s0)
            rec = out.get(k)
            if rec is None:
                sc = scenario(sid, Tamb_C=float(T), s0=float(s0), soh=float(soh))
                rec = record(sc, sc.run(), label)
                if s0 in S0_SWEEP or s0 in (0.5, 0.2):
                    out[k] = rec
            tt.append(rec["tte_h"])
            so.append(rec["s_off_pct"])
        curves[f"{label}|Tamb={T:g}|soh={soh:.1f}"] = dict(s0=list(S0_FINE), tte_h=tt, s_off_pct=so)

    # ---- S13 / S14 stochastic days (CTMC, L2 batch, 400 paths, h = 10 s, seed 7); S13 also over S15's s0
    mc = {}
    for sid, prof in (("S13", "typical"), ("S14", "heavy")):
        ctmc = CTMC(prof)
        paths = None
        for s0 in (S0_SWEEP if sid == "S13" else (1.0,)):
            b = run_batch(ctmc, 25.0, None, s0=s0, h=H_MC, n=N_PATHS, seed=SEED, paths=paths)
            paths = b["paths"]  # common random numbers across s0
            tt = b["tte_h"]
            ev = b["event"]
            q05, q50, q95 = np.percentile(tt, [5, 50, 95])
            rec = dict(label=sid, Tamb_C=25.0, soh=1.0, s0=s0, T0_C=25.0, tte_h=float(np.mean(tt)),
                       event="mixed", event_share={e: float(np.mean(ev == e)) for e in sorted(set(ev))},
                       s_off_pct=float(100 * np.mean(b["s_off"])), soe_off_pct=float(100 * np.mean(b["soe_off"])),
                       t_disp0_h=None, P_mean_W=None, T_end_C=float(np.mean(b["T_end_C"])),
                       Ts_end_C=float(np.nanmean(b["Ts_end_C"])), rho_perf=None,
                       median_h=float(q50), sd_h=float(np.std(tt)), q05_h=float(q05), q95_h=float(q95),
                       p_lt_24h=float(np.mean(tt < 24)), n_paths=N_PATHS, h_s=H_MC, seed=SEED,
                       note="tte_h is the Monte Carlo mean over CTMC paths (rk4_batch); P_mean not tracked")
            out[key(sid, 25, 1.0, s0)] = rec
            mc[(sid, s0)] = tt
    s13 = [out[key("S13", 25, 1.0, s)] for s in sorted(S0_SWEEP)]
    curves["S13|Tamb=25|soh=1.0"] = dict(s0=sorted(S0_SWEEP), tte_h=[r["tte_h"] for r in s13],
                                         q05_h=[r["q05_h"] for r in s13], q95_h=[r["q95_h"] for r in s13],
                                         s_off_pct=[r["s_off_pct"] for r in s13])
    s15 = {c: dict(v, **curvature(v)) for c, v in curves.items()}

    comp = proto_comparison(out)
    runtime = time.time() - t_start
    res = dict(meta=dict(spec=SPEC, seed=SEED, **git_meta(),
                         solver="LSODA (rtol 1e-7, atol [1e-9,1e-7,1e-5], max_step 60 s) for deterministic runs; "
                                "rk4_batch for the CTMC days S13/S14",
                         h_s=H_MC, h_s_note="h applies to the S13/S14 Monte Carlo only (LSODA is adaptive)",
                         n_paths=N_PATHS, runtime_s=round(runtime, 1),
                         key_format="{S}|Tamb={degC}|soh={SOH:.1f}|s0={s0:.2f}; S10a = S4 and S10b = S12 on an "
                                    "aged pack; S9 starts at T0 = 20 °C unless 'S9-T0=Tamb'; S7-nothrottle has "
                                    "throttling and the burst envelope disabled (as in final_proto.py)",
                         units="s_off_pct = displayed SOC D1 at shutdown; soe_off_pct = energy-fraction SOC D2; "
                               "P_mean_W = delivered energy / TTE; P_demand_W = schedule-mean demand before "
                               "throttling; rho_perf = time-averaged delivered SoC performance"),
               **{k: out[k] for k in sorted(out)}, S15=s15, proto_comparison=comp)
    RES.mkdir(exist_ok=True)
    (RES / "scenarios.json").write_text(json.dumps(clean(res), indent=1))

    fig_trajectories(trajs)
    fig_throttle(trajs)
    fig_tte_s0(s15)
    table_baseline(out)
    n_bad = sum(1 for c in comp if not c["within_tol"])
    print(f"scenarios: {len(out)} runs, {len(comp)} PROTO comparisons, {n_bad} outside tolerance, "
          f"{time.time() - t_start:.1f} s")
    for c in comp:
        if not c["within_tol"]:
            print("  OUTSIDE:", c)


# ------------------------------------------------------------------------------------------- S15 curvature
def curvature(c: dict) -> dict:
    """Shape of TTE(s0) relative to proportionality TTE(s0) = s0·TTE(100 %) (§10.5 item 4).

    * ratio_50, ratio_20: TTE(s0)/TTE(100 %) at s0 = 50 / 20 %;
    * prop_dev_*: largest deviation TTE(s0)/TTE(100 %) − s0 (percentage points of TTE(100 %)) and where;
    * chord_50_100_intercept_pct: s0 where the chord through s0 = 50 % and 100 % reaches TTE = 0 (an effective
      empty display level; ≈ s_off when TTE is linear with an offset);
    * second_diff: [TTE(100) − 2·TTE(60) + TTE(20)]/TTE(100) (0 for a straight line, > 0 convex);
    * s0_min_pct: smallest grid s0 that lasts more than 3 min (shutdown at the start below it).
    """
    s = np.asarray(c["s0"], float)
    t = np.asarray(c["tte_h"], float)
    o = np.argsort(s)
    s, t = s[o], t[o]
    t1 = float(np.interp(1.0, s, t))
    f = lambda x: float(np.interp(x, s, t))  # noqa: E731
    dev = t / t1 - s
    j = int(np.argmax(np.abs(dev)))
    slope = (t1 - f(0.5)) / 0.5
    ok = t > 0.05
    out = dict(tte_100_h=t1, ratio_50=f(0.5) / t1, ratio_20=f(0.2) / t1,
               prop_dev_max_pp=float(100 * dev[j]), prop_dev_at_s0_pct=float(100 * s[j]),
               chord_50_100_intercept_pct=float(100 * (0.5 - f(0.5) / slope)) if slope > 0 else None,
               s0_min_pct=float(100 * s[ok].min()) if ok.any() else None)
    if s.min() <= 0.2 + 1e-9 and np.any(np.isclose(s, 0.6)):
        out["second_diff"] = (t1 - 2 * f(0.6) + f(0.2)) / t1
    return out


# ------------------------------------------------------------------------------------------- PROTO table
def proto_comparison(out: dict) -> list:
    rows = []

    def add(k, qty, model, proto, tol, kind="rel"):
        if proto is None or model is None:
            return
        if isinstance(proto, str):
            ok = model == proto
            rows.append(dict(key=k, qty=qty, model=model, proto=proto, diff=None, within_tol=ok))
            return
        if kind == "rel":
            d = 100 * (model / proto - 1)
            ok = abs(d) <= 100 * tol
            rows.append(dict(key=k, qty=qty, model=model, proto=proto, diff_pct=d, tol_pct=100 * tol, within_tol=ok))
        else:
            d = model - proto
            ok = abs(d) <= tol
            rows.append(dict(key=k, qty=qty, model=model, proto=proto, diff_abs=d, tol_abs=tol, within_tol=ok))

    def main_row(k, tgt):
        r = out[k]
        tte, ev, soff, Te, Tse, P = tgt
        add(k, "tte_h", r["tte_h"], tte, 0.01)
        add(k, "event", r["event"], ev, None)
        add(k, "s_off_pct", r["s_off_pct"], soff, 0.3, "abs")
        add(k, "T_end_C", r["T_end_C"], Te, 0.06, "abs")  # PROTO printed to 0.1 K
        add(k, "Ts_end_C", r["Ts_end_C"], Tse, 0.06, "abs")
        return r

    for sid, tgt in PROTO_MAIN.items():
        r = main_row(key(sid, 25, 1.0, 1.0), tgt)
        add(key(sid, 25, 1.0, 1.0), "P_demand_W", r["P_demand_W"], tgt[5], 0.0006, "abs")
    for T, tgt in PROTO_S7.items():
        r = main_row(key("S7", T, 1.0, 1.0), tgt)
        add(key("S7", T, 1.0, 1.0), "P_end_W", r["P_end_W"], tgt[5], 0.006, "abs")
    add(key("S7-nothrottle", 25, 1.0, 1.0), "tte_h", out[key("S7-nothrottle", 25, 1.0, 1.0)]["tte_h"],
        PROTO_S7_NOTHROTTLE_25, 0.01)
    for s0, tte in PROTO_S12.items():
        add(key("S12", 25, 1.0, s0), "tte_h", out[key("S12", 25, 1.0, s0)]["tte_h"], tte, 0.01)
        add(key("S12", 25, 1.0, s0), "event", out[key("S12", 25, 1.0, s0)]["event"], "cutoff", None)
    add(key("S12", 25, 1.0, 1.0), "s_off_pct", out[key("S12", 25, 1.0, 1.0)]["s_off_pct"], PROTO_S12_SOFF, 0.3,
        "abs")
    add(key("S12", 25, 1.0, 1.0), "P_demand_W", out[key("S12", 25, 1.0, 1.0)]["P_demand_W"], PROTO_S12_PDEM,
        0.0006, "abs")
    for (T, soh), (tte, soff) in PROTO_S9.items():
        k = key("S9", T, soh, 1.0)
        tol, tol_s = (0.02, 1.0) if (T, soh) == (-20, 0.8) else (0.01, 0.3)
        add(k, "tte_h", out[k]["tte_h"], tte, tol)
        add(k, "s_off_pct", out[k]["s_off_pct"], soff, tol_s, "abs")
        add(k, "event", out[k]["event"], "cutoff", None)
    add(key("S9-T0=Tamb", -10, 0.8, 1.0), "tte_h", out[key("S9-T0=Tamb", -10, 0.8, 1.0)]["tte_h"],
        PROTO_S9_COLDSOAK, 0.01)
    for (T, soh), (soff, tte) in PROTO_S10A.items():
        k = key("S10a", T, soh, 1.0)
        add(k, "tte_h", out[k]["tte_h"], tte, 0.01)
        add(k, "s_off_pct", out[k]["s_off_pct"], soff, 0.3, "abs")
        add(k, "event", out[k]["event"], "cutoff" if (T, soh) == (25, 1.0) else "brownout", None)
    for (T, soh), (tte, ev, soff) in PROTO_S10B.items():
        k = key("S10b", T, soh, 1.0)
        add(k, "tte_h", out[k]["tte_h"], tte, 0.01)
        add(k, "event", out[k]["event"], ev, None)
        add(k, "s_off_pct", out[k]["s_off_pct"], soff, 0.3, "abs")
    for rs, tte in PROTO_S11A.items():
        add(key(f"S11a(RSRP={rs})", 25, 1.0, 1.0), "tte_h", out[key(f"S11a(RSRP={rs})", 25, 1.0, 1.0)]["tte_h"],
            tte, 0.01)
    for rs, tte in PROTO_S11C.items():
        add(key(f"S11c(RSRP={rs})", 25, 1.0, 1.0), "tte_h", out[key(f"S11c(RSRP={rs})", 25, 1.0, 1.0)]["tte_h"],
            tte, 0.01)
    add(key("S6", -10, 1.0, 1.0), "tte_h", out[key("S6", -10, 1.0, 1.0)]["tte_h"], PROTO_S6_COLD, 0.01)
    for sid, t in PROTO_MC.items():  # 3 % tolerance (§15.3, seed-dependent)
        k = key(sid, 25, 1.0, 1.0)
        r = out[k]
        add(k, "mc_mean_h", r["tte_h"], t["mean"], 0.03)
        add(k, "mc_median_h", r["median_h"], t["q"][1], 0.03)
        add(k, "mc_sd_h", r["sd_h"], t["sd"], 0.10)
        add(k, "mc_brownout_share", r["event_share"].get("brownout", 0.0), t["brown"], 0.05, "abs")
    return rows


# ------------------------------------------------------------------------------------------- figures
def _cycles(traj, period=600.0):
    """Per-cycle min / max of V for a periodic schedule (S9 burst band)."""
    t = traj["t"]
    cyc = np.floor(t / period).astype(int)
    ids = np.unique(cyc)
    tc = (ids + 0.5) * period / 3600.0
    vmin = np.array([traj["V"][cyc == i].min() for i in ids])
    vmax = np.array([traj["V"][cyc == i].max() for i in ids])
    return tc, vmin, vmax


def fig_trajectories(trajs):
    cols = [("S1", key("S1", 25, 1.0, 1.0), "S1 standby"),
            ("S5a", key("S5a", 25, 1.0, 1.0), "S5a video, Wi-Fi"),
            ("S7", key("S7", 25, 1.0, 1.0), "S7 3D gaming"),
            ("S9", key("S9", -20, 0.8, 1.0), "S9 cold hike, −20 °C")]
    ref9 = key("S9", 25, 1.0, 1.0)
    fig, axs = plt.subplots(3, 4, figsize=(W1, 3.75), sharex="col", gridspec_kw=dict(hspace=0.14, wspace=0.30))
    for j, (sid, k, title) in enumerate(cols):
        tr, sc = trajs[k]
        t = tr["t_h"]
        a0, a1, a2 = axs[:, j]
        a0.set_title(title, fontsize=8.5, loc="left", color=INK, pad=3)
        panel_label(a0, f"({'abcd'[j]})")
        if sid == "S9":
            tr2 = trajs[ref9][0]
            a0.plot(tr2["t_h"], 100 * tr2["s"], color=MUTED, lw=1.2)
            a0.plot(t, 100 * tr["s"], color=BLUE, lw=1.6)
            a0.text(14.5, 52, "25 °C, SOH 1.0", fontsize=7.0, color=INK2, ha="left", va="bottom")
            a0.text(t[-1] - 0.8, 100 * tr["s"][-1] - 5, "−20 °C,\nSOH 0.8", fontsize=7.0, color=INK2,
                    ha="right", va="top", linespacing=1.0)
            xmax = tr2["t_h"][-1] * 1.03
        else:
            a0.plot(t, 100 * tr["s"], color=BLUE, lw=1.6)
            xmax = t[-1] * 1.03
        s_off = 100 * tr["s"][-1]
        a0.plot([t[-1]], [s_off], "o", color=INK, ms=3)
        # stranded-charge label in the empty lower-left corner (the SOC curve runs top-left to bottom-right)
        a0.text(0.04, 0.05, (f"$s_{{off}}$ = {s_off:.1f} %" if s_off >= 1 else f"$s_{{off}}$ = {s_off:.2f} %"),
                transform=a0.transAxes, ha="left", va="bottom", fontsize=7.4, color=INK)
        a0.set_ylim(-4, 106)
        # terminal voltage (S9: per-cycle band between the photo-burst minimum and the tracking maximum)
        if sid == "S9":
            tc, vmin, vmax = _cycles(tr)
            a1.fill_between(tc, vmin, vmax, color=BLUE, alpha=0.18, lw=0)
            a1.plot(tc, vmax, color=BLUE, lw=0.9)
            a1.plot(tc, vmin, color=BLUE, lw=0.9)
            i, i2 = len(tc) // 2, len(tc) // 8
            a1.annotate("tracking", (tc[i], vmax[i]), xytext=(4, 3), textcoords="offset points", fontsize=7.0,
                        color=INK2)
            a1.annotate("30 s photo bursts", (tc[i2], vmin[i2]), xytext=(0.6, 3.56), textcoords="data",
                        fontsize=7.0, color=INK2, va="center",
                        arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6, shrinkA=1, shrinkB=1))
        else:
            a1.plot(t, tr["V"], color=BLUE, lw=1.4)
        a1.axhline(float(sc.p.V_cut), color=MUTED, lw=0.8, ls="--")
        a1.set_ylim(3.2, 4.5)
        if j == 0:
            a1.text(xmax * 0.03, float(sc.p.V_cut) - 0.03, "$V_{cut}$", fontsize=7.0, color=INK2, va="top")
        # temperatures (same scale for the three 25 °C columns)
        a2.plot(t, tr["T_C"], color=BLUE, lw=1.4, label="battery $T$")
        a2.plot(t, tr["Ts_C"], color=ORANGE, lw=1.4, label="skin $T_s$")
        if sid == "S7":
            a2.axhline(45.0, color=MUTED, lw=0.8, ls="--")
            a2.text(t[-1], 45.4, "$T_{lim}$", fontsize=7.0, color=INK2, ha="right", va="bottom")
        a2.set_ylim((-22, 24) if sid == "S9" else (20, 48))
        a2.set_xlabel("time (h)")
        a2.set_xlim(0, xmax)
        for a in (a0, a1, a2):
            a.tick_params(labelsize=7.2)
    axs[0, 0].set_ylabel("SOC $s$ (%)")
    axs[1, 0].set_ylabel("$V$ (V)")
    axs[2, 0].set_ylabel("$T$, $T_s$ (°C)")
    axs[2, 0].legend(loc="upper right", fontsize=7.0, handlelength=1.3)
    fig.align_ylabels(axs[:, 0])
    save(fig, "fig_trajectories")


def fig_throttle(trajs):
    """S7 at 15 / 25 / 35 °C (diverging blue-gray-red around the 25 °C reference); log time axis so that the
    7.5 min thermal transient and the multi-hour discharge are both visible."""
    tcol = {15: SEQ_BLUE[9], 25: INK2, 35: RED_ARM[3]}
    fig, axs = plt.subplots(1, 3, figsize=(W1, 2.3), gridspec_kw=dict(wspace=0.33))
    ticks = ([1 / 60, 10 / 60, 1.0, 10.0], ["1 min", "10 min", "1 h", "10 h"])
    for T in (15, 25, 35):
        tr, sc = trajs[key("S7", T, 1.0, 1.0)]
        m = tr["t_h"] > 0
        t = tr["t_h"][m]
        rp = rho_perf(tr["Psoc"][m], sc.schedule[0][1].x_c, sc.p)
        lab = f"{T} °C"
        axs[0].plot(t, tr["Ts_C"][m], color=tcol[T], lw=1.6, label=lab)
        axs[1].plot(t, tr["P"][m], color=tcol[T], lw=1.6, label=lab)
        axs[2].plot(t, rp, color=tcol[T], lw=1.6, label=lab)
        for a, y in ((axs[0], tr["Ts_C"][-1]), (axs[1], tr["P"][-1]), (axs[2], rp[-1])):
            a.plot([t[-1]], [y], "o", color=tcol[T], ms=3)
        if T == 15:
            axs[1].annotate("brown-out", (t[-1], tr["P"][-1]), xytext=(3, 3), textcoords="offset points",
                            fontsize=6.3, color=INK2)
    tr, sc = trajs[key("S7-nothrottle", 25, 1.0, 1.0)]
    m = tr["t_h"] > 0
    for a, y in ((axs[0], tr["Ts_C"]), (axs[1], tr["P"])):
        a.plot(tr["t_h"][m], y[m], color=INK2, lw=1.1, ls=(0, (3, 2)), label="25 °C, no throttling")
    axs[0].axhline(45.0, color=MUTED, lw=0.8, ls="--")
    axs[0].text(0.011, 45.6, "$T_{lim}$ = 45 °C", fontsize=6.3, color=INK2, ha="left", va="bottom")
    axs[0].set_ylabel("skin temperature $T_s$ (°C)")
    axs[1].set_ylabel("battery power $P$ (W)")
    axs[2].set_ylabel(r"SoC performance $\rho_{perf}$")
    axs[0].set_ylim(10, 52)
    axs[1].set_ylim(0, 4.6)
    axs[2].set_ylim(0, 1.08)
    for a, lab in zip(axs, "abc"):
        a.set_xscale("log")
        a.set_xlim(0.01, 13)
        a.set_xticks(*ticks)
        a.minorticks_off()
        a.set_xlabel("time since start")
        panel_label(a, f"({lab})")
    axs[1].legend(loc="lower left", fontsize=6.3, ncol=2, handlelength=1.6, columnspacing=0.8)
    save(fig, "fig_throttle")


def fig_tte_s0(s15):
    fig, axs = plt.subplots(1, 2, figsize=(W1, 2.45), gridspec_kw=dict(wspace=0.28))
    a = axs[0]
    names = [("S4", "S4 social"), ("S5a", "S5a video"), ("S7", "S7 gaming"), ("S12", "S12 day"), ("S13", "S13 typical day (MC)")]
    for i, (sid, lab) in enumerate(names):
        c = s15[f"{sid}|Tamb=25|soh=1.0"]
        s = 100 * np.asarray(c["s0"])
        if sid == "S13":
            y = np.asarray(c["tte_h"])
            a.errorbar(s, y, yerr=[y - np.asarray(c["q05_h"]), np.asarray(c["q95_h"]) - y], color=CAT[i], lw=1.4,
                       marker="o", ms=3.2, capsize=1.8, elinewidth=0.8, label=lab)
        else:
            a.plot(s, c["tte_h"], color=CAT[i], lw=1.6, label=lab)
    a.set_xlabel("initial displayed SOC $s_0$ (%)")
    a.set_ylabel("TTE (h)")
    a.set_xlim(0, 102)
    a.set_ylim(0, 56)
    a.legend(loc="upper left", fontsize=6.5, handlelength=1.5)
    a.set_title("25 °C, new battery", fontsize=8, loc="left", color=INK)
    panel_label(a, "(a)")
    # (b) TTE normalised by TTE(100 %): proportional line vs cold / aged curves
    b = axs[1]
    b.plot([0, 100], [0, 1], color=GRID, lw=3.0, zorder=1, label="proportional")
    sel = [("S4|Tamb=25|soh=1.0", "S4, 25 °C, SOH 1.0", CAT[0], "-"),
           ("S4|Tamb=-10|soh=0.8", "S4, −10 °C, SOH 0.8", CAT[0], (0, (4, 1.5))),
           ("S4|Tamb=-20|soh=0.8", "S4, −20 °C, SOH 0.8", CAT[0], (0, (1.2, 1.2))),
           ("S12|Tamb=25|soh=1.0", "S12, 25 °C, SOH 1.0", CAT[3], "-"),
           ("S9|Tamb=-20|soh=0.8", "S9, −20 °C, SOH 0.8", CAT[6], "-")]
    for kk, lab, col, ls in sel:
        c = s15[kk]
        s = 100 * np.asarray(c["s0"])
        y = np.asarray(c["tte_h"]) / c["tte_100_h"]
        b.plot(s, y, color=col, lw=1.5, ls=ls, label=lab)
    b.set_xlabel("initial displayed SOC $s_0$ (%)")
    b.set_ylabel("TTE($s_0$) / TTE(100 %)")
    b.set_xlim(0, 102)
    b.set_ylim(0, 1.03)
    b.legend(loc="upper left", fontsize=6.5, handlelength=2.0)
    panel_label(b, "(b)")
    save(fig, "fig_tte_s0")


# ------------------------------------------------------------------------------------------- LaTeX table
EV = {"cutoff": "cut-off", "brownout": "brown-out", "thermal": "thermal", "infeasible": "infeasible",
      "charge": "charge", "none": "none"}


def table_baseline(out):
    rows = [("S1", "Standby, home Wi-Fi"), ("S2a", "Standby, cellular"), ("S2b", "Standby, 5G, heavy push"),
            ("S3", "Reading"), ("S4", "Social feed, Wi-Fi"), ("S5a", "Video, Wi-Fi"), ("S5b", "Video, LTE"),
            ("S5c", "Video, 5G NSA"), ("S6", "Navigation"), ("S7", "3D gaming"), ("S8a", "Video call, LTE"),
            ("S8b", "Video call, Wi-Fi"), ("CM", "Camera"), ("S12", "Scheduled day"), ("S16", "Audio, screen off")]
    L = [r"\begin{tabular}{llrrlrr}", r"\toprule",
         r"ID & Scenario & $\bar P$ (W) & TTE (h) & Event & $s_\mathrm{off}$ (\%) & $T_\mathrm{end}$ / $T_{s,\mathrm{end}}$ ($^\circ$C) \\",
         r"\midrule"]
    for sid, name in rows:
        r = out[key(sid, 25, 1.0, 1.0)]
        P = r["P_demand_W"]
        Ps = f"{P:.3f}" if sid != "S7" else f"{P:.2f}$\\to${r['P_mean_W']:.2f}"
        tte = r["tte_h"]
        ts = f"{tte:.1f}" if tte >= 100 else f"{tte:.2f}"
        L.append(f"{sid} & {name} & {Ps} & {ts} & {EV[r['event']]} & {r['s_off_pct']:.2f} & "
                 f"{r['T_end_C']:.1f} / {r['Ts_end_C']:.1f} \\\\")
    for T in (15, 35):
        r = out[key("S7", T, 1.0, 1.0)]
        L.append(f"S7 & 3D gaming, {T}\\,$^\\circ$C & {r['P_demand_W']:.2f}$\\to${r['P_mean_W']:.2f} & "
                 f"{r['tte_h']:.2f} & {EV[r['event']]} & {r['s_off_pct']:.2f} & {r['T_end_C']:.1f} / "
                 f"{r['Ts_end_C']:.1f} \\\\")
    r = out[key("S7-nothrottle", 25, 1.0, 1.0)]
    L.append(f"S7 & 3D gaming, no throttling & {r['P_demand_W']:.2f} & {r['tte_h']:.2f} & {EV[r['event']]} & "
             f"{r['s_off_pct']:.2f} & {r['T_end_C']:.1f} / {r['Ts_end_C']:.1f} \\\\")
    L += [r"\bottomrule", r"\end{tabular}"]
    TAB.mkdir(parents=True, exist_ok=True)
    (TAB / "tab_baseline.tex").write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
