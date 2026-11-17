"""AN14: recommendations - user levers, scenario-specific levers and OS strategies O1-O9 (FINAL_SPEC §12, §14).

Writes
* results/recommendations.json
  - ``base``: S12 P̄ (quasi-steady state powers), E_u/P̄, S13 Dynkin mean
  - ``levers.{id}.{label, kind, s12_pct, ctmc_pct, s12_ci90_pct, ctmc_ci90_pct, ctmc_ci90_usage_pct, p_gain_s12,
    p_gain_ctmc, rank_ctmc_median, rank_ctmc_ci90, s12_L2_pct, proto, dev_pp}``.  Central values are the smooth metrics of
    §12.1 (``analytic.lever_effects``); intervals are the 5-95 % range over a joint Sobol' sample of the §4 parameter
    distributions (N = 1024, seed 2026; ``ctmc_ci90_usage_pct`` adds the §4.6 usage-model factors).  ``s12_L2_pct`` is
    the non-smooth L2 TTE of the S12 day (LSODA), reported to show why it is not used for ranking.
  - ``ranking``: levers sorted by the S13 metric, with rank stability over the parameter sample
  - ``scenario_levers``: L2 LSODA comparisons (video/video call over Wi-Fi vs cellular, keeping warm on the cold
    hike S9, gaming in a hot car, weak signal, standby RAT), with parameter intervals of the quasi-steady power ratio
    where the comparison is a power ratio
  - ``os.{O1..O9}``: each OS strategy quantified as in §12.2
  - ``energy_share_S12``: component energy shares of the S12 day (quasi-steady state powers)
* paper/figures/fig_levers.pdf (+ png): ranked horizontal bars with 90 % parameter intervals, (a) S12 E_u/P̄,
  (b) S13 Dynkin mean; gains in blue (OS strategies in a lighter blue), 'avoid' items in red

Run:  cd /home/user/MCM && PYTHONPATH=src:analysis python3 analysis/run_recommendations.py   (~1 min)
      (--figure-only redraws fig_levers from results/recommendations.json)
"""
from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import replace

import matplotlib.pyplot as plt
import numpy as np

from basics_common import RESULTS, meta, write_json
from phonebatt import CTMC, LEVERS, RP63, aged, run, run_batch
from phonebatt import aging as G
from phonebatt import analytic as A
from phonebatt import cell
from phonebatt.params import UQ_SPEC, USAGE_UQ_SPEC, ACTIVITIES, sample, with_theta
from phonebatt.power import BREAKDOWN_KEYS, demand_power, power_parts, steady_power, steady_skin_temperature, \
    throttle
from phonebatt.usage import PRIMITIVES, PROFILES, Schedule, s9_schedule, s12_schedule
from style import BLUE, INK, INK2, MUTED, RED, SEQ_BLUE, W1, apply, panel_label, save

SEED = 2026
N_THETA = 1024
TAMB = 25.0

# §12.1 [PROTO] lever table: (S12 E_u/P̄ %, S13 Dynkin %); None = not given
PROTO_LEVERS = {"skip_game": (34.9, 36.9), "dark_theme": (7.6, 4.6), "voice_call": (4.8, 4.0),
                "brightness_m01": (4.4, 3.1), "refresh_60": (4.1, 4.7), "vcall_wifi": (3.9, 1.6),
                "half_bg": (0.1, 0.4), "standby_5g": (-4.4, -4.5), "weak_signal": (-3.9, -5.7),
                "brightness_p02": (-10.3, -7.5), "os_coalesce": (0.1, 0.7), "os_frame_cap": (3.4, 2.7),
                "wifi_share_80": (None, 3.7)}
# §12.1 [PROTO] L2 TTE of the S12 day with a lever (final_proto.py 'recs'; lever maps of that section)
PROTO_S12_L2 = {"base": 34.92, "skip_game": 37.7, "dark_theme": 1.4, "brightness_m01": 0.8, "refresh_60": 0.7,
                "voice_call": 0.9, "half_bg": 0.0, "standby_5g": -1.5, "weak_signal": -0.8}
KIND = {k: ("avoid" if k in ("standby_5g", "weak_signal", "brightness_p02") else
            "os" if k.startswith("os_") else "gain") for k in list(LEVERS) + ["wifi_share_80"]}
SHORT = {"skip_game": "Skip 3D gaming (1 h)", "dark_theme": "Dark theme", "voice_call": "Voice call, not video",
         "brightness_m01": "Brightness −0.1", "refresh_60": "Refresh capped at 60 Hz",
         "vcall_wifi": "Video call on Wi-Fi, not LTE", "half_bg": "Halve background wake-ups",
         "standby_5g": "Avoid: 5G NSA standby", "weak_signal": "Avoid: weak signal (−110 dBm)",
         "brightness_p02": "Avoid: brightness +0.2", "os_coalesce": "OS: coalesce background ×5",
         "os_frame_cap": "OS: game frame cap", "wifi_share_80": "Prefer Wi-Fi (50 → 80 % of time)"}


# ------------------------------------------------------------------------------------------- vectorised metrics
def deviation(L):
    """Stationary distribution and deviation matrix of one generator (16 x 16) or a batch (N, 16, 16)."""
    L = np.asarray(L, float)
    if L.ndim == 2:
        pi = A.stationary(L)
        Pi = np.outer(np.ones(pi.size), pi)
        return pi, np.linalg.inv(Pi - L) - Pi
    n = L.shape[-1]
    M = np.concatenate([np.swapaxes(L, 1, 2), np.ones((L.shape[0], 1, n))], axis=1)
    b = np.zeros((n + 1,))
    b[-1] = 1.0
    pi = np.linalg.pinv(M) @ b  # (N, n), least squares as in analytic.stationary
    Pi = pi[:, None, :] * np.ones((1, n, 1))
    return pi, np.linalg.inv(Pi - L) - Pi


def dynkin(pi, Dm, P, Eu, x0=0):
    """Dynkin mean (F6) vectorised over samples.  P (16, N); pi (16,) or (N, 16); Dm (16,16) or (N,16,16)."""
    if np.ndim(pi) == 1:
        Pbar = pi @ P
        hv = Dm @ (P - Pbar)
        hP = np.sum(pi[:, None] * P * hv, axis=0) / Pbar
    else:
        Pt = P.T  # (N, 16)
        Pbar = np.sum(pi * Pt, axis=1)
        hv = np.einsum("nij,nj->ni", Dm, Pt - Pbar[:, None])
        hP = np.sum(pi * Pt * hv, axis=1) / Pbar
        hv = hv.T
    return (Eu - hv[x0] + hP) / Pbar


class Engine:
    """Smooth lever metrics for array-valued parameters (one value per sample)."""

    def __init__(self, pb, Tamb_C=TAMB):
        self.pb, self.Ta = pb, Tamb_C
        self.sch = s12_schedule()
        self.ct = CTMC("typical")
        self.us = self.ct.primitives()
        self.Eu = cell.usable_energy(pb)
        self.pi, self.Dm = deviation(self.ct.generator())
        self.pi80, self.Dm80 = deviation(CTMC("typical", mean_w=4.0, mean_c=1.0).generator())
        self.cache = {}

    def P(self, lever, name, u):
        key = (lever, name)
        if key not in self.cache:
            uu = u if lever is None else LEVERS[lever][1](name, u)
            self.cache[key] = np.atleast_1d(steady_power(uu, self.pb, self.Ta))
        return self.cache[key]

    def s12_pbar(self, lever=None):
        return sum(d * self.P(lever, u.name, u) for d, u in self.sch) / self.sch.period

    def state_P(self, lever=None):
        return np.stack([self.P(lever, n, u) for n, u in zip(self.ct.state_map, self.us)])

    def ctmc_mu(self, lever=None, pi=None, Dm=None):
        pi = self.pi if pi is None else pi
        Dm = self.Dm if Dm is None else Dm
        return dynkin(pi, Dm, self.state_P(lever), self.Eu)

    def all(self, gens=None, gens80=None):
        """% gains of every lever: (s12, ctmc) arrays over samples.  ``gens``: per-sample generators (usage UQ)."""
        pi, Dm = (self.pi, self.Dm) if gens is None else deviation(gens)
        pi80, Dm80 = (self.pi80, self.Dm80) if gens80 is None else deviation(gens80)
        P0, mu0 = self.s12_pbar(), self.ctmc_mu(None, pi, Dm)
        out = {}
        for k in LEVERS:
            out[k] = (100 * (P0 / self.s12_pbar(k) - 1), 100 * (self.ctmc_mu(k, pi, Dm) / mu0 - 1))
        mu80 = dynkin(pi80, Dm80, self.state_P(None), self.Eu)
        out["wifi_share_80"] = (np.full_like(mu80, np.nan), 100 * (mu80 / mu0 - 1))
        return out, P0, mu0


def usage_generators(th, n, mean_w=2.5, mean_c=2.5):
    """Per-sample CTMC generators from the §4.6 usage factors (typical profile)."""
    lam0, _, d0 = PROFILES["typical"]
    gens = np.empty((n, 16, 16))
    for i in range(n):
        p = {a: float(th["p_typical"][j, i]) for j, a in enumerate(ACTIVITIES)}
        d = {a: d0[a] * float(th[f"d_mul_{a}"][i]) for a in ACTIVITIES}
        c = CTMC("typical", mean_w=mean_w * float(th["conn_mul_wifi"][i]),
                 mean_c=mean_c * float(th["conn_mul_cell"][i]), lam_on=lam0 * float(th["lam_on_mul"][i]), p=p,
                 d_min=d)
        gens[i] = c.generator()
    return gens


def ci(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return [float(v) for v in np.percentile(x, (5, 95))] if x.size else [None, None]


def pct(a, b):
    return 100.0 * (a / b - 1.0)


# ------------------------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figure-only", action="store_true")
    a = ap.parse_args(argv)
    if a.figure_only:
        make_figure(json.loads((RESULTS / "recommendations.json").read_text()))
        return
    t_start = time.time()
    out: dict = {}

    # ---- central smooth metrics (core) and a vectorised cross-check
    LE = A.lever_effects(RP63, TAMB)
    eng0 = Engine(RP63)
    cen, P0, mu0 = eng0.all()
    xcheck = max(max(abs(float(cen[k][0][0]) - LE[k]["s12_pct"]) if k != "wifi_share_80" else 0.0,
                     abs(float(cen[k][1][0]) - LE[k]["ctmc_pct"])) for k in cen)
    out["base"] = dict(LE["base"], Eu_Wh=float(cell.usable_energy(RP63)), vectorised_max_abs_diff_pp=xcheck)
    print(f"base: S12 Pbar {LE['base']['s12_Pbar_W']:.4f} W, E_u/Pbar {LE['base']['s12_tte_inf_h']:.2f} h, "
          f"S13 Dynkin {LE['base']['ctmc_mu_h']:.2f} h; vectorised vs core max diff {xcheck:.2e} pp")

    # ---- parameter sample (θ: §4) and usage factors (§4.6), one joint scrambled Sobol' design
    spec = dict(UQ_SPEC)
    spec.update({k: v for k, v in USAGE_UQ_SPEC.items() if k != "p_heavy"})
    th_all = sample(spec, N_THETA, method="sobol", seed=SEED)
    theta = {k: th_all[k] for k in UQ_SPEC}
    pb = with_theta(RP63, theta)
    eng = Engine(pb)
    res_t, P0t, mu0t = eng.all()
    gens = usage_generators(th_all, N_THETA)
    gens80 = usage_generators(th_all, N_THETA, mean_w=4.0, mean_c=1.0)
    res_u, _, mu0u = eng.all(gens, gens80)
    print(f"UQ: {N_THETA} samples; S13 Dynkin base q05-q95 {ci(mu0t)} h (θ), {ci(mu0u)} h (θ + usage)")

    # ranks (1 = largest gain) per sample on the S13 metric
    ids = list(res_t)
    M = np.stack([res_t[k][1] for k in ids])  # (K, N)
    ranks = (-M).argsort(axis=0).argsort(axis=0) + 1
    Ms = np.stack([res_t[k][0] for k in ids if k != "wifi_share_80"])
    ranks_s12 = (-Ms).argsort(axis=0).argsort(axis=0) + 1
    ids_s12 = [k for k in ids if k != "wifi_share_80"]

    # ---- non-smooth L2 TTE of the S12 day with each lever (LSODA)
    sch = s12_schedule()
    base_l2 = run(sch, TAMB).tte_h
    l2 = {}
    for k, (label, fn) in LEVERS.items():
        l2[k] = run(sch.map(lambda u, fn=fn: fn(u.name, u)), TAMB).tte_h

    levers = {}
    for j, k in enumerate(ids):
        s12c = float("nan") if k == "wifi_share_80" else LE[k]["s12_pct"]
        ctc = LE[k]["ctmc_pct"]
        ps, pc = PROTO_LEVERS[k]
        e = dict(label=LE[k]["label"], short=SHORT[k], kind=KIND[k], s12_pct=s12c, ctmc_pct=ctc,
                 s12_ci90_pct=ci(res_t[k][0]), ctmc_ci90_pct=ci(res_t[k][1]),
                 ctmc_ci90_usage_pct=ci(res_u[k][1]),
                 s12_median_pct=float(np.nanmedian(res_t[k][0])) if k != "wifi_share_80" else None,
                 ctmc_median_pct=float(np.median(res_t[k][1])),
                 p_gain_s12=float(np.mean(res_t[k][0] > 0)) if k != "wifi_share_80" else None,
                 p_gain_ctmc=float(np.mean(res_t[k][1] > 0)),
                 rank_ctmc_median=float(np.median(ranks[j])),
                 rank_ctmc_ci90=[float(v) for v in np.percentile(ranks[j], (5, 95))],
                 proto=dict(s12_pct=ps, ctmc_pct=pc),
                 dev_pp=dict(s12=None if ps is None else s12c - ps, ctmc=ctc - pc))
        if k in ids_s12:
            jj = ids_s12.index(k)
            e["rank_s12_median"] = float(np.median(ranks_s12[jj]))
            e["rank_s12_ci90"] = [float(v) for v in np.percentile(ranks_s12[jj], (5, 95))]
        if k in l2:
            e["s12_L2_tte_h"] = l2[k]
            e["s12_L2_pct"] = pct(l2[k], base_l2)
            e["s12_L2_proto_pct"] = PROTO_S12_L2.get(k)
        if k == "wifi_share_80":
            e["ctmc_mu_h"] = LE[k]["ctmc_mu_h"]
        levers[k] = e
        print(f"  {k:15s} S12 {s12c:+6.1f} {e['s12_ci90_pct']}  S13 {ctc:+6.1f} {e['ctmc_ci90_pct']} "
              f"(proto {ps}, {pc})  L2 S12 {e.get('s12_L2_pct', float('nan')):+.1f}")
    out["levers"] = levers
    out["s12_L2_base_h"] = base_l2
    out["s12_L2_note"] = ("L2 TTE of the S12 day with each lever (LSODA): non-smooth, the shutdown lands in the "
                          "evening of day 2. The package lever 'standby_5g' keeps each state's background rate "
                          "(the §12.1 smooth-metric definition); final_proto 'recs' used SB-5's own rate "
                          "(1/60 s^-1), so its -1.5 % is not the same lever.")
    order = sorted(ids, key=lambda k: -levers[k]["ctmc_pct"])
    out["ranking"] = dict(by="S13 Dynkin mean (central)", order=order,
                          p_skip_game_rank1=float(np.mean(ranks[ids.index("skip_game")] == 1)),
                          note="ranks over the parameter sample: 1 = largest TTE gain among all 13 items")

    # ---- energy shares of the S12 day by component (quasi-steady state powers)
    def shares(proto_method=False):
        tot = {k: 0.0 for k in BREAKDOWN_KEYS}
        tot["soc"] = 0.0
        for d, u in sch:
            Pns, dem, _, br = power_parts(u, RP63)
            if proto_method:  # final_proto 'recs': T_s from the unthrottled demand (not a fixed point)
                Ts = TAMB + 273.15 + float(RP63.k_env * RP63.R_skin) * float(Pns + dem)
            else:
                Ts = steady_skin_temperature(u, RP63, TAMB)
            for k in BREAKDOWN_KEYS:
                tot[k] += d * float(br[k])
            tot["soc"] += d * float(throttle(dem, Ts, RP63))
        S = sum(tot.values())
        return {k: v / S for k, v in sorted(tot.items(), key=lambda x: -x[1]) if v > 0}, S / sch.period

    sh, pbar_fp = shares()
    sh_p, pbar_p = shares(proto_method=True)
    out["energy_share_S12"] = dict(
        shares=sh, implied_Pbar_W=pbar_fp, proto_method_shares=sh_p, proto_method_implied_Pbar_W=pbar_p,
        proto=dict(disp=0.321, soc=0.294, net=0.170, awake=0.085, susp=0.048, gnss=0.040, aud=0.020, cam=0.015,
                   sens=0.007),
        note="'shares' use the SoC power at the F5 fixed point, consistent with the lever metric (P̄ = 0.4954 W). "
             "The §12.1 [PROTO] shares are reproduced by 'proto_method_shares', which takes the game's skin "
             "temperature from its unthrottled demand (T_s = 46.5 °C) and so throttles the game SoC to 0.58 W "
             "instead of 2.37 W; that implies P̄ = 0.421 W, inconsistent with the [PROTO] lever base 0.4954 W.")

    # ---- scenario-specific levers (L2 LSODA)
    def tte(u, Ta=TAMB, p=RP63, **kw):
        return run(Schedule.constant(u), Ta, p=p, **kw)

    def ratio_ci(ua, ub):
        """Parameter interval of the quasi-steady TTE gain E_u/P_a vs E_u/P_b (power ratio)."""
        return ci(100 * (steady_power(ub, pb, TAMB) / steady_power(ua, pb, TAMB) - 1))

    sl = {}
    r_a, r_c = tte("VW"), tte("V5")
    sl["video_wifi_vs_5g"] = dict(tte_wifi_h=r_a.tte_h, tte_5g_h=r_c.tte_h, gain_pct=pct(r_a.tte_h, r_c.tte_h),
                                  gain_ci90_pct=ratio_ci(PRIMITIVES["VW"], PRIMITIVES["V5"]),
                                  proto=dict(tte_wifi_h=32.32, tte_5g_h=21.10, gain_pct=53.0))
    r_a, r_c = tte("VC-W"), tte("VC")
    sl["vcall_wifi_vs_lte"] = dict(tte_wifi_h=r_a.tte_h, tte_lte_h=r_c.tte_h, gain_pct=pct(r_a.tte_h, r_c.tte_h),
                                   gain_ci90_pct=ratio_ci(PRIMITIVES["VC-W"], PRIMITIVES["VC"]),
                                   proto=dict(tte_wifi_h=10.83, tte_lte_h=6.93, gain_pct=56.0))
    r_b = tte(LEVERS["voice_call"][1]("VC", PRIMITIVES["VC"]))
    sl["voice_vs_video_call_lte"] = dict(tte_voice_h=r_b.tte_h, tte_video_h=r_c.tte_h,
                                         gain_pct=pct(r_b.tte_h, r_c.tte_h),
                                         gain_ci90_pct=ratio_ci(LEVERS["voice_call"][1]("VC", PRIMITIVES["VC"]),
                                                                PRIMITIVES["VC"]))
    warm = {}
    for soh in (1.0, 0.8):
        ref = run(s9_schedule(25.0), 25.0, T0_C=20.0, p=aged(RP63, 1 - soh))
        for T in (-10.0, -20.0):
            r = run(s9_schedule(T), T, T0_C=20.0, p=aged(RP63, 1 - soh))
            warm[f"{T:g}|{soh:.1f}"] = dict(s_off_cold_pct=100 * r.s_off, tte_cold_h=r.tte_h, event=r.event,
                                            s_off_warm_pct=100 * ref.s_off, tte_warm_h=ref.tte_h,
                                            stranded_avoided_pp=100 * (r.s_off - ref.s_off))
    sl["keep_warm_S9"] = dict(cells=warm, note="S9 cold hike (TRK + 30 s photo bursts), T0 = 20 °C; 'warm' = the same "
                                               "pack at 25 °C ambient", proto_s_off_range_pct=[7.2, 44.2])
    r25, r35 = tte("GM", 25.0), tte("GM", 35.0)
    f5 = A.thermal_cap("GM", RP63, 35.0)
    sl["hot_car_gaming"] = dict(tte_25C_h=r25.tte_h, tte_35C_h=r35.tte_h, rho_perf_25C_L2=r25.rho_perf,
                                rho_perf_35C_L2=r35.rho_perf, rho_perf_35C_F5=float(f5["rho_perf"]),
                                proto=dict(tte_25C_h=5.68, tte_35C_h=10.29, rho_perf_35C_F5=0.531))
    so85 = tte(replace(PRIMITIVES["SO-L"], RSRP=-85.0))
    so110 = tte(replace(PRIMITIVES["SO-L"], RSRP=-110.0))
    vc110 = tte(replace(PRIMITIVES["VC"], RSRP=-110.0))
    vc85 = tte(PRIMITIVES["VC"])
    sl["weak_signal"] = dict(social_lte_85_h=so85.tte_h, social_lte_110_h=so110.tte_h,
                             social_loss_pct=pct(so110.tte_h, so85.tte_h), vcall_85_h=vc85.tte_h,
                             vcall_110_h=vc110.tte_h, vcall_loss_pct=pct(vc110.tte_h, vc85.tte_h),
                             social_loss_ci90_pct=ratio_ci(replace(PRIMITIVES["SO-L"], RSRP=-110.0),
                                                           PRIMITIVES["SO-L"]),
                             proto=dict(social_loss_pct=-12.0, vcall_loss_pct=-23.0))
    sb = {k: tte(k).tte_h for k in ("SB-W", "SB-L", "SB-5")}
    sl["standby_rat"] = dict(tte_h=sb, wifi_vs_lte_pct=pct(sb["SB-W"], sb["SB-L"]),
                             lte_vs_5g_pct=pct(sb["SB-L"], sb["SB-5"]),
                             proto_tte_h={"SB-W": 442.7, "SB-L": 378.3, "SB-5": 173.8})
    out["scenario_levers"] = sl

    # ---- OS strategies
    os_: dict = {}
    pa = aged(RP63, 0.2)
    gov = {}
    for mb in (1.0, 0.5, 0.0):
        r = run(Schedule.constant("SO"), -10.0, p=replace(pa, m_b=mb))
        gov[f"{mb:g}"] = dict(s_off_pct=100 * r.s_off, tte_h=r.tte_h, event=r.event, burst_delivered_frac=mb)
    for v in gov.values():
        v["gain_tte_pct"] = pct(v["tte_h"], gov["1"]["tte_h"])
    r20a = run(Schedule.constant("SO"), -20.0, p=pa)
    r20b = run(Schedule.constant("SO"), -20.0, p=replace(pa, m_b=0.0))
    O1 = dict(mechanism="cap dPb <= max(0, 0.9 dPb_safe) with F4 evaluated online", m_b_sweep_m10C_soh0p8=gov,
              m20C_soh0p8=dict(s_off_pct=[100 * r20a.s_off, 100 * r20b.s_off], tte_h=[r20a.tte_h, r20b.tte_h],
                               note="m_b = 1 vs 0 (= the online governor's TTE and s_off, see AN12)"),
              proto=dict(s_off_pct=[8.5, 2.6, 1.5], tte_h=[12.15, 12.88, 13.00]))
    try:
        st = json.loads((RESULTS / "stranded.json").read_text())
        O1["online_governor_AN12"] = dict(st["governor"]["online_O1"], source="results/stranded.json (AN12)")
    except (OSError, KeyError, ValueError):
        O1["online_governor_AN12"] = None
    os_["O1"] = O1
    gm, gm64 = PRIMITIVES["GM"], PRIMITIVES["GM"].with_(x_c=0.64)
    r7 = tte(gm)
    r7c = tte(gm64)
    os_["O2"] = dict(mechanism="game x_c 0.80 -> 0.64 (-20 % frequency demand)",
                     steady_game_power_W=[float(steady_power(gm, RP63)), float(steady_power(gm64, RP63))],
                     s12_pct=levers["os_frame_cap"]["s12_pct"], s12_ci90_pct=levers["os_frame_cap"]["s12_ci90_pct"],
                     ctmc_pct=levers["os_frame_cap"]["ctmc_pct"], ctmc_ci90_pct=levers["os_frame_cap"]["ctmc_ci90_pct"],
                     S7_L2_tte_h=[r7.tte_h, r7c.tte_h], S7_L2_gain_pct=pct(r7c.tte_h, r7.tte_h),
                     S7_perf_vs_requested_0p80=[r7.rho_perf, r7c.rho_perf * 0.64 / 0.80],
                     note="performance = time-mean delivered x_eff / 0.80; reactive throttling oscillates around "
                          "its mean, the frame cap holds x_c = 0.64 (rho relative to 0.64 shown in S7_rho_capped)",
                     S7_rho_capped=r7c.rho_perf,
                     proto=dict(steady_game_power_W=[3.116, 2.722], s12_pct=3.4, ctmc_pct=2.7))
    co = LEVERS["os_coalesce"][1]
    stb = {}
    for k, pk in (("S1", "SB-W"), ("S2a", "SB-L"), ("S2b", "SB-5")):
        b, c_ = tte(PRIMITIVES[pk]).tte_h, tte(co(pk, PRIMITIVES[pk])).tte_h
        stb[k] = dict(base_h=b, coalesced_h=c_, gain_pct=pct(c_, b),
                      gain_ci90_pct=ratio_ci(co(pk, PRIMITIVES[pk]), PRIMITIVES[pk]))
    os_["O3"] = dict(mechanism="background burst rate lambda -> lambda/5 in standby states (E14)", standby=stb,
                     day_s12_pct=levers["os_coalesce"]["s12_pct"], day_ctmc_pct=levers["os_coalesce"]["ctmc_pct"],
                     day_ctmc_ci90_pct=levers["os_coalesce"]["ctmc_ci90_pct"],
                     proto=dict(standby_gain_pct=dict(S1=2.4, S2a=25.4, S2b=66.7), day_pct=[0.1, 0.7]))
    mu_nr = A.ctmc_mean_tte(CTMC("typical", rat="nr"), RP63)
    os_["O4"] = dict(mechanism="RAT steering: LTE or Wi-Fi for low-throughput states, 5G only for bulk",
                     avoided_standby_5g_ctmc_pct=levers["standby_5g"]["ctmc_pct"],
                     avoided_standby_5g_ci90_pct=levers["standby_5g"]["ctmc_ci90_pct"],
                     streaming_5g_vs_wifi_pct=pct(sl["video_wifi_vs_5g"]["tte_5g_h"],
                                                  sl["video_wifi_vs_5g"]["tte_wifi_h"]),
                     ctmc_all_cellular_on_5g_pct=pct(mu_nr, LE["base"]["ctmc_mu_h"]),
                     ctmc_all_cellular_on_5g_mu_h=mu_nr,
                     literature="XU20 measured -24.8 % (cited, not computed here)",
                     proto=dict(standby_5g_pct=-4.5, streaming_pct=-35.0))
    os_["O5"] = dict(mechanism="defer bulk transfers while RSRP < -100 dBm (m_sig)",
                     avoided_ctmc_pct=levers["weak_signal"]["ctmc_pct"],
                     avoided_ctmc_ci90_pct=levers["weak_signal"]["ctmc_ci90_pct"],
                     avoided_social_lte_pct=sl["weak_signal"]["social_loss_pct"],
                     proto=dict(ctmc_pct=-5.7, social_pct=-12.0))
    os_["O6"] = dict(mechanism="adaptive refresh (LTPO): 120 Hz only while scrolling; bound = 60 Hz cap",
                     s12_pct=levers["refresh_60"]["s12_pct"], ctmc_pct=levers["refresh_60"]["ctmc_pct"],
                     s12_ci90_pct=levers["refresh_60"]["s12_ci90_pct"],
                     ctmc_ci90_pct=levers["refresh_60"]["ctmc_ci90_pct"], proto=[4.1, 4.7])
    # O7: cold-aware gauge, F4 prediction of the stranded charge vs the L2 shutdown SOC (S10a cells)
    P_SO = float(demand_power(PRIMITIVES["SO"], RP63))
    o7 = {}
    for T in (0.0, -10.0, -20.0):
        for soh in (1.0, 0.8):
            pp = aged(RP63, 1 - soh)
            r = run(Schedule.constant("SO"), T, p=pp)
            f4 = A.brownout_root(P_SO, 3.0, r.T_end_C + 273.15, pp)
            o7[f"{T:g}|{soh:.1f}"] = dict(s_off_L2_pct=100 * r.s_off, event=r.event, T_end_C=r.T_end_C,
                                          s_off_F4_pct=100 * float(f4["s_off"]),
                                          err_pp=100 * (float(f4["s_off"]) - r.s_off))
    os_["O7"] = dict(mechanism="remap displayed 0 % to z*(T, SOH, dPb) from F4 (battery temperature measured); "
                               "show 'X % unavailable until warm'",
                     S10a_F4_vs_L2=o7, S9_s_off_pct={k: v["s_off_cold_pct"] for k, v in warm.items()},
                     note="F4 uses the battery temperature at shutdown (T_end of the L2 run) and dPb = 3 W; the S9 "
                          "shutdowns are explicit 30 s bursts against V_cut, which F4 (sub-second envelope) does not "
                          "cover - there the gauge needs the burst size")
    # O8: IG remaining-time band and its coverage
    ct = CTMC("typical")
    Lg = ct.generator()
    Pst = ct.state_powers(RP63, TAMB)
    Eu = float(cell.usable_energy(RP63))
    o8 = {}
    for frac in (1.0, 0.5):
        m = A.ctmc_moments(Lg, Pst, frac * Eu)
        q = A.ig_quantiles(m["mu_dynkin"], m["lam_IG"])
        tl0 = A.mc_L0(Lg, Pst, frac * Eu, n=20000, seed=1)
        o8[f"E_frac={frac:g}"] = dict(mu_dynkin_h=m["mu_dynkin"], ig_q05_q50_q95_h=[float(v) for v in q],
                                      coverage_L0=float(np.mean((tl0 >= q[0]) & (tl0 <= q[2]))),
                                      mc_L0_q05_q50_q95_h=[float(v) for v in np.percentile(tl0, (5, 50, 95))])
    t0 = time.time()
    b = run_batch(ct, TAMB, None, n=400, h=10.0, seed=7, t_max_h=120.0)
    tl2 = b["tte_h"][np.isfinite(b["tte_h"])]
    q = o8["E_frac=1"]["ig_q05_q50_q95_h"]
    ratio = float(np.mean(tl2) / o8["E_frac=1"]["mu_dynkin_h"])
    o8["L2_400_paths"] = dict(n=int(tl2.size), mean_h=float(np.mean(tl2)),
                              coverage_raw_IG=float(np.mean((tl2 >= q[0]) & (tl2 <= q[2]))),
                              L2_over_dynkin=ratio,
                              coverage_IG_scaled=float(np.mean((tl2 >= ratio * q[0]) & (tl2 <= ratio * q[2]))),
                              note="L2 Monte Carlo (h = 10 s, seed 7); 'scaled' multiplies the IG quantiles by the "
                                   "in-sample L2/Dynkin mean ratio (optimistic)", runtime_s=time.time() - t0)
    os_["O8"] = dict(mechanism="EKF on (z, v1) + online P̄, σ² and the IG law (F6-F7) -> remaining time with a "
                               "90 % band", S13_typical=o8,
                     display_example=f"{q[1]:.0f} h ({q[0]:.0f}-{q[2]:.0f} h) from a full charge at 07:00")
    # O9: charging defaults (quasi-static aging driver; central values; UQ from results/aging.json if present)
    pcal = replace(RP63, **dict(zip(("k_H", "k_L"), G.calibrate(RP63))))
    soh2 = {sid: float(G.run_years(years=2, p=pcal, **G.AGING_SCENARIOS[sid])["soh"][729])
            for sid in ("A1", "A2", "A3", "A4", "A5", "A6b", "A6c")}
    O9 = dict(mechanism="charging defaults (E17-E21)", soh_2yr=soh2,
              limit80_pp=100 * (soh2["A2"] - soh2["A1"]), obc_pp=100 * (soh2["A3"] - soh2["A1"]),
              hot_climate_pp=100 * (soh2["A4"] - soh2["A1"]), car_dashboard_pp=100 * (soh2["A5"] - soh2["A1"]),
              jeita_cold_pp=100 * (soh2["A6c"] - soh2["A6b"]),
              proto=dict(limit80_pp=1.9, obc_pp=0.7, hot_climate_pp=-12.5, jeita_cold_pp=7.0))
    try:
        ag = json.loads((RESULTS / "aging.json").read_text())
        O9["uq_AN13"] = dict(diff_vs_A1_2yr_pp=ag["uq"]["diff_vs_A1_2yr_pp"],
                             diff_A6c_vs_A6b_2yr_pp=ag["uq"]["diff_A6c_vs_A6b_2yr_pp"], n=ag["uq"]["n"],
                             source="results/aging.json (AN13)")
    except (OSError, KeyError, ValueError):
        O9["uq_AN13"] = None
    os_["O9"] = O9
    out["os"] = os_

    runtime = time.time() - t_start
    out = dict(meta=meta("run_recommendations.py", seed=SEED, solver="smooth metrics (F5 fixed points, F6 Dynkin); "
                                                                     "LSODA for L2 TTE; rk4_batch h = 10 s for the "
                                                                     "O8 L2 coverage", h_s=10.0,
                         runtime_s=round(runtime, 1), n_theta=N_THETA, Tamb_C=TAMB),
               **out)
    write_json("recommendations.json", out)
    make_figure(json.loads((RESULTS / "recommendations.json").read_text()))
    print(f"done in {time.time() - t_start:.1f} s")


# ------------------------------------------------------------------------------------------- figure
def make_figure(o):
    apply()
    lv = o["levers"]
    order = sorted(lv, key=lambda k: lv[k]["ctmc_pct"])  # bottom -> top = smallest -> largest
    y = np.arange(len(order))
    col = {"gain": BLUE, "os": SEQ_BLUE[3], "avoid": RED}
    fig, axes = plt.subplots(1, 2, figsize=(W1, 3.0), sharey=True, gridspec_kw=dict(wspace=0.08))
    for ax, key, title in ((axes[0], "s12", "S12 scheduled day, $E_u/\\bar P$"),
                           (axes[1], "ctmc", "S13 typical day (CTMC), Dynkin mean")):
        for i, k in enumerate(order):
            e = lv[k]
            v = e[f"{key}_pct"]
            if v is None or not np.isfinite(v):
                ax.text(0.6, i, "n/a (no schedule)", va="center", ha="left", fontsize=6.3, color=MUTED)
                continue
            lo, hi = e[f"{key}_ci90_pct"]
            ax.barh(i, v, height=0.62, color=col[e["kind"]], zorder=2)
            ax.plot([lo, hi], [i, i], color=INK2, lw=0.8, zorder=3)
            ax.plot([lo, lo], [i - 0.17, i + 0.17], color=INK2, lw=0.8, zorder=3)
            ax.plot([hi, hi], [i - 0.17, i + 0.17], color=INK2, lw=0.8, zorder=3)
            med = e.get(f"{key}_median_pct")
            if med is not None:
                ax.plot([med], [i], marker="o", ms=3.0, mfc="white", mec=INK, mew=0.7, zorder=4)
            if v >= 0:
                ax.text(max(hi, v) + 0.8, i, f"+{v:.1f}", va="center", ha="left", fontsize=6.3, color=INK)
            else:
                ax.text(min(lo, v) - 0.8, i, f"{v:.1f}".replace("-", "−"), va="center", ha="right",
                        fontsize=6.3, color=INK)
        ax.axvline(0, color=INK2, lw=0.8, zorder=2)
        ax.set_xlim(-32, 45)
        ax.set_xticks([-20, -10, 0, 10, 20, 30, 40])
        ax.set_xticklabels(["−20", "−10", "0", "+10", "+20", "+30", "+40"])
        ax.grid(axis="y", visible=False)
        ax.set_xlabel("Change in time-to-empty (%)")
        ax.set_title(title, loc="left", fontsize=8.5)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([lv[k]["short"] for k in order], fontsize=7.2)
    axes[0].tick_params(axis="y", length=0)
    axes[1].tick_params(axis="y", length=0)
    panel_label(axes[0], "(a)")
    panel_label(axes[1], "(b)")
    from matplotlib.patches import Patch
    handles = [Patch(color=BLUE, label="user lever"), Patch(color=SEQ_BLUE[3], label="OS strategy"),
               Patch(color=RED, label="habit to avoid"),
               plt.Line2D([], [], color=INK2, lw=0.8, label="90 % param. interval"),
               plt.Line2D([], [], color=INK, lw=0, marker="o", ms=3.0, mfc="white", mew=0.7,
                          label="param. median")]
    axes[0].legend(handles=handles, loc="lower right", fontsize=6.3, handlelength=1.2, borderaxespad=0.15,
                   labelspacing=0.3)
    save(fig, "fig_levers")


if __name__ == "__main__":
    main()
