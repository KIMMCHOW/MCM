"""AN7: uncertainty quantification (FINAL_SPEC §9) -> results/uq.json, fig_uq_violins, tab_tte_uq.tex (T3).

What is computed
----------------
1. Deterministic scenarios (§9 item 1).  One joint scrambled-Sobol' sample of the §4 parameters
   (``UQ_SPEC``, N = 2048, seed 2026) is pushed through the L2 model with ``model.run_batch`` (rk4_batch):
   S1–S9, S10a (S4 at SOH 0.8), S10b (S12 at SOH 0.8), S11a/S11c (−110 dBm, 25 °C only), S12, S16 at
   T_amb = 25 / 0 / −10 °C from s0 = 100 %, and the S15 set (S4, S5a, S7, S12) from s0 = 50 / 20 % at the same
   three temperatures.  Extra cells: S9 at −10 °C and −20 °C on an SOH-0.8 pack (the §5.4 stranded-charge cases).
   Per cell: median, 5–95 %, mean, SD, CV, Pr(TTE < t*), s_off quantiles, event shares, and the default-θ
   LSODA value ("nominal").  t* = 16 h for the S12 day from 100 % (§9 item 5, "dies before 23:00"); otherwise
   t* = 0.9 × nominal (probability of falling more than 10 % short of the point prediction).
   Compute subset (stated deviation from §9, see ``deviations`` in the JSON):
   * constant-load scenarios use h = 10 s (no input switches; checked against h = 2 s on the first 256 Sobol'
     points of five cells, ``step_check``); S9, S12, S10b (switching schedules) use the §9 default h = 2 s;
     standby S1/S2 use h = 10 s as §9 allows.
   * standby S1/S2a/S2b at 0 and −10 °C use the first 512 points (a balanced Sobol' subset) instead of 2048.
2. Delta-method check (§9 item 2): CV_TTE² ≈ CV_Eu² + Σ s_i² CV_Pi² (L0, quasi-steady component powers) vs the
   exact L0 sample E_u(θ)/P̄(θ) and the L2 QMC CV, for every s0 = 100 % cell.
3. Stochastic scenarios S13/S14 (§9 item 3): law of total variance with the inner moments from F6–F7 per θ
   sample (μ_Dynkin, σ²_CLT·E_u/P̄³), N_θ = 2048: (A) hardware θ only (the §9 statement), (B) usage-model
   factors only (§4.6 USAGE_UQ_SPEC, independent pseudo-random stream, seed 2027), (C) both, and a pick-freeze
   group decomposition of (C) into usage path / hardware / usage model / interaction.  Check on 50 θ × 200 paths of
   full L2 (h = 10 s): per-θ mean within 3 % and SD within 10 % of F6/F7 after applying the L2/L0 mean ratio
   (pooled over the check, measured on the same paths).
4. Bonn-constrained posterior (§9 item 4, ABC): the same 2048 θ on the iPhone 14 Pro black-screen idle case
   (Q = 3.32 Ah measured, R·Q scaling), accept |TTE/55.03 h − 1| ≤ 5 %; prior vs posterior of the floor
   parameters and the narrowing of every 25 °C scenario interval and of the S13/S14 mixture.

Outputs: results/uq.json with the §14 keys ``"{S}|s0=..|Tamb=..[|soh=..]".{median_h, q05_h, q95_h, cv,
p_lt_target, s_off_q05_q50_q95}`` (s_off as a fraction; ``..._pct`` duplicates in %), ``variance_split.{S13,S14}
.{usage_frac, param_frac}``, ``delta_method_check``, ``abc_bonn``; the digest for the paper writer is results/uq_sens_summary.md.
Per-cell raw samples are cached in results/cache/uq/*.npz (``--fresh`` recomputes everything).

Run:  cd /home/user/MCM && PYTHONPATH=src:analysis python3 analysis/run_uq.py [--workers 3] [--fresh]
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import math
import multiprocessing as mp
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import lru_cache

import numpy as np

from basics_common import RESULTS, meta, write_json, write_tex
from phonebatt import cell
from phonebatt.analytic import ctmc_moments
from phonebatt.model import run, run_batch
from phonebatt.params import RP63, UQ_SPEC, USAGE_UQ_SPEC, ACTIVITIES, iphone14pro, sample, with_theta
from phonebatt.power import demand_power, power_parts, steady_power
from phonebatt.usage import BONN_CASES, CTMC, PROFILES, Schedule, scenario
from phonebatt.validation import bonn_observations

SEED = 2026
SEED_USAGE = 2027
SEED_B, SEED_USAGE_B = 2028, 2029  # second, independent samples for the pick-freeze variance decomposition
N_THETA = 2048
N_COLD_STANDBY = 512
N_STEP_CHECK = 256
CACHE = RESULTS / "cache" / "uq"
CONST = ("S1", "S2a", "S2b", "S3", "S4", "S5a", "S5b", "S5c", "S6", "S7", "S8a", "S8b", "S16")
STANDBY = ("S1", "S2a", "S2b")
S15_SET = ("S4", "S5a", "S7", "S12")
TEMPS = (25.0, 0.0, -10.0)
LABELS = {"S1": "S1 standby, home Wi-Fi", "S2a": "S2a standby, LTE", "S2b": "S2b standby, 5G heavy push",
          "S3": "S3 reading", "S4": "S4 social, Wi-Fi", "S5a": "S5a video, Wi-Fi", "S5b": "S5b video, LTE",
          "S5c": "S5c video, 5G NSA", "S6": "S6 navigation", "S7": "S7 3D gaming", "S8a": "S8a video call, LTE",
          "S8b": "S8b video call, Wi-Fi", "S9": "S9 cold hike (T0 = 20 C)", "S10a": "S10a social, SOH 0.8",
          "S10b": "S10b day, SOH 0.8", "S11a": "S11a social LTE, -110 dBm", "S11c": "S11c video call, -110 dBm",
          "S12": "S12 scheduled day", "S16": "S16 audio, screen off"}
CTMC_CHECK = dict(n_theta=50, n_paths=200, chunk_theta=10, h=10.0, t_max_h=150.0)
STEP_CHECK_CELLS = (("S4", 25.0, 1.0, 1.0), ("S7", 25.0, 1.0, 1.0), ("S10a", -10.0, 1.0, 0.8),
                    ("S6", -10.0, 1.0, 1.0), ("S8a", 0.0, 1.0, 1.0))
CODES = ("none", "charge", "cutoff", "brownout", "thermal", "infeasible")


# ------------------------------------------------------------------------------------------- samples
@lru_cache(maxsize=1)
def theta_all() -> dict:
    """Joint §4 parameter sample (scrambled Sobol', seed 2026, N = 2048), dict name -> (N,)."""
    return sample(UQ_SPEC, N_THETA, method="sobol", seed=SEED)


@lru_cache(maxsize=1)
def usage_all() -> dict:
    """§4.6 usage-model factors: multipliers (N,) and Dirichlet p (7, N), paired by index with ``theta_all``.

    Drawn from an independent pseudo-random stream (seed 2027), not from a second scrambled Sobol' sequence:
    two scrambled Sobol' sequences paired by index are deterministically coupled dimension by dimension (the
    same radical inverse of the index underlies both), so factor j of one would be a function of factor j of
    the other (2-D histogram chi-square 14336 vs about 63 for independent pairs on 8 x 8 bins, N = 2048).
    """
    return sample(USAGE_UQ_SPEC, N_THETA, method="random", seed=SEED_USAGE)


def theta_first(n: int, idx=None) -> dict:
    th = theta_all()
    if idx is not None:
        return {k: v[idx] for k, v in th.items()}
    return {k: v[:n] for k, v in th.items()}


# ------------------------------------------------------------------------------------------- cells
def make_cell(sid, Ta, s0, soh, h, N, kind="det"):
    key = f"{sid}|s0={int(round(100 * s0))}|Tamb={Ta:g}" + ("" if soh == 1.0 else f"|soh={soh:g}")
    return dict(kind=kind, key=key, sid=sid, Tamb=float(Ta), s0=float(s0), soh=float(soh), h=float(h), N=int(N))


def cell_list() -> list[dict]:
    out = []
    for Ta in TEMPS:
        for sid in CONST:
            N = N_COLD_STANDBY if (sid in STANDBY and Ta != 25.0) else N_THETA
            out.append(make_cell(sid, Ta, 1.0, 1.0, 10.0, N))
        out.append(make_cell("S9", Ta, 1.0, 1.0, 2.0, N_THETA))
        out.append(make_cell("S12", Ta, 1.0, 1.0, 2.0, N_THETA))
        out.append(make_cell("S10a", Ta, 1.0, 0.8, 10.0, N_THETA))
        out.append(make_cell("S10b", Ta, 1.0, 0.8, 2.0, N_THETA))
        for s0 in (0.5, 0.2):
            for sid in S15_SET:
                out.append(make_cell(sid, Ta, s0, 1.0, 2.0 if sid == "S12" else 10.0, N_THETA))
    out.append(make_cell("S9", -10.0, 1.0, 0.8, 2.0, N_THETA))
    out.append(make_cell("S9", -20.0, 1.0, 0.8, 2.0, N_THETA))
    out.append(make_cell("S11a", 25.0, 1.0, 1.0, 10.0, N_THETA))
    out.append(make_cell("S11c", 25.0, 1.0, 1.0, 10.0, N_THETA))
    return out


def build(c: dict):
    """Scenario object of a cell (schedule, base parameters, T0)."""
    sid, Ta, s0, soh = c["sid"], c["Tamb"], c["s0"], c["soh"]
    if sid == "S10a":
        return scenario("S4", Tamb_C=Ta, s0=s0, soh=soh)
    if sid == "S10b":
        return scenario("S12", Tamb_C=Ta, s0=s0, soh=soh)
    if sid == "S11a":
        return scenario("S11", Tamb_C=Ta, s0=s0, soh=soh, variant="a", RSRP=-110.0)
    if sid == "S11c":
        return scenario("S11", Tamb_C=Ta, s0=s0, soh=soh, variant="c", RSRP=-110.0)
    return scenario(sid, Tamb_C=Ta, s0=s0, soh=soh)


def est_cost(c: dict) -> float:
    """Rough relative cost (simulated hours × samples / step) for longest-first scheduling."""
    if c["kind"] == "ctmc":
        return 25.0 * (40.0 if c["profile"] == "typical" else 22.0)
    if c["kind"] == "bonn":
        return 55.0
    sc = build(c)
    P = sc.schedule.mean(lambda u: float(demand_power(u, sc.p)))
    tte = float(cell.usable_energy(sc.p)) * c["s0"] / P
    return tte * (c["N"] / N_THETA) ** 0.7 * (4.0 if c["h"] == 2.0 else 1.0)


def cache_path(c: dict):
    name = c["kind"] + "_" + c["key"].replace("|", "__").replace("=", "-").replace(".", "p")
    return CACHE / f"{name}.npz"


def _codes(ev) -> np.ndarray:
    return np.array([CODES.index(str(e)) for e in ev], dtype=np.int8)


# ------------------------------------------------------------------------------------------- workers
def job_det(c: dict) -> dict:
    """One deterministic cell: run_batch over the first N θ samples + the default-θ LSODA reference."""
    t0 = time.time()
    sc = build(c)
    th = theta_first(c["N"])
    r = run_batch(sc.schedule, c["Tamb"], th, s0=c["s0"], h=c["h"], T0_C=sc.T0_C, p=sc.p)
    t1 = time.time()
    nom = run(sc.schedule, c["Tamb"], s0=c["s0"], T0_C=sc.T0_C, p=sc.p)
    out = dict(tte_h=r["tte_h"], event=_codes(r["event"]), s_off=r["s_off"], soe_off=r["soe_off"],
               T_end_C=r["T_end_C"], nom_tte_h=nom.tte_h, nom_event=CODES.index(nom.event), nom_s_off=nom.s_off,
               nom_soe_off=nom.soe_off, nom_P_mean_W=nom.P_mean_W, runtime_s=t1 - t0, h=c["h"], N=c["N"])
    return out


def job_step(c: dict) -> dict:
    """h = 2 s reference on the first 256 θ of a constant-load cell (compared with the h = 10 s run)."""
    t0 = time.time()
    sc = build(c)
    th = theta_first(N_STEP_CHECK)
    r = run_batch(sc.schedule, c["Tamb"], th, s0=c["s0"], h=2.0, T0_C=sc.T0_C, p=sc.p)
    return dict(tte_h=r["tte_h"], event=_codes(r["event"]), s_off=r["s_off"], runtime_s=time.time() - t0)


def bonn_theta(th: dict) -> dict:
    """§4.7 device scaling applied to θ: Q fixed at the measured 3.32 Ah, R0/R1 scaled by 4.70/3.32 (A12)."""
    q = 3.32
    out = {k: v for k, v in th.items() if k != "Q_n"}
    out["R0_ref"] = th["R0_ref"] * 4.70 / q
    out["R1_ref"] = th["R1_ref"] * 4.70 / q
    return out


def job_bonn(c: dict) -> dict:
    t0 = time.time()
    th = bonn_theta(theta_first(N_THETA))
    p14 = iphone14pro()
    sch = Schedule.constant(BONN_CASES["idle"])
    r = run_batch(sch, 25.0, th, h=10.0, p=p14)
    nom = run(sch, 25.0, p=p14)
    return dict(tte_h=r["tte_h"], event=_codes(r["event"]), nom_tte_h=nom.tte_h, runtime_s=time.time() - t0)


def l0_on_paths(paths, P_sp: np.ndarray, Eu: np.ndarray) -> np.ndarray:
    """Exact L0 TTE [h] on stored CTMC paths with per-path state powers P_sp (16, n) and energies Eu (n,)."""
    t = paths.times / 3600.0
    st = paths.states
    with np.errstate(invalid="ignore"):
        dt = np.diff(t, axis=1)
    dt = np.where(np.isfinite(dt), dt, 0.0)
    rows = np.arange(st.shape[0])
    Pp = P_sp[st, rows[:, None]]  # (n, M)
    cum = np.cumsum(Pp[:, :-1] * dt, axis=1)
    hit = cum >= Eu[:, None]
    any_hit = hit.any(axis=1)
    i = np.argmax(hit, axis=1)
    before = np.where(i > 0, cum[rows, np.maximum(i - 1, 0)], 0.0)
    out = t[rows, i] + (Eu - before) / Pp[rows, i]
    return np.where(any_hit, out, np.nan)


def job_ctmc(c: dict) -> dict:
    """Full-L2 check chunk: θ indices c['idx'] × n_paths independent CTMC paths (h = 10 s) + L0 on the same paths."""
    t0 = time.time()
    cfg = CTMC_CHECK
    idx = np.asarray(c["idx"])
    m = cfg["n_paths"]
    th = {k: np.repeat(v, m) for k, v in theta_first(0, idx).items()}
    n = idx.size * m
    ct = CTMC(c["profile"])
    paths = ct.sample_paths(n, cfg["t_max_h"], h=cfg["h"], seed=c["seed"])
    r = run_batch(ct, 25.0, th, h=cfg["h"], paths=paths, t_max_h=cfg["t_max_h"])
    pb = with_theta(RP63, theta_first(0, idx))
    P_sp = np.broadcast_to(ct.state_powers(pb, 25.0, capped=True), (16, idx.size))
    Eu = np.broadcast_to(cell.usable_energy(pb), (idx.size,))
    l0 = l0_on_paths(paths, np.repeat(P_sp, m, axis=1), np.repeat(Eu, m))
    return dict(tte_h=r["tte_h"], event=_codes(r["event"]), s_off=r["s_off"], l0_h=l0, idx=idx,
                runtime_s=time.time() - t0)


JOBS = {"det": job_det, "step": job_step, "bonn": job_bonn, "ctmc": job_ctmc}


def run_job(c: dict) -> tuple[str, dict]:
    path = cache_path(c)
    out = JOBS[c["kind"]](c)
    np.savez(path, **{k: np.asarray(v) for k, v in out.items()})
    return c["key"], out


def load_cached(c: dict) -> dict | None:
    path = cache_path(c)
    if not path.exists():
        return None
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def all_jobs() -> list[dict]:
    jobs = cell_list()
    for sid, Ta, s0, soh in STEP_CHECK_CELLS:
        jobs.append(make_cell(sid, Ta, s0, soh, 2.0, N_STEP_CHECK, kind="step"))
    jobs.append(dict(kind="bonn", key="bonn_14pro_idle"))
    cfg = CTMC_CHECK
    for pi, prof in enumerate(("typical", "heavy")):
        for ch in range(cfg["n_theta"] // cfg["chunk_theta"]):
            idx = list(range(ch * cfg["chunk_theta"], (ch + 1) * cfg["chunk_theta"]))
            jobs.append(dict(kind="ctmc", key=f"ctmc_{prof}_{ch}", profile=prof, idx=idx,
                             seed=SEED + 100 * pi + ch))
    return jobs


def compute(workers: int, only: str | None = None) -> dict:
    """Run every job not in the cache (longest first, process pool) and return {key: result}.

    ``only``: regular expression on "kind:key"; jobs that do not match are not run (cached results are still
    loaded), so a full run can be split into several shorter invocations.
    """
    import re
    CACHE.mkdir(parents=True, exist_ok=True)
    jobs = all_jobs()
    res, todo = {}, []
    for c in jobs:
        r = load_cached(c)
        if r is None:
            if only is None or re.search(only, c["kind"] + ":" + c["key"]):
                todo.append(c)
        else:
            res[c["kind"] + ":" + c["key"]] = r
    if todo:
        todo.sort(key=est_cost, reverse=True)
        print(f"[uq] {len(todo)} jobs to run, {len(res)} cached, {len(jobs) - len(todo) - len(res)} skipped "
              f"(--only), {workers} workers", flush=True)
        t0 = time.time()
        if workers <= 1:
            for c in todo:
                _, r = run_job(c)
                res[c["kind"] + ":" + c["key"]] = r
                print(f"[uq] {c['kind']:5s} {c['key']:40s} {r['runtime_s']:7.1f} s  (t = {time.time() - t0:.0f} s)",
                      flush=True)
        else:
            with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("fork")) as ex:
                futs = {ex.submit(run_job, c): c for c in todo}
                for f in as_completed(futs):
                    c = futs[f]
                    _, r = f.result()
                    res[c["kind"] + ":" + c["key"]] = r
                    print(f"[uq] {c['kind']:5s} {c['key']:40s} {float(r['runtime_s']):7.1f} s  "
                          f"(t = {time.time() - t0:.0f} s)", flush=True)
    return res


# ------------------------------------------------------------------------------------------- statistics
def qs(x, q=(5, 50, 95)):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return [float(v) for v in np.percentile(x, q)] if x.size else [math.nan] * len(q)


def cell_summary(c: dict, r: dict) -> dict:
    tte = np.asarray(r["tte_h"], float)
    ev = np.asarray(r["event"])
    fin = np.isfinite(tte)
    x = tte[fin]
    nom = float(r["nom_tte_h"])
    if c["sid"] in ("S12", "S10b") and c["s0"] == 1.0:
        target, rule = 16.0, "16 h: an S12 day started at 07:00 that dies before 23:00 (FINAL_SPEC §9 item 5)"
    else:
        target, rule = 0.9 * nom, "0.9 x nominal (default-theta LSODA) TTE"
    so = np.asarray(r["s_off"], float)
    soe = np.asarray(r["soe_off"], float)
    q05, q50, q95 = qs(x)
    shares = {CODES[k]: float(np.mean(ev == k)) for k in range(len(CODES)) if np.any(ev == k)}
    return dict(
        scenario=c["sid"], label=LABELS[c["sid"]], Tamb_C=c["Tamb"], s0_pct=100 * c["s0"], soh=c["soh"],
        n=int(tte.size), n_no_event=int((~fin).sum()), n_dead_at_start=int(np.sum(tte == 0.0)), h_s=c["h"],
        median_h=q50, q05_h=q05, q95_h=q95, mean_h=float(x.mean()), sd_h=float(x.std(ddof=1)),
        cv=float(x.std(ddof=1) / x.mean()) if x.mean() > 0 else math.nan,
        rel_width90=(q95 - q05) / q50 if q50 > 0 else math.nan,
        p_lt_target=float(np.mean(tte < target)), target_h=target, target_rule=rule,
        p_lt_0p9_nominal=float(np.mean(tte < 0.9 * nom)),
        nominal_h=nom, nominal_event=CODES[int(r["nom_event"])], nominal_s_off_pct=100 * float(r["nom_s_off"]),
        median_over_nominal=q50 / nom if nom > 0 else math.nan, p_dead_at_start=float(np.mean(tte == 0.0)),
        s_off_q05_q50_q95=qs(so), s_off_q05_q50_q95_pct=[100 * v for v in qs(so)],
        s_off_mean_pct=100 * float(np.nanmean(so)), soe_off_q05_q50_q95_pct=[100 * v for v in qs(soe)],
        p_s_off_gt_5pct=float(np.mean(so > 0.05)), event_share=shares,
        T_end_q05_q50_q95_C=qs(r["T_end_C"]), runtime_s=float(r["runtime_s"]))


# ------------------------------------------------------------------------------------------- delta method (§9.2)
COMP_GROUPS = {"floor": ("susp", "awake", "wake"), "display": ("disp",), "radio": ("net",), "gnss_sens": ("gnss", "sens"),
               "camera": ("cam",), "audio_bt": ("aud",), "extra": ("extra",)}


def components(sch: Schedule, p, Tamb_C: float) -> dict:
    """Time-averaged quasi-steady component powers of a schedule (SoC at the F5 fixed point) [W]."""
    acc: dict[str, np.ndarray] = {}
    for d, u in sch:
        Pns, dem, _, br = power_parts(u, p)
        Pst = steady_power(u, p, Tamb_C)
        parts = {g: sum(np.asarray(br[k], float) for k in ks) for g, ks in COMP_GROUPS.items()}
        parts["soc"] = np.asarray(Pst - Pns, float)
        for g, v in parts.items():
            acc[g] = acc.get(g, 0.0) + d * v
    return {g: v / sch.period for g, v in acc.items()}


def delta_method(c: dict, summ: dict) -> dict:
    sc = build(c)
    pn = sc.p
    pb = with_theta(sc.p, theta_first(c["N"]))
    z0n = cell.z_from_s0(c["s0"], pn)
    Eu_n = float(cell.usable_energy(pn, z0n))
    Eu = np.asarray(cell.usable_energy(pb, cell.z_from_s0(c["s0"], pb)), float) * np.ones(c["N"])
    Cn = {g: float(np.asarray(v)) for g, v in components(sc.schedule, pn, c["Tamb"]).items()}
    Cb = {g: np.asarray(v, float) * np.ones(c["N"]) for g, v in components(sc.schedule, pb, c["Tamb"]).items()}
    Pn = sum(Cn.values())
    Pb = sum(Cb.values())
    cv = lambda a: float(np.std(a, ddof=1) / np.mean(a)) if np.mean(a) > 0 else 0.0
    shares = {g: Cn[g] / Pn for g in Cn}
    terms = {g: (shares[g] * cv(Cb[g])) ** 2 for g in Cn}
    cv_eu = cv(Eu)
    cv_delta = math.sqrt(cv_eu ** 2 + sum(terms.values()))
    l0 = Eu / Pb
    return dict(cv_delta=cv_delta, cv_L0_sample=cv(l0), cv_L2=summ["cv"], ratio_L2_to_delta=summ["cv"] / cv_delta,
                cv_Eu=cv_eu, energy_share=shares, cv_component={g: cv(Cb[g]) for g in Cb},
                variance_term=terms, top_term=max(terms, key=terms.get), L0_nominal_h=Eu_n / Pn,
                L2_nominal_h=summ["nominal_h"], L0_median_h=float(np.median(l0)))


# ------------------------------------------------------------------------------------------- CTMC split (§9.3)
def ctmc_inner(profile: str, theta: dict | None, usage: dict | None, n: int) -> dict:
    """Per-sample F6/F7 inner moments: μ_Dynkin(θ), Var(TTE|θ) = σ²E_u/P̄³ (L0, 25 °C, start in (idle, Wi-Fi))."""
    lam0, p0, d0 = PROFILES[profile]
    pb = with_theta(RP63, theta) if theta is not None else RP63
    P_sp = np.broadcast_to(np.asarray(CTMC(profile).state_powers(pb, 25.0, capped=True), float).reshape(16, -1),
                           (16, n))
    Eu = np.broadcast_to(np.asarray(cell.usable_energy(pb), float), (n,))
    mu, var, pbar = np.empty(n), np.empty(n), np.empty(n)
    L_nom = CTMC(profile).generator()
    for i in range(n):
        if usage is None:
            L = L_nom
        else:
            pv = usage["p_" + profile][:, i]
            ct = CTMC(profile, lam_on=lam0 * usage["lam_on_mul"][i],
                      p={a: float(pv[j]) for j, a in enumerate(ACTIVITIES)},
                      d_min={a: d0[a] * usage[f"d_mul_{a}"][i] for a in ACTIVITIES},
                      mean_w=2.5 * usage["conn_mul_wifi"][i], mean_c=2.5 * usage["conn_mul_cell"][i])
            L = ct.generator()
        m = ctmc_moments(L, P_sp[:, i], float(Eu[i]))
        mu[i], var[i], pbar[i] = m["mu_dynkin"], m["sd_clt"] ** 2, m["Pbar"]
    return dict(mu=mu, var=var, Pbar=pbar)


def mixture_quantiles(mu, var, q=(5, 50, 95), t_star=16.0, k=64, seed=11):
    """Quantiles and Pr(TTE < t*) of the normal mixture Σ_θ N(μ_θ, Var_θ)/N (AN5: normal fits better than IG)."""
    rng = np.random.default_rng(seed)
    z = rng.standard_normal((mu.size, k))
    x = (mu[:, None] + np.sqrt(var)[:, None] * z).ravel()
    from scipy.stats import norm
    p = float(np.mean(norm.cdf((t_star - mu) / np.sqrt(var))))
    return [float(v) for v in np.percentile(x, q)], p


def split(inner: dict) -> dict:
    mu, var = inner["mu"], inner["var"]
    Ev, Vm = float(np.mean(var)), float(np.var(mu, ddof=1))
    tot = Ev + Vm
    qq, p16 = mixture_quantiles(mu, var)
    return dict(total_var_h2=tot, E_var_usage_h2=Ev, var_E_param_h2=Vm, usage_frac=Ev / tot, param_frac=Vm / tot,
                mean_h=float(np.mean(mu)), sd_total_h=math.sqrt(tot), cv_total=math.sqrt(tot) / float(np.mean(mu)),
                mu_q05_q50_q95_h=qs(mu), sd_path_median_h=float(np.median(np.sqrt(var))),
                mixture_q05_q50_q95_h=qq, mixture_p_lt_16h=p16, Pbar_q05_q50_q95_W=qs(inner["Pbar"]))


def group_decomposition(fA: dict, fB: dict, f_hwB: dict, f_usB: dict, n_boot: int = 500, seed: int = 5) -> dict:
    """Pick-freeze split of Var(TTE) into usage path / hardware θ / usage-model θ / interaction.

    A = (hw, us), B = (hw', us') independent samples; A_B^hw = (hw', us), A_B^us = (hw, us').  On the inner mean
    μ = E[TTE | θ] the total effects VT_g = mean((f_A − f_{A_B^g})²)/2 (Jansen 1999) are precise; with two groups
    the first-order effects follow exactly as V_hw = V − VT_us and V_us = V − VT_hw, and the interaction is
    VT_hw + VT_us − V.  The direct Saltelli (2010) first-order estimator mean(f_B (f_{A_B^g} − f_A)) is reported as
    a cross-check (it is noisy here because the mean of μ is large compared with its SD).  The path term is
    E_θ[Var(TTE | θ)] over A and B.  95 % bootstrap CIs (resampling rows).
    """
    a, b, ah, au = (np.asarray(f["mu"], float) for f in (fA, fB, f_hwB, f_usB))
    ev = np.concatenate([fA["var"], fB["var"]])

    def est(ix):
        V = np.var(np.concatenate([a[ix], b[ix]]), ddof=1)
        E = float(np.mean(np.concatenate([fA["var"][ix], fB["var"][ix]])))
        VTh = 0.5 * np.mean((a[ix] - ah[ix]) ** 2)
        VTu = 0.5 * np.mean((a[ix] - au[ix]) ** 2)
        Vh_s = np.mean(b[ix] * (ah[ix] - a[ix]))
        Vu_s = np.mean(b[ix] * (au[ix] - a[ix]))
        T = E + V
        return np.array([E / T, (V - VTu) / T, (V - VTh) / T, (VTh + VTu - V) / T, VTh / T, VTu / T,
                         Vh_s / T, Vu_s / T, T])

    n = a.size
    full = est(np.arange(n))
    rng = np.random.default_rng(seed)
    boot = np.array([est(rng.integers(0, n, n)) for _ in range(n_boot)])
    lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)
    names = ("path_frac", "hw_frac", "usage_model_frac", "interaction_frac", "hw_total_frac", "usage_model_total_frac",
             "hw_frac_saltelli", "usage_model_frac_saltelli")
    out = {k: float(full[i]) for i, k in enumerate(names)}
    out.update({k + "_ci95": [float(lo[i]), float(hi[i])] for i, k in enumerate(names)})
    out.update(total_var_h2=float(full[8]), sd_total_h=math.sqrt(full[8]), E_var_path_h2=float(np.mean(ev)),
               mean_h=float(np.mean(np.concatenate([a, b]))),
               method="pick-freeze on F6-F7 inner means: Jansen (1999) totals, first order by the two-group "
                      "identity V_g = V - VT_other (Saltelli 2010 direct first order as a check); hw = §4 hardware "
                      "factors (UQ_SPEC), usage model = §4.6 factors (USAGE_UQ_SPEC)")
    return out


def l2_direct_split(m: np.ndarray, s: np.ndarray, n_paths: int, n_boot: int = 2000, seed: int = 6) -> dict:
    """Law of total variance from per-θ Monte Carlo means m and SDs s (n_paths each), with the Var(E) term
    corrected for the Monte Carlo noise of the per-θ means (subtract mean(s²)/n_paths); bootstrap over θ."""
    def est(ix):
        Ev = float(np.mean(s[ix] ** 2))
        Vm = float(np.var(m[ix], ddof=1) - np.mean(s[ix] ** 2) / n_paths)
        return Ev / (Ev + Vm), Ev, Vm
    f, Ev, Vm = est(np.arange(m.size))
    rng = np.random.default_rng(seed)
    bs = [est(rng.integers(0, m.size, m.size))[0] for _ in range(n_boot)]
    return dict(usage_frac=f, param_frac=1 - f, E_var_h2=Ev, var_E_h2=Vm, usage_frac_ci95=qs(bs, (2.5, 97.5)),
                n_theta=int(m.size), n_paths=n_paths)


def ctmc_check(res: dict, profile: str, inner_A: dict) -> dict:
    cfg = CTMC_CHECK
    chunks = [v for k, v in res.items() if k.startswith(f"ctmc:ctmc_{profile}_")]
    m = cfg["n_paths"]
    rows = []
    for r in chunks:
        idx = np.asarray(r["idx"])
        T2 = np.asarray(r["tte_h"], float).reshape(idx.size, m)
        T0 = np.asarray(r["l0_h"], float).reshape(idx.size, m)
        ev = np.asarray(r["event"]).reshape(idx.size, m)
        for j, i in enumerate(idx):
            rows.append(dict(i=int(i), m2=float(np.mean(T2[j])), s2=float(np.std(T2[j], ddof=1)),
                             m0=float(np.nanmean(T0[j])), s0=float(np.nanstd(T0[j], ddof=1)),
                             bo=float(np.mean(ev[j] == CODES.index("brownout"))),
                             nnan=int(np.sum(~np.isfinite(T2[j])) + np.sum(~np.isfinite(T0[j])))))
    rows.sort(key=lambda d: d["i"])
    i = np.array([d["i"] for d in rows])
    m2, s2, m0, s0 = (np.array([d[k] for d in rows]) for k in ("m2", "s2", "m0", "s0"))
    muD, sdC = inner_A["mu"][i], np.sqrt(inner_A["var"][i])
    r_pool = float(m2.sum() / m0.sum())
    e_mean = m2 / (r_pool * muD) - 1
    e_sd = s2 / (r_pool * sdC) - 1
    e_mean_L0 = m0 / muD - 1
    e_sd_L0 = s0 / sdC - 1
    se_mean = s2 / math.sqrt(m) / m2
    return dict(n_theta=int(i.size), n_paths_per_theta=m, h_s=cfg["h"], ratio_L2_L0_pooled=r_pool,
                ratio_L2_L0_q05_q50_q95=qs(m2 / m0),
                mean_err_pct=dict(max_abs=100 * float(np.max(np.abs(e_mean))),
                                  rms=100 * float(np.sqrt(np.mean(e_mean ** 2))),
                                  frac_within_3pct=float(np.mean(np.abs(e_mean) < 0.03))),
                sd_err_pct=dict(max_abs=100 * float(np.max(np.abs(e_sd))), rms=100 * float(np.sqrt(np.mean(e_sd ** 2))),
                                frac_within_10pct=float(np.mean(np.abs(e_sd) < 0.10)),
                                mean_signed=100 * float(np.mean(e_sd))),
                L0_vs_F6_mean_err_pct=dict(max_abs=100 * float(np.max(np.abs(e_mean_L0))),
                                           mean_signed=100 * float(np.mean(e_mean_L0))),
                L0_vs_F7_sd_err_pct=dict(max_abs=100 * float(np.max(np.abs(e_sd_L0))),
                                         mean_signed=100 * float(np.mean(e_sd_L0))),
                mc_se_mean_pct_median=100 * float(np.median(se_mean)),
                mc_se_sd_pct_approx=100 / math.sqrt(2 * (m - 1)),
                brownout_share_mean=float(np.mean([d["bo"] for d in rows])),
                n_censored=int(sum(d["nnan"] for d in rows)),
                pass_mean=bool(np.all(np.abs(e_mean) < 0.03)), pass_sd=bool(np.all(np.abs(e_sd) < 0.10)),
                split_L2_direct=l2_direct_split(m2, s2, m), split_L0_same_paths=l2_direct_split(m0, s0, m),
                split_F6F7_same_theta=dict(usage_frac=float(np.mean(sdC ** 2) / (np.mean(sdC ** 2) +
                                                                                 np.var(muD, ddof=1)))),
                per_theta=dict(idx=i, L2_mean_h=m2, L2_sd_h=s2, L0_mean_h=m0, L0_sd_h=s0, mu_dynkin_h=muD,
                               sd_clt_h=sdC))


def ratio_diagnostic(res: dict, chk: dict) -> dict:
    """Why the full-L2 check fails: the per-θ L2/L0 ratio of the CTMC day vs the constant-load S4 L2/L0 ratio
    of the same θ (25 °C, from 100 %), and the check errors when the pooled ratio is replaced by the S4 ratio
    rescaled to the same pooled mean (no fitted constants)."""
    pt = chk["per_theta"]
    i = np.asarray(pt["idx"])
    c = make_cell("S4", 25.0, 1.0, 1.0, 10.0, N_THETA)
    L2 = np.asarray(res["det:" + c["key"]]["tte_h"], float)[i]
    sc = build(c)
    pb = with_theta(sc.p, theta_first(0, i))
    P = sum(np.asarray(v, float) * np.ones(i.size) for v in components(sc.schedule, pb, 25.0).values())
    r_det = L2 / (np.asarray(cell.usable_energy(pb), float) / P)
    r_ct = pt["L2_mean_h"] / pt["L0_mean_h"]
    scaled = r_det * np.mean(r_ct) / np.mean(r_det)
    e = pt["L2_mean_h"] / (scaled * pt["mu_dynkin_h"]) - 1
    return dict(corr_ratio_ctmc_vs_S4=float(np.corrcoef(r_det, r_ct)[0, 1]), ratio_ctmc_min_max=qs(r_ct, (0, 100)),
                mean_err_with_S4_ratio_pct=dict(max_abs=100 * float(np.max(np.abs(e))),
                                                frac_within_3pct=float(np.mean(np.abs(e) < 0.03))))


# ------------------------------------------------------------------------------------------- step check, ABC
def step_check(res: dict) -> dict:
    out = {}
    for sid, Ta, s0, soh in STEP_CHECK_CELLS:
        c = make_cell(sid, Ta, s0, soh, 2.0, N_STEP_CHECK)
        r2 = res["step:" + c["key"]]
        r10 = res["det:" + c["key"]]
        T2 = np.asarray(r2["tte_h"], float)
        T10 = np.asarray(r10["tte_h"], float)[:N_STEP_CHECK]
        pos = T2 > 0
        rel = np.abs(T10[pos] / T2[pos] - 1)
        out[c["key"]] = dict(n=int(T2.size), max_rel_diff=float(rel.max()), median_rel_diff=float(np.median(rel)),
                             max_abs_diff_s=float(3600 * np.max(np.abs(T10 - T2))),
                             event_agreement=float(np.mean(np.asarray(r2["event"]) ==
                                                           np.asarray(r10["event"])[:N_STEP_CHECK])),
                             max_s_off_diff_pp=float(100 * np.nanmax(np.abs(np.asarray(r2["s_off"]) -
                                                                         np.asarray(r10["s_off"])[:N_STEP_CHECK]))),
                             runtime_h2_s=float(r2["runtime_s"]))
    return out


ABC_PARAMS = ("P_awake", "P_d0", "P_susp", "E_wake", "dU", "V_cut", "R0_ref", "P_SoC_max", "p_c")


def abc_bonn(res: dict, summaries: dict, inner: dict) -> tuple[dict, np.ndarray]:
    b = res["bonn:bonn_14pro_idle"]
    o = bonn_observations()[("14 Pro", "idle")]
    obs = float(o["obs_h"])
    tte = np.asarray(b["tte_h"], float)
    acc = np.abs(tte / obs - 1) <= 0.05
    th = theta_first(N_THETA)
    post = {k: dict(prior_q05_q50_q95=qs(th[k]), post_q05_q50_q95=qs(th[k][acc]),
                    width_ratio=(qs(th[k][acc])[2] - qs(th[k][acc])[0]) / (qs(th[k])[2] - qs(th[k])[0]))
            for k in ABC_PARAMS}
    narrowing = {}
    for key, s in summaries.items():
        r = res["det:" + key]
        if s["Tamb_C"] != 25.0 or s["s0_pct"] != 100 or s["n"] != N_THETA:
            continue
        T = np.asarray(r["tte_h"], float)
        a, b_ = qs(T), qs(T[acc])
        narrowing[key] = dict(prior_q05_q50_q95_h=a, post_q05_q50_q95_h=b_, width_ratio=(b_[2] - b_[0]) / (a[2] - a[0]),
                              median_shift_pct=100 * (b_[1] / a[1] - 1))
    ctmc_post = {}
    for prof, sid in (("typical", "S13"), ("heavy", "S14")):
        A = inner[prof]["A"]
        full, sub = split(A), split({k: v[acc] for k, v in A.items()})
        ctmc_post[sid] = dict(prior=dict(mixture_q05_q50_q95_h=full["mixture_q05_q50_q95_h"], param_frac=full["param_frac"],
                                         sd_total_h=full["sd_total_h"]),
                              post=dict(mixture_q05_q50_q95_h=sub["mixture_q05_q50_q95_h"], param_frac=sub["param_frac"],
                                        sd_total_h=sub["sd_total_h"]))
    wr = [v["width_ratio"] for v in narrowing.values()]
    out = dict(case="iPhone 14 Pro, black-screen idle (Bonn V1), Q = 3.32 Ah measured, R*Q scaling, h = 10 s",
               obs_h=obs, obs_device_min_max_h=[float(o["min_h"]), float(o["max_h"])], tolerance=0.05,
               nominal_model_h=float(b["nom_tte_h"]), prior_predictive_q05_q50_q95_h=qs(tte),
               obs_range_inside_prior_90=bool(qs(tte)[0] <= o["min_h"] and o["max_h"] <= qs(tte)[2]),
               n_prior=int(tte.size), n_accepted=int(acc.sum()), acceptance=float(acc.mean()),
               posterior=post, narrowing=narrowing, width_ratio_median=float(np.median(wr)),
               width_ratio_min_max=[float(min(wr)), float(max(wr))], ctmc=ctmc_post)
    return out, acc


# FINAL_SPEC §5.4 / §10 [PROTO] default-parameter TTEs (h) of cells run here: regression check of "nominal"
PROTO_NOMINAL = {
    "S1|s0=100|Tamb=25": 442.7, "S2a|s0=100|Tamb=25": 378.3, "S2b|s0=100|Tamb=25": 173.8, "S3|s0=100|Tamb=25": 38.88,
    "S4|s0=100|Tamb=25": 16.64, "S5a|s0=100|Tamb=25": 32.32, "S5b|s0=100|Tamb=25": 25.49, "S5c|s0=100|Tamb=25": 21.10,
    "S6|s0=100|Tamb=25": 11.83, "S7|s0=100|Tamb=25": 5.68, "S8a|s0=100|Tamb=25": 6.93, "S8b|s0=100|Tamb=25": 10.83,
    "S12|s0=100|Tamb=25": 34.92, "S12|s0=50|Tamb=25": 11.79, "S12|s0=20|Tamb=25": 8.70, "S16|s0=100|Tamb=25": 79.20,
    "S9|s0=100|Tamb=25": 26.99, "S9|s0=100|Tamb=-10|soh=0.8": 17.00, "S9|s0=100|Tamb=-20|soh=0.8": 12.33,
    "S10a|s0=100|Tamb=25|soh=0.8": 13.22, "S10a|s0=100|Tamb=0|soh=0.8": 12.80, "S10a|s0=100|Tamb=-10|soh=0.8": 12.15,
    "S4|s0=100|Tamb=0": 16.43, "S4|s0=100|Tamb=-10": 16.18, "S10b|s0=100|Tamb=25|soh=0.8": 29.19,
    "S10b|s0=100|Tamb=0|soh=0.8": 24.78, "S12|s0=100|Tamb=0": 34.50, "S11a|s0=100|Tamb=25": 12.73,
    "S11c|s0=100|Tamb=25": 5.33,
}


def proto_check(summaries: dict) -> dict:
    out = {k: dict(nominal_h=summaries[k]["nominal_h"], proto_h=v,
                   diff_pct=100 * (summaries[k]["nominal_h"] / v - 1)) for k, v in PROTO_NOMINAL.items()}
    return dict(cells=out, max_abs_diff_pct=max(abs(v["diff_pct"]) for v in out.values()))


# ------------------------------------------------------------------------------------------- figure, table
TABLE_ROWS = ("S1", "S2a", "S2b", "S3", "S4", "S5a", "S5b", "S5c", "S6", "S7", "S8a", "S8b", "S9", "S10a", "S10b",
              "S11a", "S11c", "S12", "S16")


def fmt_cell(s: dict | None) -> str:
    if s is None:
        return "--"
    m, a, b = s["median_h"], s["q05_h"], s["q95_h"]
    d = 0 if m >= 100 else 1
    return f"{m:.{d}f} [{a:.{d}f},\\,{b:.{d}f}]"


def write_table(summaries: dict) -> str:
    def get(sid, Ta, s0):
        soh = "|soh=0.8" if sid in ("S10a", "S10b") else ""
        return summaries.get(f"{sid}|s0={s0}|Tamb={Ta:g}{soh}")

    cols = [(25.0, 100), (25.0, 50), (25.0, 20), (0.0, 100), (-10.0, 100)]
    lines = [r"\setlength{\tabcolsep}{4pt}%", r"\begin{tabular}{l" + "r" * len(cols) + "}", r"\toprule",
             r" & \multicolumn{3}{c}{$T_{\mathrm{amb}} = 25\,^\circ$C} & $0\,^\circ$C & $-10\,^\circ$C \\",
             r"\cmidrule(lr){2-4}\cmidrule(lr){5-5}\cmidrule(lr){6-6}",
             r"Scenario & $s_0 = 100$\,\% & 50\,\% & 20\,\% & 100\,\% & 100\,\% \\", r"\midrule"]
    names = {"S1": "S1 standby, Wi-Fi", "S2a": "S2a standby, LTE", "S2b": "S2b standby, 5G push",
             "S8a": "S8a video call, LTE", "S8b": "S8b video call, Wi-Fi", "S9": "S9 cold hike$^{a}$",
             "S10a": "S10a social, SOH 0.8", "S10b": "S10b day, SOH 0.8", "S11a": "S11a social, $-110$\\,dBm",
             "S11c": "S11c call, $-110$\\,dBm", "S12": "S12 scheduled day", "S16": "S16 audio, screen off"}
    for sid in TABLE_ROWS:
        lab = names.get(sid, LABELS[sid]).replace("C)", "\\,$^\\circ$C)")
        row = [fmt_cell(get(sid, Ta, s0)) for Ta, s0 in cols]
        if sid in STANDBY:
            row[3] = row[3] + "$^{b}$" if row[3] != "--" else row[3]
            row[4] = row[4] + "$^{b}$" if row[4] != "--" else row[4]
        lines.append(f"{lab} & " + " & ".join(row) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    tex = "\n".join(lines) + "\n"
    tex = ("% T3 (FINAL_SPEC §9): L2 TTE median [5 %, 95 %] in hours over the joint §4 parameter sample\n"
           "% (scrambled Sobol', N = 2048, seed 2026; rk4_batch).  a: T0 = 20 C (carried out from indoors).\n"
           "% b: N = 512 (balanced Sobol' subset).  '--': not run (S15 starts only for S4, S5a, S7, S12).\n"
           "% A lower bound of 0.0 means that more than 5 % of the samples brown out at the first burst (TTE = 0).\n"
           "% \\tabcolsep is set to 4pt (local to the float that \\inputs this file).\n"
           "% Generated by analysis/run_uq.py; do not edit by hand.\n") + tex
    write_tex("tab_tte_uq.tex", tex)
    return tex


def make_figure(res: dict, summaries: dict, splits: dict) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from style import (BLUE, GRID, INK, INK2, MUTED, ORANGE, SEQ_BLUE, W1, apply, panel_label, save)

    apply()
    fig = plt.figure(figsize=(W1, 4.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], width_ratios=[1.45, 1.0], hspace=0.55, wspace=0.42)
    ax = fig.add_subplot(gs[0, :])
    # (a) all scenarios at 25 °C from 100 %, log TTE, ordered by nominal TTE
    keys = [k for k, s in summaries.items() if s["Tamb_C"] == 25.0 and s["s0_pct"] == 100 and s["scenario"] != "S9"
            or (s["scenario"] == "S9" and s["Tamb_C"] == 25.0 and s["soh"] == 1.0)]
    keys.sort(key=lambda k: -summaries[k]["nominal_h"])
    data = []
    for k in keys:  # TTE = 0 (brown-out at the first burst) cannot be shown on a log axis; counted in the JSON
        T = np.asarray(res["det:" + k]["tte_h"], float)
        data.append(np.log10(T[np.isfinite(T) & (T > 0)]))
    pos = np.arange(len(keys))
    vp = ax.violinplot(data, positions=pos, widths=0.8, showextrema=False, points=120)
    for b in vp["bodies"]:
        b.set_facecolor(BLUE)
        b.set_edgecolor(BLUE)
        b.set_alpha(0.35)
        b.set_linewidth(0.6)
    for i, k in enumerate(keys):
        s = summaries[k]
        ax.plot([i, i], np.log10([s["q05_h"], s["q95_h"]]), color=BLUE, lw=1.4, solid_capstyle="butt")
        ax.plot(i, np.log10(s["median_h"]), "o", ms=3.2, color=BLUE, mec="white", mew=0.5)
        ax.plot(i, np.log10(s["nominal_h"]), "_", ms=9, mew=1.3, color=INK)
    ticks = [5, 10, 20, 50, 100, 200, 500]
    ax.set_yticks(np.log10(ticks))
    ax.set_yticklabels([str(t) for t in ticks])
    ax.set_ylabel("TTE (h, log scale)")
    ax.set_xticks(pos)
    ax.set_xticklabels([summaries[k]["scenario"] for k in keys], fontsize=7.2)
    ax.set_xlim(-0.6, len(keys) - 0.4)
    ax.grid(axis="x", visible=False)
    ax.set_title("25 °C, from 100 % (S10a/b: SOH 0.8): TTE over the §4 parameter sample, N = 2048", loc="left",
                 fontsize=8.5, fontweight="normal", color=INK2)
    hand = [Patch(facecolor=BLUE, alpha=0.35, edgecolor=BLUE, label="density"),
            plt.Line2D([], [], color=BLUE, lw=1.4, label="5–95 %"),
            plt.Line2D([], [], color=BLUE, marker="o", ls="", ms=3.2, label="median"),
            plt.Line2D([], [], color=INK, marker="_", ls="", ms=9, mew=1.3, label="default θ")]
    ax.legend(handles=hand, loc="upper right", ncol=4, fontsize=7, handlelength=1.4, columnspacing=1.0)
    ax.set_ylim(np.log10(3.0), np.log10(1500))
    panel_label(ax, "(a)")

    # (b) cold: TTE relative to the 25 °C default-θ value, per scenario, colour = T_amb (sequential, darker = colder)
    ax = fig.add_subplot(gs[1, 0])
    sel = [("S4", 1.0), ("S6", 1.0), ("S7", 1.0), ("S12", 1.0), ("S9", 1.0), ("S10a", 0.8), ("S10b", 0.8)]
    tcol = {25.0: SEQ_BLUE[3], 0.0: SEQ_BLUE[7], -10.0: SEQ_BLUE[11]}
    for j, (sid, soh) in enumerate(sel):
        tag = "" if soh == 1.0 else f"|soh={soh:g}"
        ref = summaries[f"{sid}|s0=100|Tamb=25{tag}"]["nominal_h"]
        for t_i, Ta in enumerate(TEMPS):
            k = f"{sid}|s0=100|Tamb={Ta:g}{tag}"
            x = j + (t_i - 1) * 0.27
            T = np.asarray(res["det:" + k]["tte_h"], float) / ref
            v = ax.violinplot([T], positions=[x], widths=0.26, showextrema=False, points=100)
            for bd in v["bodies"]:
                bd.set_facecolor(tcol[Ta])
                bd.set_edgecolor(tcol[Ta])
                bd.set_alpha(0.75)
                bd.set_linewidth(0.5)
            s = summaries[k]
            ax.plot([x, x], [s["q05_h"] / ref, s["q95_h"] / ref], color=INK2, lw=0.7)
            ax.plot(x, s["median_h"] / ref, "o", ms=2.2, color=INK, mec="none")
    ax.axhline(1.0, color=MUTED, lw=0.7, ls=(0, (3, 2)))
    ax.set_xticks(range(len(sel)))
    ax.set_xticklabels([sid for sid, _ in sel], fontsize=7.2)
    ax.set_ylabel("TTE / default at 25 °C")
    ax.grid(axis="x", visible=False)
    ax.set_title("Cold ambient (S10a/b: SOH 0.8; S9: T0 = 20 °C)", loc="left", fontsize=8.5, fontweight="normal",
                 color=INK2)
    hand = [Patch(facecolor=tcol[Ta], edgecolor=tcol[Ta], alpha=0.75, label=f"{Ta:g} °C") for Ta in TEMPS]
    ax.legend(handles=hand, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3, fontsize=7,
              handlelength=1.0, columnspacing=1.0)
    panel_label(ax, "(b)")

    # (c) S13/S14 variance split (law of total variance)
    ax = fig.add_subplot(gs[1, 1])
    from style import AQUA
    comps = [("usage path", "path", BLUE), ("hardware θ", "hw", ORANGE), ("usage-model θ", "um", AQUA),
             ("interaction", "ia", GRID)]
    rows = []
    for sid, lab in (("S13", "S13"), ("S14", "S14")):
        sp = splits[sid]
        l2 = sp["l2_check"]["split_L2_direct"]
        dc = sp["decomposition"]
        rows.append((f"{lab} F6–F7", dict(path=sp["usage_frac"], hw=sp["param_frac"])))
        rows.append((f"{lab} L2, 50 θ", dict(path=l2["usage_frac"], hw=l2["param_frac"])))
        rows.append((f"{lab} + usage model", dict(path=dc["path_frac"], hw=dc["hw_frac"], um=dc["usage_model_frac"],
                                                ia=max(dc["interaction_frac"], 0.0))))
    y = np.array([7.2, 6.2, 5.2, 3.0, 2.0, 1.0])[:len(rows)]
    for (lab, d), yy in zip(rows, y):
        left = 0.0
        for name, key, col in comps:
            w = d.get(key, 0.0)
            if w <= 0:
                continue
            ax.barh(yy, w, left=left, color=col, height=0.72, edgecolor="white", lw=0.6)
            if w > 0.1:
                ax.text(left + w / 2, yy, f"{100 * w:.0f}", ha="center", va="center", fontsize=6.6,
                        color="white" if col in (BLUE, ORANGE) else INK)
            left += w
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=6.8)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(["0", "25", "50", "75", "100 %"])
    ax.set_xlabel("share of Var(TTE)")
    ax.grid(axis="y", visible=False)
    ax.set_title("Stochastic days (S13 typical, S14 heavy)", loc="left", fontsize=8.5, fontweight="normal",
                 color=INK2)
    hand = [Patch(facecolor=col, label=name) for name, _, col in comps]
    ax.legend(handles=hand, loc="upper center", bbox_to_anchor=(0.32, -0.3), ncol=2, fontsize=6.6,
              handlelength=0.9, columnspacing=0.8)
    panel_label(ax, "(c)")
    save(fig, "fig_uq_violins")


# ------------------------------------------------------------------------------------------- main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=int(os.environ.get("UQ_WORKERS", "3")))
    ap.add_argument("--fresh", action="store_true", help="delete the per-cell cache first")
    ap.add_argument("--only", default=None, help="regex on 'kind:key': run only matching jobs (see compute)")
    ap.add_argument("--no-aggregate", action="store_true", help="compute/cache only; skip JSON, figure, table")
    a = ap.parse_args(argv)
    T0 = time.time()
    if a.fresh and CACHE.exists():
        shutil.rmtree(CACHE)
    res = compute(a.workers, a.only)
    t_compute = time.time() - T0
    missing = [c["kind"] + ":" + c["key"] for c in all_jobs() if c["kind"] + ":" + c["key"] not in res]
    if a.no_aggregate or missing:
        print(f"[uq] computed/cached; {len(missing)} jobs still missing -> no aggregation"
              + (f" (e.g. {missing[0]})" if missing else ""))
        return 0

    # 1. deterministic cells
    summaries = {}
    for c in cell_list():
        summaries[c["key"]] = cell_summary(c, res["det:" + c["key"]])

    # 2. delta method (every s0 = 100 % cell)
    dm = {}
    for c in cell_list():
        if c["s0"] == 1.0:
            dm[c["key"]] = delta_method(c, summaries[c["key"]])
    dm25 = {k: v for k, v in dm.items() if "|Tamb=25" in k}
    dm_digest = dict(rule="CV_TTE^2 ~ CV_Eu^2 + sum_i s_i^2 CV_Pi^2 (L0, quasi-steady powers; FINAL_SPEC §9 item 2)",
                     ratio_L2_to_delta_25C_q05_q50_q95=qs([v["ratio_L2_to_delta"] for v in dm25.values()]),
                     ratio_L0sample_to_delta_25C_q05_q50_q95=qs([v["cv_L0_sample"] / v["cv_delta"]
                                                                 for v in dm25.values()]),
                     cells=dm)

    # 3. S13/S14 variance split
    th = theta_first(N_THETA)
    us = usage_all()
    # independent second sample B for the pick-freeze decomposition (pseudo-random: see usage_all)
    th_b = sample(UQ_SPEC, N_THETA, method="random", seed=SEED_B)
    us_b = sample(USAGE_UQ_SPEC, N_THETA, method="random", seed=SEED_USAGE_B)
    inner, splits = {}, {}
    for prof, sid in (("typical", "S13"), ("heavy", "S14")):
        A = ctmc_inner(prof, th, None, N_THETA)
        B = ctmc_inner(prof, None, us, N_THETA)
        C = ctmc_inner(prof, th, us, N_THETA)
        inner[prof] = dict(A=A, B=B, C=C)
        sA, sB, sC = split(A), split(B), split(C)
        nom = ctmc_moments(CTMC(prof).generator(), CTMC(prof).state_powers(RP63, 25.0), float(cell.usable_energy(RP63)))
        chk = ctmc_check(res, prof, A)
        fB = ctmc_inner(prof, th_b, us_b, N_THETA)
        f_hwB = ctmc_inner(prof, th_b, us, N_THETA)
        f_usB = ctmc_inner(prof, th, us_b, N_THETA)
        C3 = group_decomposition(C, fB, f_hwB, f_usB)
        r = chk["ratio_L2_L0_pooled"]
        splits[sid] = dict(
            usage_frac=sA["usage_frac"], param_frac=sA["param_frac"], design="A: hardware theta only (FINAL_SPEC §9)",
            nominal=dict(mu_dynkin_h=nom["mu_dynkin"], sd_clt_h=nom["sd_clt"]), A_hardware=sA, B_usage_model=sB,
            C_both=sC, decomposition=C3,
            L2_equivalent=dict(ratio_L2_L0=r, mean_h=r * sA["mean_h"], mixture_q05_q50_q95_h=[r * v for v in
                                                                                             sA["mixture_q05_q50_q95_h"]],
                               note="F6-F7 L0 mixture scaled by the pooled L2/L0 mean ratio of the full-L2 check"),
            l2_check={k: v for k, v in chk.items() if k != "per_theta"}, l2_check_per_theta=chk["per_theta"],
            l2_check_ratio_diagnostic=ratio_diagnostic(res, chk))

    # 4. step check and ABC
    sc_chk = step_check(res)
    abc, acc = abc_bonn(res, summaries, inner)

    # 5. outputs
    make_figure(res, summaries, splits)
    write_table(summaries)
    dev = [
        "constant-load cells use h = 10 s instead of the §9 default 2 s (no input switches); step_check compares "
        "h = 10 s with h = 2 s on the first 256 Sobol' points of five cells",
        "standby S1/S2a/S2b at 0 and -10 C use N = 512 (first 512 Sobol' points, a balanced subset), not 2048",
        "S15 starts (50 / 20 %) only for the §5.2 S15 set S4, S5a, S7, S12 (S13 is stochastic); S11 only at 25 C",
        "S9 is run at the T3 temperatures 25 / 0 / -10 C (T0 = 20 C), plus -10 and -20 C at SOH 0.8",
        "variance split inner moments are L0 (F6-F7); the L2 check rescales by the pooled L2/L0 ratio measured on "
        "the same paths",
        "usage-model factors (§4.6) and the second pick-freeze sample B are pseudo-random (seeds 2027-2029), "
        "paired by index with the Sobol' θ sample: two index-paired scrambled Sobol' sequences are not independent",
        "S13/S14 hardware-vs-usage-model split (an addition to §9) uses a pick-freeze group decomposition (Jansen "
        "totals, first order by the two-group identity) on the F6-F7 inner means; the §9 two-way split (usage path "
        "vs hardware) is reported as usage_frac/param_frac, and a direct full-L2 split is estimated from the "
        "50 x 200 check",
    ]
    out = dict(
        meta=meta("run_uq.py", seed=dict(theta=SEED, usage_factors=SEED_USAGE, pick_freeze_B=[SEED_B, SEED_USAGE_B],
                                     ctmc_paths=f"{SEED}+100*profile+chunk"),
                  solver="rk4_batch via model.run_batch (event-synced; LSODA for the default-theta 'nominal')",
                  h_s=dict(constant_load=10.0, standby=10.0, S9_S12_S10b=2.0, ctmc=10.0, bonn=10.0),
                  n_theta=N_THETA, sampling="scipy.stats.qmc.Sobol(scramble=True) over all UQ_SPEC factors",
                  runtime_compute_s=t_compute, workers=a.workers,
                  cpu_jobs_s=float(sum(float(r["runtime_s"]) for r in res.values()))),
        deviations=dev, **summaries,
        variance_split=splits, delta_method_check=dm_digest, abc_bonn=abc, step_check=sc_chk,
        proto_nominal_check=proto_check(summaries),
        runtime_s=time.time() - T0)
    write_json("uq.json", out)
    print(f"[uq] done in {time.time() - T0:.0f} s (compute {t_compute:.0f} s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
