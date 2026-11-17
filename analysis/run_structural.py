"""AN10: structural sensitivity (model-form uncertainty), FINAL_SPEC §10.4: variants X1–X18 against X0.

Every variant is implemented in analysis code; ``src/phonebatt`` is not modified.

* X2–X14, X18: the vectorised union model of ``analysis/structural_variants.py`` (a superset of the core equations
  whose extra terms vanish for X0), integrated with the core ``rk4_batch`` and the event logic of
  ``phonebatt.model.run_batch``; h = 2 s (h = 10 s for the standby-like scenarios S1, S2a, S2b, S16 and the CTMC
  days S13, S14, §6.3).  X0 is computed by the same union model at the same h, so ΔTTE compares like with like;
  X0 itself is checked against the core LSODA reference (and bit-for-bit against ``run_batch``).
* X1 (L0 energy bucket): exact; usable energy E_u(z0) depleted at the F5 quasi-steady state powers (the L0 used by
  the paper's analytic layer, §5.3/§7); ends at displayed 0 %.
* X15 (lognormal sessions, CV 1.5) and X16 (diurnal CTMC): CTMC days only, on usage paths from
  ``stochastic_common.UsageSim`` (X0 paths use the exponential law with the same seed: common random numbers).
* X17 (calendar law t^0.75): aging layer only.  A1 is re-run with q_cal = k'·(∫ stress^{4/3} dt)^{3/4}
  (k' so that q_cal(365 d, 25 °C, g = 1) = 2.7 %, k_H/k_L recalibrated by protocol C0); reported as ΔSOH and
  as the TTE change at the SOH reached after 2 years (core LSODA ``run``).

Variant definitions (choices beyond the §10.4 one-liners are stated here and in results/structural.json meta):
X2 L1 = series R0 + R1, no RC state, isothermal at T_amb (so T_s = T_amb: no throttling), no burst envelope;
X3 = L2 with T ≡ T_amb and no throttling; X4 I = P/3.87 V, V = Ũ − R0 I, burst I_pk = (P + ΔP_b)/3.87 V;
X5 = 2RC 2.5 mΩ/6.7 s + 14.1 mΩ/47 s scaled by R·Q (×2.28/4.70), Arrhenius as R1/τ1; X5n = same split
normalised to R_a + R_b = R1,ref (isolates the RC order); X6 R0 = R0∞[1 + A_R e^(−z/z0)](1 + r)·[f + (1 − f)a_R],
f = 35/60; X7 no burst envelope; X8 Δ_sys = 0; X9 KiBaM (c 0.166, κ_K 40 h⁻¹, E_a,K 35 kJ/mol) replacing the RC
state with R1 kept as a series resistance (as AN6); X10 cellular tail occupancy dr/dt = λ(1 − r) − r/T_tail with
tail power m(P_tail − P_idle) r and promotions λ(1 − r) m E_pro (Wi-Fi tail stays quasi-steady; during
real-time traffic r is driven to 1 and m·P_rt is the whole radio power); X11a/b
γ_b = 1 / 2.9; X12 s0 and s_off in D2; X13 P_leak = 0.1 W [exp(k_lk(T_s − 25 °C)) − 1] while awake, k_lk = 0.04
K⁻¹ (upper end; 0 is X0); X14 AU/AU-L with a = 0; X15 lognormal sessions CV 1.5; X16 λ_on = 0 at 00–07 h,
5.5 h⁻¹ otherwise (start 07:00); X17 t^0.75 calendar law; X18a/b OCV stretch (0.192 V, 5.44) / (0.272 V, 2.90).

Writes results/structural.json (key "{X}|{S}" -> {dtte_pct, ds_off_pct, ...}), paper/tables/tab_structural.tex
(T4), paper/figures/fig_structural.pdf (+ PNG) and results/global_summary.md (digest of AN9 + AN10).

Run:  cd /home/user/MCM && PYTHONPATH=src:analysis python3 analysis/run_structural.py [--paths 200]
"""
from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import replace

import numpy as np

from phonebatt import RP63, aged, cell, run, run_batch
from phonebatt.aging import AgingState, day_profile, run_years
from phonebatt.cell import a_deg, a_R
from phonebatt.power import steady_power
from phonebatt.usage import CTMC, PRIMITIVES, Paths, Schedule, s9_schedule, s12_schedule
from basics_common import RESULTS, meta, write_json, write_tex
from stochastic_common import UsageSim, l0_from_paths
from structural_variants import X9_DEFAULT, Column, UnionBatch, Variant, schedule_path

SEED = 2026
H_A, H_B, H_C = 2.0, 10.0, 10.0
N_PATHS = 400
RQ = 2.28 / 4.70  # R·Q scaling from the Ai2020 cell (2.28 Ah) to RP-63 (4.70 Ah)


# ------------------------------------------------------------------------------------------- variants
def _audio_dsp(u):
    return replace(u, awake=0.0) if u.name in ("AU", "AU-L") else u


VARIANTS = {
    "X0": Variant("X0", "baseline L2/L3"),
    "X2": Variant("X2", "L1 quasi-static ECM, isothermal", series=True, iso=True, throttle=False, burst=False),
    "X3": Variant("X3", "L2 isothermal, no throttling", iso=True, throttle=False),
    "X4": Variant("X4", "constant current I = P/3.87 V", cc=True),
    "X5": Variant("X5", "2RC (R·Q-scaled fit)", rc2=(2.5e-3 * RQ, 6.7, 14.1e-3 * RQ, 47.0)),
    "X5n": Variant("X5n", "2RC, R_a + R_b = R1,ref", rc2=(0.010 * 2.5 / 16.6, 6.7, 0.010 * 14.1 / 16.6, 47.0)),
    "X6": Variant("X6", "R0 split: 35 mΩ PCM + 25 mΩ cell", pcm_frac=35.0 / 60.0),
    "X7": Variant("X7", "no burst envelope", burst=False),
    "X8": Variant("X8", "single threshold V_sys = V_cut", params=dict(d_sys=0.0)),
    "X9": Variant("X9", "KiBaM instead of RC", kibam=dict(X9_DEFAULT)),
    "X10": Variant("X10", "RRC tail as ODE state", tail_ode=True),
    "X11a": Variant("X11a", "brightness law γ_b = 1", params=dict(gamma_b=1.0)),
    "X11b": Variant("X11b", "brightness law γ_b = 2.9", params=dict(gamma_b=2.9)),
    "X12": Variant("X12", "SOC output D2 (s0, s_off)", soc_mode="D2"),
    "X13": Variant("X13", "SoC leakage, k_lk = 0.04 /K", k_lk=0.04),
    "X14": Variant("X14", "audio on DSP (a = 0 in AU)", input_map=_audio_dsp),
    "X18a": Variant("X18a", "OCV stretch 4.40 V / 3.85 V", params=dict(dV_hv=0.192, p_hv=5.44)),
    "X18b": Variant("X18b", "OCV stretch 4.48 V / 3.89 V", params=dict(dV_hv=0.272, p_hv=2.90)),
}
LABELS = {"X1": "L0 energy bucket", "X15": "lognormal sessions, CV 1.5", "X16": "diurnal CTMC",
          "X17": "calendar law t^0.75", **{k: v.label for k, v in VARIANTS.items()}}
ORDER = ["X1", "X2", "X3", "X4", "X5", "X5n", "X6", "X7", "X8", "X9", "X10", "X11a", "X11b", "X12", "X13", "X14",
         "X15", "X16", "X17", "X18a", "X18b"]

# ------------------------------------------------------------------------------------------- scenarios
SO, GM = Schedule.constant("SO"), Schedule.constant("GM")
DET = {  # key: (label, schedule, Tamb_C, soh, T0_C, s0, group)
    "S1": ("standby, home Wi-Fi", Schedule.constant("SB-W"), 25.0, 1.0, None, 1.0, "B"),
    "S2a": ("standby, cellular", Schedule.constant("SB-L"), 25.0, 1.0, None, 1.0, "B"),
    "S2b": ("standby, 5G NSA, heavy push", Schedule.constant("SB-5"), 25.0, 1.0, None, 1.0, "B"),
    "S3": ("reading", Schedule.constant("RD"), 25.0, 1.0, None, 1.0, "A"),
    "S4": ("social, Wi-Fi", SO, 25.0, 1.0, None, 1.0, "A"),
    "S4@s0=50": ("social, from 50 %", SO, 25.0, 1.0, None, 0.5, "A"),
    "S5a": ("video, Wi-Fi", Schedule.constant("VW"), 25.0, 1.0, None, 1.0, "A"),
    "S5b": ("video, LTE", Schedule.constant("VL"), 25.0, 1.0, None, 1.0, "A"),
    "S6": ("navigation", Schedule.constant("NV"), 25.0, 1.0, None, 1.0, "A"),
    "S7@15": ("gaming, 15 °C", GM, 15.0, 1.0, None, 1.0, "A"),
    "S7": ("gaming, 25 °C", GM, 25.0, 1.0, None, 1.0, "A"),
    "S7@35": ("gaming, 35 °C", GM, 35.0, 1.0, None, 1.0, "A"),
    "S8a": ("video call, LTE", Schedule.constant("VC"), 25.0, 1.0, None, 1.0, "A"),
    "CM": ("camera", Schedule.constant("CM"), 25.0, 1.0, None, 1.0, "A"),
    "S9@-10/0.8": ("cold hike, −10 °C, SOH 0.8", s9_schedule(-10.0), -10.0, 0.8, 20.0, 1.0, "A"),
    "S9@-20/0.8": ("cold hike, −20 °C, SOH 0.8", s9_schedule(-20.0), -20.0, 0.8, 20.0, 1.0, "A"),
    "S10a@-10/0.8": ("social, −10 °C, SOH 0.8", SO, -10.0, 0.8, None, 1.0, "A"),
    "S10a@-20/0.8": ("social, −20 °C, SOH 0.8", SO, -20.0, 0.8, None, 1.0, "A"),
    "S10b@0/0.8": ("scheduled day, 0 °C, SOH 0.8", s12_schedule(), 0.0, 0.8, None, 1.0, "A"),
    "S11a@-110": ("social, LTE, −110 dBm", Schedule.constant(replace(PRIMITIVES["SO-L"], RSRP=-110.0)), 25.0, 1.0,
                  None, 1.0, "A"),
    "S12": ("scheduled day", s12_schedule(), 25.0, 1.0, None, 1.0, "A"),
    "S12@s0=20": ("scheduled day, from 20 %", s12_schedule(), 25.0, 1.0, None, 0.2, "A"),
    "S16": ("audio, screen off", Schedule.constant("AU"), 25.0, 1.0, None, 1.0, "B"),
}
CTMC_SC = {"S13": ("typical day (CTMC)", "typical"), "S14": ("heavy day (CTMC)", "heavy")}
SC_ORDER = ["S1", "S2a", "S2b", "S3", "S4", "S4@s0=50", "S5a", "S5b", "S6", "S7@15", "S7", "S7@35", "S8a", "CM",
            "S11a@-110", "S16", "S12", "S12@s0=20", "S13", "S14", "S9@-10/0.8", "S9@-20/0.8", "S10a@-10/0.8",
            "S10a@-20/0.8", "S10b@0/0.8"]


def _p(soh):
    return aged(RP63, 1.0 - soh) if soh < 1.0 else RP63


# ------------------------------------------------------------------------------------------- X16 diurnal usage
class DiurnalSim(UsageSim):
    """UsageSim with a diurnal session-start rate: λ_on = 0 from 00:00 to 07:00, ``lam_day`` otherwise.

    t = 0 is 07:00.  Idle sojourns are drawn by inverting the cumulative rate (the same standard-exponential draw
    as UsageSim, so the random streams line up with the homogeneous X0 paths until the first night).
    """

    def __init__(self, ctmc, lam_day=5.5, off=(0.0, 7.0), **kw):
        super().__init__(ctmc, **kw)
        self.lam_day = float(lam_day)
        self.off = off
        self.active = 24.0 - (off[1] - off[0])

    def _A(self, c):  # active hours from midnight of day 0 to clock time c (h)
        return self.active * np.floor(c / 24.0) + np.clip(np.mod(c, 24.0) - self.off[1], 0.0, self.active)

    def _Ainv(self, a):
        return 24.0 * np.floor(a / self.active) + self.off[1] + np.mod(a, self.active)

    def idle_at(self, rng, t):
        e = rng.exponential(1.0, t.size)
        c0 = t + 7.0
        return self._Ainv(self._A(c0) + e / self.lam_day) - c0

    def run(self, n, seed, P=None, Eu=None, store=True, t_max_h=400.0, x0=0):
        """Path sampler (store=True only); a copy of UsageSim.run with the idle draws made time dependent."""
        if not store:
            raise NotImplementedError
        rng = np.random.default_rng(seed)
        self._rng_j = np.random.default_rng([seed, 7919])
        a = np.full(n, x0 // 2, dtype=np.int64)
        c = np.full(n, x0 % 2, dtype=np.int64)
        k = self._kidx(rng, n)
        t = np.zeros(n)
        ra = np.where(a == 0, self.idle_at(rng, t), 0.0)
        rc = self._conn(rng, c)
        idx = np.arange(n)
        times, states = [np.zeros(n)], [self._index(a, c, k)]
        while idx.size:
            dt = np.minimum(ra, rc)
            t = t + dt
            act = ra <= rc
            ia = np.nonzero(act)[0]
            if ia.size:
                was_idle = a[ia] == 0
                st = ia[was_idle]
                if st.size:
                    j = np.searchsorted(self.cum_p, rng.random(st.size), side="right")
                    j = np.minimum(j, self.p.size - 1)
                    a[st] = j + 1
                    ra[st] = self._session(rng, j)
                en = ia[~was_idle]
                if en.size:
                    a[en] = 0
                    ra[en] = self.idle_at(rng, t[en])
                k[ia] = self._kidx(rng, ia.size)
                rc[ia] -= dt[ia]
            ic = np.nonzero(~act)[0]
            if ic.size:
                c[ic] = 1 - c[ic]
                rc[ic] = self._conn(rng, c[ic])
                ra[ic] -= dt[ic]
            tt = np.full(n, np.inf)
            tt[idx] = t
            ss = np.zeros(n, dtype=np.int64)
            ss[idx] = self._index(a, c, k)
            times.append(tt * 3600.0)
            states.append(ss)
            keep = t < t_max_h
            if not keep.all():
                idx, a, c, k, ra, rc, t = (v[keep] for v in (idx, a, c, k, ra, rc, t))
        return Paths(np.stack(times, axis=1), np.stack(states, axis=1), 10.0)


# ------------------------------------------------------------------------------------------- X1 (L0)
def l0_schedule(sch: Schedule, Tamb_C, p, s0=1.0):
    """Exact L0 TTE [h]: E_u(z0) depleted at the F5 quasi-steady segment powers (segment ambients honoured)."""
    Eu = float(cell.usable_energy(p, cell.z_from_s0(s0, p)))
    P = np.array([float(steady_power(u, p, Tamb_C if u.T_amb_C is None else u.T_amb_C)) for _, u in sch])
    d_h = sch.durations / 3600.0
    Eper = float(np.sum(P * d_h))
    n_full = math.floor(Eu / Eper)
    E = Eu - n_full * Eper
    t = n_full * sch.period / 3600.0
    for Pk, dk in zip(P, d_h):
        if Pk * dk >= E:
            return t + E / Pk
        E -= Pk * dk
        t += dk
    return t


# ------------------------------------------------------------------------------------------- X17 (aging law)
K75 = 0.027 / 365.0 ** 0.75  # d^-0.75: q_cal(365 d, 25 °C, g = 1) = 2.7 %, as k_c for the √t law


def _advance_x17(st, prof, p, kH, kL):
    """E17–E19 over one day with the calendar law q = k'(∫ s^{4/3} dt)^{3/4} (st.Q_cal holds W = q^{4/3})."""
    from phonebatt.aging import g_stress
    W, QH, QL = st
    for s in prof:
        W += float(K75 * a_deg(s.T_K, p.E_cal) * g_stress(s.z, p)) ** (4.0 / 3.0) * s.dt_h / 24.0
        QH += float(kH * a_deg(s.T_K, p.E_H)) ** 2 * s.dn
        if s.dn_ch > 0:
            QL += float(kL * a_R(s.T_K, p.E_L) * math.exp(p.beta_L * (s.c_ch - p.c_ref))) ** 2 * s.dn_ch
    return AgingState(W, QH, QL)


def _qloss_x17(st):
    return st.Q_cal ** 0.75 + math.sqrt(st.Q_H) + math.sqrt(st.Q_L)


def calibrate_x17(p=RP63):
    N = int(round(float(p.N_80)))
    prof = day_profile("full", 25.0, 1.0, p, cal_mode=True)
    st = AgingState()
    for _ in range(N):
        st = _advance_x17(st, prof, p, 1.0, 1.0)
    qc = st.Q_cal ** 0.75
    need = 0.2 - qc
    return (1.0 - p.w_L) * need / math.sqrt(st.Q_H), p.w_L * need / math.sqrt(st.Q_L), qc


def run_years_x17(years=3.0, p=RP63, P_day=0.440):
    kH, kL, qc1000 = calibrate_x17(p)
    En = float(cell.nominal_energy(p))
    st = AgingState()
    soh = []
    cache = {}
    for d in range(int(round(years * 365))):
        q = _qloss_x17(st)
        ds = round(min(0.95, 24.0 * P_day / (En * (1.0 - q))), 12)
        if ds not in cache:
            cache = {ds: day_profile("full", 25.0, ds, p, P_day=P_day)}
        st = _advance_x17(st, cache[ds], p, kH, kL)
        soh.append(1.0 - _qloss_x17(st))
    soh = np.array(soh)
    below = np.nonzero(soh <= 0.8)[0]
    return dict(soh=soh, years_to_80=(below[0] + 1) / 365.0 if below.size else math.inf, k_H=kH, k_L=kL,
                q_cal_1000=qc1000)


def x17_block() -> dict:
    t0 = time.time()
    base = run_years("full", 25.0, 3.0)
    alt = run_years_x17(3.0)
    out = dict(soh_X0=[float(base["soh"][365 * y - 1]) for y in (1, 2, 3)],
               soh_X17=[float(alt["soh"][365 * y - 1]) for y in (1, 2, 3)],
               years_to_80_X0=float(base["years_to_80"]), years_to_80_X17=float(alt["years_to_80"]),
               k_H_X17=alt["k_H"], k_L_X17=alt["k_L"], q_cal_1000_X17=alt["q_cal_1000"], k_cal_X17=K75)
    out["dsoh_2y_pp"] = 100 * (out["soh_X17"][1] - out["soh_X0"][1])
    tte = {}
    for key, (sch, Ta) in {"S4@2y": (SO, 25.0), "S12@2y": (s12_schedule(), 25.0),
                           "S10a@-10@2y": (SO, -10.0)}.items():
        r0 = run(sch, Ta, p=aged(RP63, 1.0 - out["soh_X0"][1]))
        r1 = run(sch, Ta, p=aged(RP63, 1.0 - out["soh_X17"][1]))
        tte[key] = dict(tte_X0_h=r0.tte_h, tte_h=r1.tte_h, dtte_pct=100 * (r1.tte_h / r0.tte_h - 1),
                        s_off_X0_pct=100 * r0.s_off, s_off_pct=100 * r1.s_off,
                        ds_off_pct=100 * (r1.s_off - r0.s_off), event=r1.event, event_X0=r0.event)
    out["tte_at_2y"] = tte
    out["runtime_s"] = round(time.time() - t0, 1)
    return out


# ------------------------------------------------------------------------------------------- runs
def det_columns(keys, variants, t_max_h):
    cols = []
    for k in keys:
        lab, sch, Ta, soh, T0, s0, _ = DET[k]
        t, s = schedule_path(sch, t_max_h[k] * 3600.0)
        for v in variants:
            cols.append(Column(v, k, [u for _, u in sch], t, s, Ta, soh, T0, s0, RP63))
    return cols


def ctmc_columns(key, prof, variants, paths, extra_paths=None):
    ct = CTMC(prof)
    us = ct.primitives()
    cols = []
    for v in variants:
        P = paths if (extra_paths is None or v.id not in extra_paths) else extra_paths[v.id]
        for i in range(P.n):
            tt = P.times[i]
            m = int(np.sum(np.isfinite(tt)))
            cols.append(Column(v, key, us, tt[:m].copy(), P.states[i, :m].copy(), 25.0, 1.0, None, 1.0, RP63,
                               table_key=key))
    return cols


def summarise_det(res, cols, out_rows):
    for j, c in enumerate(cols):
        out_rows.setdefault(c.key, {})[c.variant.id] = dict(tte_h=float(res["tte_h"][j]), event=str(res["event"][j]),
                                                           s_off_pct=100 * float(res["s_off"][j]),
                                                           s_off_D1_pct=100 * float(res["s_off_D1"][j]))


def compute(n_paths=N_PATHS) -> dict:
    t_all = time.time()
    timing = {}
    vlist = list(VARIANTS.values())

    # ---- X0 reference with the core LSODA solver (also sets the horizons)
    t0 = time.time()
    lsoda = {}
    for k, (lab, sch, Ta, soh, T0, s0, _) in DET.items():
        r = run(sch, Ta, s0=s0, T0_C=T0, p=_p(soh))
        lsoda[k] = dict(tte_h=r.tte_h, event=r.event, s_off_pct=100 * r.s_off)
    timing["lsoda_X0"] = round(time.time() - t0, 1)
    t_max = {k: max(3.0 * v["tte_h"], v["tte_h"] + 24.0) for k, v in lsoda.items()}

    # ---- deterministic groups
    rows: dict = {}
    for grp, h in (("A", H_A), ("B", H_B)):
        keys = [k for k in DET if DET[k][6] == grp]
        cols = det_columns(keys, vlist, t_max)
        t0 = time.time()
        ub = UnionBatch(cols, h)
        res = ub.run(max(t_max[k] for k in keys))
        timing[f"group_{grp}"] = round(time.time() - t0, 1)
        print(f"group {grp}: {len(cols)} columns, {timing[f'group_{grp}']} s", flush=True)
        summarise_det(res, cols, rows)

    # ---- h convergence of the variant effects (h = 1 s vs 2 s on three group-A scenarios)
    conv_keys = [k for k in ("S4", "S9@-20/0.8", "S10a@-20/0.8") if k in DET]
    t0 = time.time()
    cols = det_columns(conv_keys, vlist, t_max)
    res1 = UnionBatch(cols, 1.0).run(max(t_max[k] for k in conv_keys))
    rows_h1: dict = {}
    summarise_det(res1, cols, rows_h1)
    timing["h_convergence"] = round(time.time() - t0, 1)

    # ---- core run_batch consistency (X0 union model == core batch path, same h)
    t0 = time.time()
    core_check = {}
    for k in [k for k in ("S4", "S9@-10/0.8") if k in DET]:
        lab, sch, Ta, soh, T0, s0, _ = DET[k]
        b = run_batch(sch, Ta, None, s0=s0, h=H_A, T0_C=T0, p=_p(soh), n=1)
        core_check[k] = dict(run_batch_tte_h=float(b["tte_h"][0]), union_tte_h=rows[k]["X0"]["tte_h"],
                             abs_diff_h=abs(float(b["tte_h"][0]) - rows[k]["X0"]["tte_h"]),
                             event_equal=str(b["event"][0]) == rows[k]["X0"]["event"])
    timing["core_check"] = round(time.time() - t0, 1)

    # ---- CTMC days (common random numbers across variants)
    ctmc_rows = {}
    for key, (lab, prof) in CTMC_SC.items():
        t0 = time.time()
        ct = CTMC(prof)
        sim = UsageSim(ct, law="exp")
        P0 = sim.run(n_paths, SEED, store=True, t_max_h=160.0)
        P15 = UsageSim(ct, law="lognormal", cv=1.5).run(n_paths, SEED, store=True, t_max_h=160.0)
        P16 = DiurnalSim(ct, law="exp").run(n_paths, SEED, store=True, t_max_h=160.0)
        vX15 = Variant("X15", LABELS["X15"])
        vX16 = Variant("X16", LABELS["X16"])
        cols = ctmc_columns(key, prof, vlist + [vX15, vX16], P0, extra_paths={"X15": P15, "X16": P16})
        ub = UnionBatch(cols, H_C)
        res = ub.run(160.0)
        timing[f"ctmc_{key}"] = round(time.time() - t0, 1)
        print(f"{key}: {len(cols)} columns, {timing[f'ctmc_{key}']} s", flush=True)
        vids = [c.variant.id for c in cols]
        per = {}
        for vid in dict.fromkeys(vids):
            m = np.array([x == vid for x in vids])
            per[vid] = dict(tte=res["tte_h"][m], s_off=res["s_off"][m], s_off_D1=res["s_off_D1"][m],
                            event=res["event"][m])
        # X1: exact L0 on the X0 paths with F5 state powers
        Pss = ct.state_powers(RP63, 25.0, capped=True)
        per["X1"] = dict(tte=l0_from_paths(P0, Pss, float(cell.usable_energy(RP63))), s_off=np.zeros(n_paths),
                         s_off_D1=np.zeros(n_paths), event=np.full(n_paths, "empty (L0)"))
        ctmc_rows[key] = per

    # ---- X17
    x17 = x17_block()
    timing["X17"] = x17["runtime_s"]

    # ---- assemble
    out = {}
    x0_check = {}
    for k in DET:
        r0 = rows[k]["X0"]
        x0_check[k] = dict(union_tte_h=r0["tte_h"], lsoda_tte_h=lsoda[k]["tte_h"],
                           rel_diff_pct=100 * (r0["tte_h"] / lsoda[k]["tte_h"] - 1), event_union=r0["event"],
                           event_lsoda=lsoda[k]["event"], s_off_union_pct=r0["s_off_pct"],
                           s_off_lsoda_pct=lsoda[k]["s_off_pct"], h_s=H_A if DET[k][6] == "A" else H_B)
        lab, sch, Ta, soh, T0, s0, _ = DET[k]
        # X1
        t1 = l0_schedule(sch, Ta, _p(soh), s0)
        out[f"X1|{k}"] = dict(dtte_pct=100 * (t1 / r0["tte_h"] - 1), ds_off_pct=-r0["s_off_pct"], tte_h=t1,
                              s_off_pct=0.0, event="empty (L0)", tte_X0_h=r0["tte_h"], s_off_X0_pct=r0["s_off_pct"])
        for vid, r in rows[k].items():
            if vid == "X0":
                continue
            s0X = r0["s_off_pct"] if vid != "X12" else r0["s_off_pct"]  # X12: D2 output vs the D1 baseline output
            out[f"{vid}|{k}"] = dict(dtte_pct=100 * (r["tte_h"] / r0["tte_h"] - 1),
                                     ds_off_pct=r["s_off_pct"] - s0X, tte_h=r["tte_h"], s_off_pct=r["s_off_pct"],
                                     event=r["event"], tte_X0_h=r0["tte_h"], s_off_X0_pct=r0["s_off_pct"],
                                     event_X0=r0["event"])
    for key, per in ctmc_rows.items():
        b = per["X0"]
        mb = float(np.nanmean(b["tte"]))
        for vid, r in per.items():
            if vid == "X0":
                continue
            dif = r["tte"] - b["tte"]
            ok = np.isfinite(dif)
            out[f"{vid}|{key}"] = dict(
                dtte_pct=100 * (float(np.nanmean(r["tte"])) / mb - 1),
                dtte_se_pct=100 * float(np.std(dif[ok], ddof=1) / math.sqrt(ok.sum())) / mb,
                ds_off_pct=100 * float(np.nanmean(r["s_off"]) - np.nanmean(b["s_off"])),
                tte_h=float(np.nanmean(r["tte"])), sd_h=float(np.nanstd(r["tte"], ddof=1)),
                s_off_pct=100 * float(np.nanmean(r["s_off"])), tte_X0_h=mb, sd_X0_h=float(np.nanstd(b["tte"], ddof=1)),
                s_off_X0_pct=100 * float(np.nanmean(b["s_off"])),
                brownout_share=float(np.mean(r["event"] == "brownout")),
                brownout_share_X0=float(np.mean(b["event"] == "brownout")), n_paths=int(ok.sum()),
                paths="same as X0 (common random numbers)" if vid not in ("X15", "X16") else
                "own law, same seed")
    for k, v in x17["tte_at_2y"].items():
        out[f"X17|{k}"] = dict(dtte_pct=v["dtte_pct"], ds_off_pct=v["ds_off_pct"], tte_h=v["tte_h"],
                               s_off_pct=v["s_off_pct"], tte_X0_h=v["tte_X0_h"], s_off_X0_pct=v["s_off_X0_pct"],
                               event=v["event"], event_X0=v["event_X0"], dsoh_2y_pp=x17["dsoh_2y_pp"])
    # h convergence: variant effects at h = 1 s vs 2 s
    hconv = {}
    for k in conv_keys:
        for vid in rows_h1[k]:
            if vid == "X0":
                continue
            d2 = 100 * (rows[k][vid]["tte_h"] / rows[k]["X0"]["tte_h"] - 1)
            d1 = 100 * (rows_h1[k][vid]["tte_h"] / rows_h1[k]["X0"]["tte_h"] - 1)
            hconv[f"{vid}|{k}"] = dict(dtte_h2_pct=d2, dtte_h1_pct=d1, diff_pp=d2 - d1,
                                       tte_h1_h=rows_h1[k][vid]["tte_h"], tte_h2_h=rows[k][vid]["tte_h"])
    x0_ctmc = {key: dict(mean_tte_h=float(np.nanmean(per["X0"]["tte"])), sd_h=float(np.nanstd(per["X0"]["tte"],
                                                                                                   ddof=1)),
                         s_off_pct=100 * float(np.nanmean(per["X0"]["s_off"])),
                         brownout_share=float(np.mean(per["X0"]["event"] == "brownout")),
                         n_paths=int(np.sum(np.isfinite(per["X0"]["tte"]))))
               for key, per in ctmc_rows.items()}
    return dict(results=out, x0_check=x0_check, x0_ctmc=x0_ctmc, core_check=core_check, h_convergence=hconv,
                X17=x17, timing=timing, runtime_s=round(time.time() - t_all, 1))


# ------------------------------------------------------------------------------------------- table T4
T4_COLS = [("S4", "dtte_pct"), ("S7", "dtte_pct"), ("S12", "dtte_pct"), ("S13", "dtte_pct"),
           ("S9@-20/0.8", "dtte_pct"), ("S9@-20/0.8", "ds_off_pct"), ("S10a@-20/0.8", "dtte_pct"),
           ("S10a@-20/0.8", "ds_off_pct")]
TEX_LABEL = {
    "X1": "L0 energy bucket", "X2": "L1, isothermal", "X3": "L2 isothermal, no throttling",
    "X4": "constant current", "X5": r"2RC ($R\!\cdot\!Q$-scaled fit)", "X5n": r"2RC, $R_a{+}R_b=R_{1}$",
    "X6": r"$R_0$ split: PCM + cell", "X7": "no burst envelope", "X8": r"$V_\mathrm{sys}=V_\mathrm{cut}$",
    "X9": "KiBaM instead of RC", "X10": "RRC tail ODE", "X11a": r"$\gamma_b=1$", "X11b": r"$\gamma_b=2.9$",
    "X12": "SOC output D2", "X13": r"leakage, $k_\mathrm{lk}=0.04$", "X14": "audio on DSP",
    "X15": "lognormal sessions", "X16": "diurnal CTMC", "X17": r"calendar law $t^{0.75}$",
    "X18a": r"OCV stretch 4.40\,V", "X18b": r"OCV stretch 4.48\,V"}


def _fmt(v, nd=1):
    if v is None or not np.isfinite(v):
        return "--"
    s = f"{v:+.{nd}f}"
    return "0.0" if s in ("+0.0", "-0.0") else s.replace("-", "$-$")


def make_table(js: dict) -> str:
    """T4: ΔTTE (%) and Δs_off (pp) of every variant on six key scenarios (compact, booktabs)."""
    R = js["results"]
    L = [r"\begin{tabular}{@{}llrrrrrrrr@{}}", r"\toprule",
         r" & & \multicolumn{4}{c}{25\,\textdegree C, new: $\Delta$TTE (\%)} & "
         r"\multicolumn{4}{c}{$-20$\,\textdegree C, SOH 0.8} \\",
         r"\cmidrule(lr){3-6}\cmidrule(l){7-10}",
         r" & & & & & & \multicolumn{2}{c}{S9 hike} & \multicolumn{2}{c}{S10a social} \\",
         r"\cmidrule(lr){7-8}\cmidrule(l){9-10}",
         r"ID & Variant & S4 & S7 & S12 & S13 & $\Delta$TTE & $\Delta s_\mathrm{off}$ & $\Delta$TTE & "
         r"$\Delta s_\mathrm{off}$ \\", r"\midrule"]
    for x in ORDER:
        cells = []
        for sk, fld in T4_COLS:
            key = f"{x}|{sk}"
            if x == "X17":
                key = {"S4": "X17|S4@2y", "S12": "X17|S12@2y"}.get(sk, "")
            r = R.get(key)
            cells.append(_fmt(r[fld]) if r else "--")
        L.append(f"{x} & {TEX_LABEL[x]} & " + " & ".join(cells) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}", ""]
    return "\n".join(L)


# ------------------------------------------------------------------------------------------- figure
FIG_SC = ["S1", "S3", "S4", "S5b", "S6", "S7@15", "S7", "S7@35", "S8a", "CM", "S11a@-110", "S16", "S12", "S13",
          "S14", "S9@-10/0.8", "S9@-20/0.8", "S10a@-10/0.8", "S10a@-20/0.8", "S10b@0/0.8"]
FIG_SC_LABEL = {"S7@15": "S7 15°C", "S7": "S7 25°C", "S7@35": "S7 35°C", "S11a@-110": "S11a −110",
                "S9@-10/0.8": "S9 −10°C", "S9@-20/0.8": "S9 −20°C", "S10a@-10/0.8": "S10a −10°C",
                "S10a@-20/0.8": "S10a −20°C", "S10b@0/0.8": "S10b 0°C"}
FIG_ROWS = [x for x in ORDER if x not in ("X12", "X17")]  # X12 changes no TTE at s0 = 100 %; X17 is aging-only


def make_figure(js: dict):
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    from style import INK, INK2, MUTED, W1, apply, cmap_div, save

    apply()
    R = js["results"]
    M = np.full((len(FIG_ROWS), len(FIG_SC)), np.nan)
    for i, x in enumerate(FIG_ROWS):
        for j, s in enumerate(FIG_SC):
            r = R.get(f"{x}|{s}")
            if r is not None:
                M[i, j] = r["dtte_pct"]
    lim = 20.0
    fig, ax = plt.subplots(figsize=(W1, 3.55))
    cm = cmap_div.copy()
    cm.set_bad("#ffffff")
    im = ax.imshow(np.ma.masked_invalid(np.clip(M, -lim, lim)), cmap=cm, norm=TwoSlopeNorm(0.0, -lim, lim),
                   aspect="auto", interpolation="nearest")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            if not np.isfinite(v):
                ax.text(j, i, "·", ha="center", va="center", fontsize=5, color=MUTED)
                continue
            if abs(v) < 0.5:  # selective labels: |ΔTTE| < 0.5 % is left to the colour
                continue
            txt = f"{v:.0f}" if abs(v) >= 9.5 else f"{v:.1f}"
            ax.text(j, i, txt.replace("-", "−"), ha="center", va="center", fontsize=4.9,
                    color="white" if abs(v) > 0.6 * lim else INK)
    ax.set_xticks(range(len(FIG_SC)))
    ax.set_xticklabels([FIG_SC_LABEL.get(s, s) for s in FIG_SC], rotation=55, ha="right", fontsize=6.4)
    ax.set_yticks(range(len(FIG_ROWS)))
    ax.set_yticklabels([f"{x}  {LABELS[x]}" for x in FIG_ROWS], fontsize=6.2)
    ax.tick_params(length=0)
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
    for x in (FIG_SC.index("S12") - 0.5, FIG_SC.index("S9@-10/0.8") - 0.5):
        ax.axvline(x, color="white", lw=2.0)
    cb = fig.colorbar(im, ax=ax, fraction=0.022, pad=0.01, extend="both")
    cb.ax.tick_params(labelsize=6.5, length=2)
    cb.outline.set_visible(False)
    cb.set_label("ΔTTE vs X0 (%)", fontsize=7, color=INK2)
    ax.set_title("ΔTTE of each structural variant relative to the baseline L2/L3 model (X0)", fontsize=8.5,
                 loc="left", color=INK)
    npth = js.get("meta", {}).get("n_paths", N_PATHS)
    ax.text(-0.42, -0.215, f"Colour clipped at ±{lim:.0f} %; numbers exact, shown where |ΔTTE| ≥ 0.5 %. Column blocks: "
                         f"constant loads (25 °C) | schedule and CTMC days (S13/S14: mean of {npth} paths) |\n"
                         "cold, aged packs (SOH 0.8). · = not applicable (usage-law variants). X12 (SOC output "
                         "definition) and X17 (aging law) are reported in Table T4.",
             fontsize=5.8, color=MUTED, ha="left", va="top", linespacing=1.3, transform=ax.transAxes)
    save(fig, "fig_structural")


# ------------------------------------------------------------------------------------------- digest
def write_summary(js: dict):
    """results/global_summary.md: AN9 + AN10 digest for the paper writer (reads results/sobol.json if present)."""
    try:
        from summary_global import build_summary
    except ImportError:  # digest module not present (development)
        return
    (RESULTS / "global_summary.md").write_text(build_summary(js))


# ------------------------------------------------------------------------------------------- entry
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--paths", type=int, default=N_PATHS, help="CTMC paths per profile (default 200)")
    ap.add_argument("--figure-only", action="store_true", help="re-draw the figure/table from structural.json")
    a = ap.parse_args(argv)
    if a.figure_only:
        js = json.loads((RESULTS / "structural.json").read_text())
    else:
        js = compute(a.paths)
        js = dict(meta=meta("run_structural.py", seed=SEED, solver="union model (analysis/structural_variants.py) "
                                                                    "on rk4_batch; X1 exact; X17 core LSODA",
                            h_s=H_A, h_s_standby=H_B, h_s_ctmc=H_C, n_paths=a.paths, runtime_s=js.pop("runtime_s"),
                            timing_s=js.pop("timing"), variant_definitions={k: LABELS[k] for k in ORDER},
                            key_format="'{X}|{S}' -> dtte_pct = 100 (TTE_X/TTE_X0 − 1) at the same solver and h; "
                                       "ds_off_pct = s_off_X − s_off_X0 in percentage points (X12: D2 minus D1); "
                                       "CTMC days: means over paths, dtte_se_pct from paired differences",
                            scenarios={k: DET[k][0] for k in DET} | {k: v[0] for k, v in CTMC_SC.items()}),
                  **js)
        write_json("structural.json", js)
    write_tex("tab_structural.tex", make_table(js))
    make_figure(js)
    write_summary(js)
    print("wrote results/structural.json, paper/tables/tab_structural.tex, fig_structural")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
