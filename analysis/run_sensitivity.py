"""AN8: local elasticities and one-at-a-time tornado (FINAL_SPEC §10.1, §10.2) -> results/sensitivity.json,
fig_tornado.pdf, paper/tables/tab_elasticities.tex.

Local elasticities (§10.1)
    ε_θ = ∂ln TTE/∂ln θ by central differences at ±5 % (the [PROTO] formula
    [ln TTE(1.05θ) − ln TTE(0.95θ)]/[ln 1.05 − ln 0.95]).  Voltages (δU, V_cut, Δ_sys) are reported as % TTE per
    +10 mV (central ±10 mV) and the skin limit T_lim as % TTE per +1 K (central ±1 K), because a ±5 % change of a
    voltage level or of a Celsius temperature is not meaningful.  Cases: S4 at 25 °C (new), S4 at −10 °C on an
    SOH-0.8 pack, S7 at 25 °C: L2 with the LSODA reference solver (§6.2), the same solver as [PROTO].  Also the two
    smooth metrics of §12.1: E_u/P̄ of the S12 day (quasi-steady powers) and the S13 Dynkin mean (F6).
    Every UQ_SPEC factor that acts within one discharge is included (aging-rate constants are not); the 12 factors
    of the [PROTO] table are compared with it (``proto_check``).

One-at-a-time tornado (§10.2)
    Each factor moved to its §4 range endpoints (others at default): S4 (25 °C, new), S7 (25 °C), S10a (S4 at
    −10 °C, SOH 0.8) on the L2 TTE (LSODA); S12 on E_u/P̄; S13 on the Dynkin mean.  User-behaviour and environment
    factors on the same axis: brightness b ± 0.2 (screen-on states), dark theme (APL × 0.25 on RD, SO, SO-L, NV),
    refresh 60 / 120 Hz (screen-on states), background wake-up rate × 0.5 / × 5 (standby states), RSRP −115 /
    −85 dBm, data RAT Wi-Fi / 5G NSA (LTE stored as the middle level; navigation and tracking stay cellular),
    T_amb −10 / 35 °C (soaked start) and enclosure k_env 0.8 / 1.5.  "lo"/"hi" are the lower/upper factor values
    (for RSRP lo = −115 dBm; for the RAT lo = Wi-Fi, hi = 5G NSA; for the dark theme lo = dark, hi = default).

JSON keys (§14): ``elasticities."{S}|{T}|{SOH}".{param}``, ``oat."{S}".{factor}.{lo_h, hi_h}`` (plus lo/hi values,
ΔTTE %, kind, event), ``elasticity_units``, ``proto_check``, ``base``.

Run:  cd /home/user/MCM && PYTHONPATH=src:analysis python3 analysis/run_sensitivity.py [--workers 3]
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import math
import multiprocessing as mp
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace

import numpy as np

from basics_common import meta, write_json, write_tex
from phonebatt import cell
from phonebatt.analytic import ctmc_mean_tte, schedule_mean_power
from phonebatt.model import run
from phonebatt.params import RP63, UQ_SPEC, aged
from phonebatt.usage import CTMC, PRIMITIVES, UI_STATES, Schedule, s12_schedule

# --------------------------------------------------------------------------------------------- factors
# Hardware factors: every UQ_SPEC entry that acts within one discharge (aging-rate constants, charging and
# calendar factors are excluded; k_env is treated as an environment/behaviour factor below).
EXCLUDE = {"eta_c", "dT_ch", "eta_charger", "k_c", "E_cal", "g_hi", "k_HV", "w_L", "E_H", "E_L", "beta_L", "k_env"}
HW = [k for k in UQ_SPEC if k not in EXCLUDE]
VOLT = {"dU", "V_cut", "d_sys"}  # % TTE per +10 mV
TEMP = {"T_lim_C"}  # % TTE per +1 K
PROTO_ELAST = {  # FINAL_SPEC §10.1 [PROTO] table
    "S4|25|1.0": dict(Q_n=1.000, P_SoC_max=-0.406, p_c=-0.119, P_d0=-0.114, P_awake=-0.086, R0_ref=-0.007,
                      Ea_R0=0.001, R1_ref=-0.001, m_b=0.000, R_skin=0.000, R_th=0.000, C_th=0.000),
    "S4|-10|0.8": dict(Q_n=1.000, P_SoC_max=-0.431, p_c=-0.127, P_d0=-0.121, P_awake=-0.091, R0_ref=-0.274,
                       Ea_R0=-0.236, R1_ref=-0.018, m_b=-0.190, R_skin=0.000, R_th=0.030, C_th=0.000),
    "S7|25|1.0": dict(Q_n=1.009, P_SoC_max=-0.180, p_c=-0.012, P_d0=-0.010, P_awake=-0.007, R0_ref=-0.012,
                      Ea_R0=0.003, R1_ref=-0.001, m_b=0.000, R_skin=0.767, R_th=-0.007, C_th=-0.009),
}

# Cases: (kind, primitive or None, T_amb, SOH)
CASES = {
    "S4|25|1.0": ("lsoda", "SO", 25.0, 1.0),
    "S4|-10|0.8": ("lsoda", "SO", -10.0, 0.8),
    "S7|25|1.0": ("lsoda", "GM", 25.0, 1.0),
    "S12-EuPbar|25|1.0": ("s12", None, 25.0, 1.0),
    "S13-Dynkin|25|1.0": ("s13", None, 25.0, 1.0),
}
OAT_CASES = {"S4": "S4|25|1.0", "S7": "S7|25|1.0", "S10a": "S4|-10|0.8", "S12": "S12-EuPbar|25|1.0",
             "S13": "S13-Dynkin|25|1.0"}
OAT_TITLES = {"S4": "S4 social, 25 °C, new", "S7": "S7 gaming, 25 °C, new", "S10a": "S10a social, −10 °C, SOH 0.8",
              "S12": "S12 day, E_u/P̄ (25 °C)", "S13": "S13 typical day, Dynkin mean (25 °C)"}

# Behaviour / environment factors: id -> (label, lo description, hi description)
BEHAV = {
    "brightness": ("Brightness", "−0.2", "+0.2"),
    "dark_theme": ("Dark theme", "APL × 0.25", "off"),
    "refresh": ("Refresh rate", "60 Hz", "120 Hz"),
    "bg_rate": ("Background wake-ups", "× 0.5", "× 5"),
    "rsrp": ("Cellular RSRP", "−115 dBm", "−85 dBm"),
    "rat": ("Data RAT", "Wi-Fi", "5G NSA"),
    "tamb": ("Ambient", "−10 °C", "35 °C"),
    "k_env": ("Enclosure k_env", "0.8", "1.5"),
}


def _scr(fn):
    return lambda n, u: fn(u) if np.all(np.asarray(u.s_scr) > 0) else u


def _rat(level: str):
    def f(n, u):
        if u.rho == "none" or n in ("NV", "TRK"):
            return u
        if level == "wifi":
            return replace(u, w=1, rho="wifi", c="lte")
        return replace(u, w=0, rho="cell", c="lte" if level == "lte" else "nr")
    return f


def behaviour(fid: str, level: str, Tamb: float):
    """(map_fn or None, T_amb, parameter overrides) of a behaviour/environment factor at level lo/mid/hi."""
    lo = level == "lo"
    if fid == "brightness":
        d = -0.2 if lo else 0.2
        return _scr(lambda u: replace(u, b=float(np.clip(u.b + d, 0.0, 1.0)))), Tamb, {}
    if fid == "dark_theme":
        return ((lambda n, u: replace(u, m_APL=0.25) if n in UI_STATES else u) if lo else None), Tamb, {}
    if fid == "refresh":
        return _scr(lambda u: replace(u, f_r=60.0 if lo else 120.0)), Tamb, {}
    if fid == "bg_rate":
        m = 0.5 if lo else 5.0
        return (lambda n, u: replace(u, lam=u.lam * m) if n.startswith("SB") else u), Tamb, {}
    if fid == "rsrp":
        return (lambda n, u: replace(u, RSRP=-115.0 if lo else -85.0)), Tamb, {}
    if fid == "rat":
        return _rat({"lo": "wifi", "mid": "lte", "hi": "nr"}[level]), Tamb, {}
    if fid == "tamb":
        return None, (-10.0 if lo else 35.0), {}
    if fid == "k_env":
        return None, Tamb, {"k_env": 0.8 if lo else 1.5}
    raise KeyError(fid)


# --------------------------------------------------------------------------------------------- evaluation
def evaluate(job: dict) -> dict:
    """Metric of one case under parameter overrides and/or one behaviour factor level."""
    kind, prim, Tamb, soh = CASES[job["case"]]
    p = RP63 if soh == 1.0 else aged(RP63, 1.0 - soh)
    ov = dict(job.get("ov", {}))
    map_fn = None
    if job.get("fid"):
        map_fn, Tamb, ov2 = behaviour(job["fid"], job["level"], Tamb)
        ov.update(ov2)
    if ov:
        p = replace(p, **ov)
    t0 = time.time()
    if kind == "lsoda":
        u = PRIMITIVES[prim]
        if map_fn is not None:
            u = map_fn(prim, u)
        r = run(Schedule.constant(u), Tamb, p=p)
        out = dict(value=r.tte_h, event=r.event, s_off_pct=100 * r.s_off)
    elif kind == "s12":
        Pb = schedule_mean_power(s12_schedule(), p, Tamb, map_fn=map_fn)
        out = dict(value=float(cell.usable_energy(p)) / Pb, event="", s_off_pct=math.nan, Pbar_W=Pb)
    else:
        out = dict(value=float(ctmc_mean_tte(CTMC("typical"), p, Tamb, map_fn=map_fn)), event="", s_off_pct=math.nan)
    out["runtime_s"] = time.time() - t0
    return out


def default_of(name: str) -> float:
    return float(getattr(RP63, name))


def elasticity_jobs() -> list[dict]:
    jobs = []
    for case in CASES:
        jobs.append(dict(tag=("E", case, "base"), case=case))
        for k in HW:
            v = default_of(k)
            if k in VOLT:
                pair = (v - 0.010, v + 0.010)
            elif k in TEMP:
                pair = (v - 1.0, v + 1.0)
            else:
                pair = (0.95 * v, 1.05 * v)
            for s, val in zip(("m", "p"), pair):
                jobs.append(dict(tag=("E", case, k, s), case=case, ov={k: val}))
    return jobs


def oat_jobs() -> list[dict]:
    jobs = []
    for sid, case in OAT_CASES.items():
        for k in HW:
            _, lo, _, hi = UQ_SPEC[k]
            jobs.append(dict(tag=("O", sid, k, "lo"), case=case, ov={k: float(lo)}))
            jobs.append(dict(tag=("O", sid, k, "hi"), case=case, ov={k: float(hi)}))
        for fid in BEHAV:
            for lev in ("lo", "hi") + (("mid",) if fid == "rat" else ()):
                jobs.append(dict(tag=("O", sid, fid, lev), case=case, fid=fid, level=lev))
    return jobs


def _run(job):
    return job["tag"], evaluate(job)


def run_all(jobs: list[dict], workers: int) -> dict:
    if workers <= 1:
        return dict(_run(j) for j in jobs)
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("fork")) as ex:
        return dict(ex.map(_run, jobs, chunksize=4))


# --------------------------------------------------------------------------------------------- figure, table
SYM = {"Q_n": r"$Q_n$", "R0_ref": r"$R_{0,\mathrm{ref}}$", "A_R": r"$A_R$", "Ea_R0": r"$E_{a,R0}$",
       "R1_ref": r"$R_{1,\mathrm{ref}}$", "tau1_ref": r"$\tau_{1,\mathrm{ref}}$", "Ea_R1": r"$E_{a,R1}$",
       "Ea_tau": r"$E_{a,\tau}$", "k_S": r"$k_S$", "dU": r"$\delta U$", "V_cut": r"$V_{\mathrm{cut}}$",
       "d_sys": r"$\Delta_{\mathrm{sys}}$", "C_th": r"$C_{\mathrm{th}}$", "R_th": r"$R_{\mathrm{th}}$",
       "R_skin": r"$R_{\mathrm{skin}}$", "T_lim_C": r"$T_{\mathrm{lim}}$", "w_thr": r"$w_{\mathrm{thr}}$",
       "P_fl": r"$P_{\mathrm{fl}}$", "P_susp": r"$P_{\mathrm{susp}}$", "P_awake": r"$\Delta P_{\mathrm{awake}}$",
       "E_wake": r"$E_{\mathrm{wake}}$", "P_d0": r"$P_{d0}$", "p_c": r"$p_c$", "gamma_b": r"$\gamma_b$",
       "k_RR": r"$\kappa_{RR}$", "P_SoC_max": r"$P_{\mathrm{SoC,max}}$", "v0": r"$v_0$", "P_G": r"$P_G$",
       "P_cam": r"$\Delta P_{\mathrm{cam}}$", "P_camf": r"$\Delta P_{\mathrm{camf}}$", "P_aud": r"$\Delta P_{\mathrm{aud}}$",
       "P_BT": r"$\Delta P_{\mathrm{BT}}$", "m_b": r"$m_b$", "k_ul": r"$\kappa_{ul}$", "gamma_sig": r"$\gamma_{\mathrm{sig}}$",
       "gamma_R": r"$\gamma_R$"}
RADIO_SYM = {"P_idle": r"$P_{\mathrm{idle}}$", "P_tail": r"$P_{\mathrm{tail}}$", "T_tail": r"$T_{\mathrm{tail}}$",
             "P_act": r"$P_{\mathrm{act}}$", "T_act": r"$T_{\mathrm{act}}$", "E_pro": r"$E_{\mathrm{pro}}$",
             "e_b": r"$e_b$", "P_rt": r"$P_{\mathrm{rt}}$"}
UNIT_FMT = {"R0_ref": (1e3, "mΩ", 0), "R1_ref": (1e3, "mΩ", 0), "Ea_R0": (1e-3, "kJ/mol", 0),
            "Ea_R1": (1e-3, "kJ/mol", 0), "Ea_tau": (1e-3, "kJ/mol", 0), "V_cut": (1, "V", 2), "d_sys": (1e3, "mV", 0),
            "dU": (1e3, "mV", 0), "P_SoC_max": (1, "W", 0), "p_c": (97.4e3, "mW", 0), "P_d0": (1e3, "mW", 0),
            "P_awake": (1e3, "mW", 0), "P_susp": (1e3, "mW", 0), "C_th": (1, "J/K", 0), "R_th": (1, "K/W", 0),
            "R_skin": (1, "K/W", 1), "T_lim_C": (1, "°C", 0), "Q_n": (1, "Ah", 2), "m_b": (1, "", 1),
            "A_R": (1, "", 1), "tau1_ref": (1, "s", 0), "w_thr": (1, "K", 1), "P_fl": (1, "W", 1),
            "gamma_b": (1, "", 1), "v0": (1, "", 2), "gamma_R": (1, "", 1), "k_S": (1, "", 1), "E_wake": (1, "J", 2),
            "k_RR": (1e3, "mW/Hz", 1), "P_G": (1e3, "mW", 0), "P_cam": (1, "W", 2), "gamma_sig": (1, "", 2)}


def hw_label(k: str) -> str:
    _, lo, _, hi = UQ_SPEC[k]
    sym = SYM.get(k)
    if sym is None:
        rat, name = k.split("_", 1)
        sym = {"wifi": "Wi-Fi ", "lte": "LTE ", "nr": "5G "}[rat] + RADIO_SYM[name]
    sc, unit, d = UNIT_FMT.get(k, (1, "", 3))
    if k == "p_c":
        rng = f"{lo * sc:.0f}–{hi * sc:.0f} mW"
    else:
        rng = f"{lo * sc:.{d}f}–{hi * sc:.{d}f}" + (f" {unit}" if unit else "")
    return f"{sym}  {rng}"


def make_figure(oat: dict, base: dict) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from style import BLUE, INK, INK2, MUTED, RED, W1, apply, panel_label, save

    apply()
    fig, axes = plt.subplots(1, 2, figsize=(W1, 3.15))
    n_show = 13
    for ax, sid, tag in zip(axes, ("S4", "S10a"), ("(a)", "(b)")):
        d = oat[sid]
        rows = sorted(d.items(), key=lambda kv: -max(abs(kv[1]["lo_pct"]), abs(kv[1]["hi_pct"])))[:n_show]
        rows = rows[::-1]
        y = np.arange(len(rows))
        for yy, (k, v) in zip(y, rows):
            ax.barh(yy + 0.19, v["lo_pct"], height=0.36, color=BLUE, lw=0)
            ax.barh(yy - 0.19, v["hi_pct"], height=0.36, color=RED, lw=0)
            for val, off in ((v["lo_pct"], 0.19), (v["hi_pct"], -0.19)):
                if abs(val) >= 6:
                    ax.text(val + (1.0 if val > 0 else -1.0), yy + off, f"{val:+.0f}", va="center",
                            ha="left" if val > 0 else "right", fontsize=6.3, color=INK2)
        labels = []
        for k, v in rows:
            if v["kind"] == "behaviour":
                lab = BEHAV[k][0] + f"  {BEHAV[k][1]} / {BEHAV[k][2]}"
            else:
                lab = hw_label(k)
            labels.append(lab)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=6.6)
        for tl, (k, v) in zip(ax.get_yticklabels(), rows):
            if v["kind"] == "behaviour":
                tl.set_fontstyle("italic")
                tl.set_color(INK)
        ax.axvline(0, color=MUTED, lw=0.8)
        ax.grid(axis="y", visible=False)
        vals = [x for _, v in rows for x in (v["lo_pct"], v["hi_pct"])] + [0.0]
        span = max(vals) - min(vals)
        ax.set_xlim(min(vals) - 0.13 * span - 2, max(vals) + 0.13 * span + 2)
        ax.set_xlabel("ΔTTE vs default (%)")
        ax.set_ylim(-0.7, len(rows) - 0.3)
        ax.set_title(f"{OAT_TITLES[sid]}: {base[sid]:.1f} h", loc="left", fontsize=8.3, fontweight="normal",
                     color=INK2)
        panel_label(ax, tag)
    hand = [Patch(facecolor=BLUE, label="factor at its low end"), Patch(facecolor=RED, label="factor at its high end")]
    fig.legend(handles=hand, loc="lower left", ncol=2, bbox_to_anchor=(0.0, -0.02), fontsize=7.2)
    fig.text(1.0, 0.0, "italic: user behaviour / environment (low / high value shown)\n"
             "upright: hardware parameter at its §4 range ends", ha="right", va="bottom", fontsize=6.6, color=INK2)
    fig.subplots_adjust(wspace=1.05, bottom=0.22)
    save(fig, "fig_tornado")


def write_table(el: dict) -> None:
    rows = ["Q_n", "P_SoC_max", "v0", "p_c", "P_d0", "k_RR", "P_awake", "R0_ref", "Ea_R0", "R1_ref", "A_R",
            "gamma_R", "m_b", "R_skin", "R_th", "C_th", "T_lim_C", "V_cut", "d_sys", "dU"]
    cases = ["S4|25|1.0", "S4|-10|0.8", "S7|25|1.0", "S12-EuPbar|25|1.0", "S13-Dynkin|25|1.0"]
    sym = {k: SYM.get(k, k) for k in rows}
    lines = [r"\begin{tabular}{lrrrrr}", r"\toprule",
             r"Parameter & S4, 25\,$^\circ$C & S4, $-10\,^\circ$C, SOH 0.8 & S7, 25\,$^\circ$C & S12 $E_u/\bar P$ & "
             r"S13 Dynkin \\", r"\midrule"]
    for k in rows:
        unit = r" (\%/10\,mV)" if k in VOLT else (r" (\%/K)" if k in TEMP else "")
        vals = []
        for c in cases:
            v = el[c][k]
            vals.append(f"${v:+.3f}$" if abs(v) >= 5e-4 else "$0.000$")
        lines.append(sym[k] + unit + " & " + " & ".join(vals) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    tex = ("% Local elasticities d ln TTE / d ln theta (central +-5 %; voltages % per +10 mV, T_lim % per +1 K).\n"
           "% L2 with LSODA (S4, S7); smooth metrics of §12.1 (S12 E_u/Pbar, S13 Dynkin mean).  analysis/run_sensitivity.py\n"
           + "\n".join(lines) + "\n")
    write_tex("tab_elasticities.tex", tex)


# --------------------------------------------------------------------------------------------- main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=int(os.environ.get("SENS_WORKERS", "3")))
    a = ap.parse_args(argv)
    T0 = time.time()
    jobs = elasticity_jobs() + oat_jobs()
    res = run_all(jobs, a.workers)
    t_comp = time.time() - T0

    # elasticities
    el, el_detail = {}, {}
    units = {k: ("pct_TTE_per_+10mV" if k in VOLT else "pct_TTE_per_+1K" if k in TEMP else "dlnTTE/dlntheta")
             for k in HW}
    base = {c: res[("E", c, "base")]["value"] for c in CASES}
    for c in CASES:
        el[c], el_detail[c] = {}, {}
        T0c = base[c]
        for k in HW:
            tm, tp = res[("E", c, k, "m")]["value"], res[("E", c, k, "p")]["value"]
            if k in VOLT or k in TEMP:
                e = 100.0 * (tp - tm) / (2.0 * T0c)
            else:
                e = (math.log(tp) - math.log(tm)) / (math.log(1.05) - math.log(0.95))
            el[c][k] = e
            el_detail[c][k] = dict(minus_h=tm, plus_h=tp, ev_minus=res[("E", c, k, "m")]["event"],
                                   ev_plus=res[("E", c, k, "p")]["event"])
    proto = {}
    for c, d in PROTO_ELAST.items():
        proto[c] = {k: dict(model=el[c][k], proto=v, diff=el[c][k] - v) for k, v in d.items()}
    max_dev = max(abs(v["diff"]) for d in proto.values() for v in d.values())

    # OAT
    oat, oat_base = {}, {}
    for sid, case in OAT_CASES.items():
        b = base[case]
        oat_base[sid] = b
        d = {}
        for k in HW:
            lo, hi = res[("O", sid, k, "lo")], res[("O", sid, k, "hi")]
            d[k] = dict(kind="hardware", lo=float(UQ_SPEC[k][1]), hi=float(UQ_SPEC[k][3]), default=default_of(k),
                        lo_h=lo["value"], hi_h=hi["value"], lo_pct=100 * (lo["value"] / b - 1),
                        hi_pct=100 * (hi["value"] / b - 1), lo_event=lo["event"], hi_event=hi["event"],
                        lo_s_off_pct=lo["s_off_pct"], hi_s_off_pct=hi["s_off_pct"])
        for fid, (lab, lod, hid) in BEHAV.items():
            lo, hi = res[("O", sid, fid, "lo")], res[("O", sid, fid, "hi")]
            d[fid] = dict(kind="behaviour", label=lab, lo=lod, hi=hid, lo_h=lo["value"], hi_h=hi["value"],
                          lo_pct=100 * (lo["value"] / b - 1), hi_pct=100 * (hi["value"] / b - 1),
                          lo_event=lo["event"], hi_event=hi["event"], lo_s_off_pct=lo["s_off_pct"],
                          hi_s_off_pct=hi["s_off_pct"])
            if fid == "rat":
                mid = res[("O", sid, fid, "mid")]
                d[fid].update(mid="LTE", mid_h=mid["value"], mid_pct=100 * (mid["value"] / b - 1))
        oat[sid] = d
    ranking = {}
    for sid, d in oat.items():
        rk = sorted(d, key=lambda k: -max(abs(d[k]["lo_pct"]), abs(d[k]["hi_pct"])))
        ranking[sid] = [dict(factor=k, kind=d[k]["kind"], lo_pct=d[k]["lo_pct"], hi_pct=d[k]["hi_pct"],
                             swing_pct=abs(d[k]["hi_pct"] - d[k]["lo_pct"])) for k in rk[:15]]

    make_figure(oat, oat_base)
    write_table(el)
    out = dict(
        meta=meta("run_sensitivity.py", seed=None, solver="LSODA reference (model.run, rtol 1e-7) for S4/S7/S10a; "
                  "analytic quasi-steady metrics for S12 (E_u/Pbar) and S13 (Dynkin mean, F6)", h_s=None,
                  runtime_compute_s=t_comp, workers=a.workers, n_runs=len(jobs)),
        base=dict(cases={c: dict(value_h=base[c], kind=CASES[c][0], T_amb_C=CASES[c][2], soh=CASES[c][3],
                                 event=res[("E", c, "base")]["event"]) for c in CASES}, oat=oat_base),
        elasticities=el, elasticity_units=units, elasticity_detail=el_detail,
        proto_check=dict(table=proto, max_abs_diff=max_dev), oat=oat, oat_ranking=ranking,
        behaviour_factors={k: dict(label=v[0], lo=v[1], hi=v[2]) for k, v in BEHAV.items()},
        notes=["elasticity step: central +-5 % (PROTO formula); voltages +-10 mV reported as % TTE per +10 mV; "
               "T_lim +-1 K reported as % TTE per +1 K",
               "S12/S13 metrics are smooth L0 quantities with quasi-steady (F5) powers: resistance and cut-off "
               "parameters have no effect on them by construction",
               "RAT factor: every state with a data path is forced to Wi-Fi (lo), LTE (mid) or 5G NSA (hi); "
               "navigation (NV) and tracking (TRK) stay cellular"],
        runtime_s=time.time() - T0)
    write_json("sensitivity.json", out)
    print(f"[sens] {len(jobs)} runs, compute {t_comp:.0f} s; max |elasticity - PROTO| = {max_dev:.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
