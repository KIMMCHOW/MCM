"""AN15: generalisation to other portable devices - dimensionless groups Π1-Π5 (FINAL_SPEC §13, §14).

The groups are written per cell-equivalent, which makes them independent of the pack size and series count
(A12: R ∝ (N_s/N_p)/Q_cell, so R·Q per cell-equivalent is the material constant):

    Π1 = P R_Σ / Ū²                   = c · RQ_Σ / Ū                     (c = P/E, the C-rate in h⁻¹)
    Π2 = P_dem / P_sust,              P_sust = (T_lim − ΔT_thr − T_amb)/(k_env R_skin)
    Π3 = τ_th / TTE                   (TTE under the heaviest sustained load, E_u/min(P_dem, P_sust))
    Π4 = ΔP_b / (P_cut,sys − P),      P_cut,sys = V_sys (U(z,T) − V_sys)/R0  →  c_b Ū / (V_sys(U − V_sys)/RQ0 − c Ū)
    Π5 = √(2 τ_c/μ) · CV_P            (= CV of the TTE, F7)

RP-63 values are computed from the model (S4 social load for Π1 and Π4, gaming for Π2 and Π3, the typical CTMC for
Π5).  P_cut,sys follows the §7 F4 table convention (no RC polarisation, Ũ = U), which reproduces the §13 [PROTO]
values 0.11 and 1.2; the variant with v1 = R1 P/U is reported as well.  The resistance and OCV *shapes* (SOC, T and
aging multipliers) of RP-63 are applied to every class (same Li-ion chemistry), so the cold-aged Π4 of a class uses
R0(z, T, SOH)/R0_ref and U(z, T) of RP-63.

Device classes (tablet, laptop, smartwatch, earbuds, e-reader, drone) are ILLUSTRATIVE ORDERS OF MAGNITUDE
[unverified]: no device data were collected for them (§13).  Each input has a (mid, lo, hi) range; each Π is
evaluated at the mid values and at every lo/hi corner of its inputs, and the min/max over corners is the range.

Writes
* results/generalization.json: ``RP63.{Pi1..Pi5, Pi4_cold_aged, ...}``, ``classes.{name}.{Pi1..Pi5 (mid, lo, hi),
  Pi4_cold_aged, floor_elasticity, inputs, notes}``, ``scaling_rules``, ``ereader_example``, ``meta``
* paper/figures/fig_generalization.pdf (+ png): (a) (Π1, Π4) map, (b) (Π2, Π5) map

Run:  cd /home/user/MCM && PYTHONPATH=src:analysis python3 analysis/run_generalization.py   (~10 s)
"""
from __future__ import annotations

import itertools
import json
import math
import time

import matplotlib.pyplot as plt
import numpy as np

from basics_common import RESULTS, meta, write_json
from phonebatt import CTMC, RP63, aged, run
from phonebatt import analytic as A
from phonebatt import cell
from phonebatt.params import K0, T_REF
from phonebatt.power import demand_power
from phonebatt.usage import PRIMITIVES, Schedule
from style import BLUE, INK, INK2, MUTED, NEUTRAL, RED, W1, apply, panel_label, save

T_COLD_K = -20.0 + K0
SOH_COLD = 0.8
Z_NOM, Z_COLD = 0.5, 0.2
PROTO = dict(Pi1=0.005, Pi2=1.19, Pi3=0.02, Pi4=0.11, Pi4_cold_aged=1.2, Pi5=0.175)

# Illustrative device classes: (mid, lo, hi) per input [unverified orders of magnitude].
# E_Wh pack energy; RQ_sigma = R_Σ·Q per cell-equivalent (Ω·Ah; RP-63 0.329 incl. PCM); P_typ typical active power;
# P_dem heaviest sustained demand; P_sust thermally sustainable power; dPb sub-second burst above P_typ;
# tau_th_min thermal time constant; CV_P, tau_c_h, mu_h usage variability, correlation time and typical TTE;
# floor_share = standby/always-on share of typical-use energy (F1 elasticity of TTE to the floor = −share).
CLASSES = {
    "tablet": dict(
        inputs=dict(E_Wh=(30, 20, 40), RQ_sigma=(0.33, 0.2, 0.5), P_typ=(3.0, 2.0, 5.0), P_dem=(8, 6, 12),
                    P_sust=(9, 6, 13), dPb=(4, 3, 6), tau_th_min=(10, 5, 20), CV_P=(1.5, 1.0, 2.0),
                    tau_c_h=(0.5, 0.25, 1.0), mu_h=(40, 20, 80), floor_share=(0.2, 0.1, 0.35)),
        notes="display-dominated; 2.5-4x the phone's surface relaxes throttling (smaller Π2); same cell chemistry "
              "and similar C-rates as a phone, so Π1 and Π4 stay phone-like"),
    "laptop": dict(
        inputs=dict(E_Wh=(55, 40, 100), RQ_sigma=(0.2, 0.1, 0.33), P_typ=(8, 5, 15), P_dem=(45, 25, 90),
                    P_sust=(25, 15, 45), dPb=(20, 10, 40), tau_th_min=(15, 5, 30), CV_P=(0.5, 0.3, 0.8),
                    tau_c_h=(0.5, 0.2, 1.0), mu_h=(8, 5, 12), floor_share=(0.05, 0.02, 0.1)),
        notes="multi-cell (2-4 S; N_s drops out of the per-cell groups); large cells with a small PCM share give a "
              "small Π1; CPU/GPU turbo exceeds the sustainable power (Π2 > 1); continuous sessions give a small "
              "CV_P"),
    "smartwatch": dict(
        inputs=dict(E_Wh=(1.1, 0.8, 2.0), RQ_sigma=(0.33, 0.2, 0.6), P_typ=(0.1, 0.05, 0.2),
                    P_dem=(0.6, 0.4, 1.0), P_sust=(0.5, 0.3, 0.8), dPb=(0.5, 0.3, 1.0), tau_th_min=(1.0, 0.5, 3.0),
                    CV_P=(2.0, 1.0, 3.0), tau_c_h=(0.1, 0.03, 0.3), mu_h=(30, 18, 48), floor_share=(0.6, 0.4, 0.8)),
        notes="tiny C_th (τ_th about a minute); standby and BLE floor dominate, so the floor elasticity approaches "
              "-1; LTE bursts are large relative to the cell (large c_b), so Π4 is large and exceeds 1 in the cold"),
    "earbuds": dict(
        inputs=dict(E_Wh=(0.2, 0.1, 0.3), RQ_sigma=(0.33, 0.2, 0.6), P_typ=(0.025, 0.015, 0.04),
                    P_dem=(0.04, 0.03, 0.06), P_sust=(0.3, 0.2, 0.5), dPb=(0.03, 0.015, 0.06),
                    tau_th_min=(1.0, 0.5, 2.0), CV_P=(0.3, 0.1, 0.5), tau_c_h=(0.05, 0.01, 0.2), mu_h=(6, 4, 8),
                    floor_share=(0.1, 0.05, 0.2)),
        notes="one bud; playback is nearly constant power (small CV_P, small Π5), no thermal limit (Π2 << 1); "
              "Π4 is set by BLE/ANC peaks on a ~50 mAh cell"),
    "e-reader": dict(
        inputs=dict(E_Wh=(5.0, 3.5, 7.5), RQ_sigma=(0.33, 0.2, 0.5), P_typ=(0.1, 0.05, 0.3), P_dem=(0.5, 0.3, 1.0),
                    P_sust=(3.0, 2.0, 4.0), dPb=(1.0, 0.5, 2.0), tau_th_min=(10, 5, 20), CV_P=(2.5, 1.5, 4.0),
                    tau_c_h=(0.5, 0.25, 1.0), mu_h=(1000, 500, 1700), floor_share=(0.7, 0.5, 0.85)),
        notes="TTE ≈ E_u/(P_sleep + ν_page E_page + P_light); weeks of runtime make the TTE very predictable "
              "(small Π5) although single sessions are bursty; sleep floor dominates"),
    "drone": dict(
        inputs=dict(E_Wh=(40, 18, 77), RQ_sigma=(0.05, 0.03, 0.1), P_typ=(100, 50, 180), P_dem=(150, 100, 250),
                    P_sust=(500, 300, 1000), dPb=(100, 50, 200), tau_th_min=(5, 2, 10), CV_P=(0.15, 0.1, 0.3),
                    tau_c_h=(0.01, 0.003, 0.03), mu_h=(0.4, 0.25, 0.6), floor_share=(0.02, 0.01, 0.05)),
        notes="near-constant hover power at 2-3C: the largest Π1; manoeuvre transients give a large Π4 that "
              "exceeds 1 on a cold, aged pack (hence pre-warming advice); battery heats through the whole flight "
              "(Π3 ~ 0.3)"),
}
PI_INPUTS = {"Pi1": ("P_typ", "E_Wh", "RQ_sigma"), "Pi2": ("P_dem", "P_sust"),
             "Pi3": ("tau_th_min", "E_Wh", "P_dem", "P_sust"), "Pi4": ("P_typ", "dPb", "E_Wh", "RQ_sigma"),
             "Pi5": ("tau_c_h", "mu_h", "CV_P")}


# ------------------------------------------------------------------------------------------- RP-63 shapes
P0 = RP63
U_BAR = float(cell.mean_ocv(P0))
RQ_SIGMA_RP = float(cell.R_sigma(Z_NOM, T_REF, P0)) * float(P0.Q_n)  # Ω·Ah
F_R0 = float(P0.R0_ref) / (float(P0.R0_ref) + float(P0.R1_ref))  # R0 share of R_Σ at the reference point
EU_EN = float(cell.usable_energy(P0) / cell.nominal_energy(P0))


def r0_mult(z, T_K, soh):
    """R0(z, T, SOH)/R0_ref of the RP-63 chemistry (A_R low-SOC rise, Arrhenius, 1 + γ_R q)."""
    return float(cell.R0(z, T_K, aged(P0, 1 - soh)) / P0.R0_ref)


def u_ocv(z, T_K):
    return float(cell.ocv(z, T_K, P0))


def pi_values(x, z=Z_NOM, T_K=T_REF, soh=1.0):
    """Π1-Π5 of one input set (dict of scalars), per cell-equivalent."""
    c = x["P_typ"] / x["E_Wh"]
    cb = x["dPb"] / x["E_Wh"]
    rq0 = F_R0 * x["RQ_sigma"] * r0_mult(z, T_K, soh)
    Vs = float(P0.V_sys)
    margin = Vs * (u_ocv(z, T_K) - Vs) / rq0 - c * U_BAR  # (P_cut,sys − P)/Q in W per Ah
    pi4 = cb * U_BAR / margin if margin > 0 else math.inf
    p_eff = min(x["P_dem"], x["P_sust"])
    tte_heavy_h = EU_EN * x["E_Wh"] / p_eff
    return dict(Pi1=c * x["RQ_sigma"] / U_BAR, Pi2=x["P_dem"] / x["P_sust"], Pi3=x["tau_th_min"] / 60.0 / tte_heavy_h,
                Pi4=pi4, Pi5=math.sqrt(2 * x["tau_c_h"] / x["mu_h"]) * x["CV_P"])


def class_pis(inputs):
    mid = {k: v[0] for k, v in inputs.items()}
    out = {}
    for cond, kw in (("", {}), ("_cold_aged", dict(z=Z_COLD, T_K=T_COLD_K, soh=SOH_COLD))):
        m = pi_values(mid, **kw)
        for pi, names in PI_INPUTS.items():
            if cond and pi != "Pi4":
                continue
            vals = []
            for corner in itertools.product(*[(inputs[n][1], inputs[n][2]) for n in names]):
                x = dict(mid)
                x.update(dict(zip(names, corner)))
                vals.append(pi_values(x, **kw)[pi])
            out[pi + cond] = dict(mid=m[pi], lo=min(vals), hi=max(vals))
    return out


# ------------------------------------------------------------------------------------------- main
def main():
    t0 = time.time()
    out = {}
    # ---- RP-63 from the model
    P_SO = float(demand_power(PRIMITIVES["SO"], P0))
    P_GM = float(demand_power(PRIMITIVES["GM"], P0))
    P_sust = float((P0.T_lim_C - P0.dT_thr - 25.0) / (P0.k_env * P0.R_skin))
    tau_th_s = float(P0.C_th * P0.R_th * P0.k_env)
    r_gm = run(Schedule.constant("GM"), 25.0)
    Rs = float(cell.R_sigma(Z_NOM, T_REF, P0))
    pi4 = 3.0 / float(A.safe_burst(Z_NOM, T_REF, P_SO, P0, v1=0.0))
    pa = aged(P0, 1 - SOH_COLD)
    pi4c = 3.0 / float(A.safe_burst(Z_COLD, T_COLD_K, P_SO, pa, v1=0.0))
    pi4c_v1 = 3.0 / float(A.safe_burst(Z_COLD, T_COLD_K, P_SO, pa))
    Rs_c = float(cell.R_sigma(Z_COLD, T_COLD_K, pa))
    ct = {}
    for prof in ("typical", "heavy"):
        c = CTMC(prof)
        ct[prof] = A.ctmc_moments(c.generator(), c.state_powers(P0, 25.0), float(cell.usable_energy(P0)))
    rp = dict(
        Pi1=1.0 * Rs / U_BAR ** 2, Pi1_S4=P_SO * Rs / U_BAR ** 2, Pi1_gaming=r_gm.P_mean_W * Rs / U_BAR ** 2,
        Pi1_cold_aged_S4=P_SO * Rs_c / U_BAR ** 2,
        Pi2=P_GM / P_sust, Pi2_S4=P_SO / P_sust, Pi3=tau_th_s / 3600.0 / r_gm.tte_h, Pi4=pi4, Pi4_cold_aged=pi4c,
        Pi4_cold_aged_with_v1=pi4c_v1, Pi5=ct["typical"]["cv_tte"], Pi5_heavy=ct["heavy"]["cv_tte"],
        inputs=dict(U_bar_V=U_BAR, R_sigma_ohm=Rs, R_sigma_cold_aged_ohm=Rs_c, RQ_sigma_ohm_Ah=RQ_SIGMA_RP,
                    P_S4_W=P_SO, P_gaming_demand_W=P_GM, P_sust_W=P_sust, tau_th_min=tau_th_s / 60,
                    tte_gaming_L2_h=r_gm.tte_h, dPb_W=3.0, E_Wh=float(cell.nominal_energy(P0)),
                    tau_c_h=ct["typical"]["tau_c"], cv_P=ct["typical"]["cv_P"], mu_dynkin_h=ct["typical"]["mu_dynkin"],
                    floor_share=0.045),
        conditions=dict(Pi1="P = 1 W, 25 °C, z = 0.5 (Pi1_S4: social load 1.06 W)",
                        Pi2="gaming demand 3.906 W vs sharp-cap P_sust at 25 °C (Pi2_S4: social)",
                        Pi3="τ_th = k_env R_th C_th vs the L2 gaming TTE at 25 °C",
                        Pi4="S4 (P = 1.06 W, ΔP_b = 3 W) at 25 °C, new, z = 0.5; cold_aged: −20 °C, SOH 0.8, z = 0.2 "
                            "(P_cut,sys without RC polarisation, §7 F4 table)",
                        Pi5="typical-day CTMC (S13), F7"),
        proto=PROTO,
    )
    rp["dev"] = {k: rp[k] - v for k, v in PROTO.items()}
    rp["per_cell_check"] = pi_values(dict(P_typ=P_SO, E_Wh=float(cell.nominal_energy(P0)), RQ_sigma=RQ_SIGMA_RP,
                                          dPb=3.0, P_dem=P_GM, P_sust=P_sust, tau_th_min=tau_th_s / 60,
                                          CV_P=ct["typical"]["cv_P"], tau_c_h=ct["typical"]["tau_c"],
                                          mu_h=ct["typical"]["mu_dynkin"]))
    rp["per_cell_check_cold_aged_Pi4"] = pi_values(
        dict(P_typ=P_SO, E_Wh=float(cell.nominal_energy(P0)), RQ_sigma=RQ_SIGMA_RP, dPb=3.0, P_dem=P_GM,
             P_sust=P_sust, tau_th_min=tau_th_s / 60, CV_P=1, tau_c_h=1, mu_h=1), z=Z_COLD, T_K=T_COLD_K,
        soh=SOH_COLD)["Pi4"]
    out["RP63"] = rp
    print("RP-63: " + ", ".join(f"{k}={rp[k]:.4g}" for k in ("Pi1", "Pi1_S4", "Pi2", "Pi3", "Pi4", "Pi4_cold_aged",
                                                            "Pi5")))
    print(f"  per-cell formula check: Pi1(S4) {rp['per_cell_check']['Pi1']:.5f} vs {rp['Pi1_S4']:.5f}; "
          f"Pi4 {rp['per_cell_check']['Pi4']:.4f} vs {pi4:.4f}; cold {rp['per_cell_check_cold_aged_Pi4']:.4f} "
          f"vs {pi4c:.4f}")

    # ---- device classes
    cls = {}
    for name, d in CLASSES.items():
        pis = class_pis(d["inputs"])
        fs = d["inputs"]["floor_share"]
        cls[name] = dict(**pis, floor_elasticity=dict(mid=-fs[0], lo=-fs[2], hi=-fs[1]),
                         inputs={k: dict(mid=v[0], lo=v[1], hi=v[2]) for k, v in d["inputs"].items()},
                         notes=d["notes"], status="illustrative order of magnitude [unverified]; no device data")
        print(f"  {name:10s} " + ", ".join(f"{k}={v['mid']:.3g} [{v['lo']:.3g}, {v['hi']:.3g}]"
                                           for k, v in pis.items()))
    out["classes"] = cls
    out["scaling_rules"] = {
        "energy": "E_u = N_s Q Ū; TTE = E_u/P̄ still holds",
        "pack_resistance": "R ∝ (N_s/N_p)/Q_cell (A12); RP-63 R0·Q = 0.28 Ω·Ah incl. PCM",
        "thermal": "C_th ∝ mass; R_th, R_skin ∝ 1/(h A_surface); P_sust ∝ (T_lim − T_amb) h A_surface",
        "display": "OLED ∝ A g(b) APL; LCD P_panel + c1 b + c2 b²; e-paper energy per refresh",
        "usage": "device-specific activity set; same PDMP formulas F6-F9",
        "size_invariance": "Π1 = c RQ_Σ/Ū and Π4 = c_b Ū/(V_sys(U − V_sys)/RQ0 − c Ū) depend on C-rates and the "
                           "specific resistance R·Q only, not on pack size or N_s",
    }
    # e-reader TTE example (TTE ≈ E_u/(P_sleep + ν_page E_page + P_light))
    E_u = 5.0 * EU_EN
    P_sleep, nu, E_page, P_light = 0.003, 30 / 86400, 0.6, 0.05 * 0.5 / 24
    tte = E_u / (P_sleep + nu * E_page + P_light)
    out["ereader_example"] = dict(E_u_Wh=E_u, P_sleep_W=P_sleep, pages_per_day=30, E_page_J=E_page,
                                  P_light_avg_W=P_light, tte_h=tte, tte_days=tte / 24,
                                  note="illustrative inputs [unverified]: 3 mW sleep, 30 page turns/day at 0.6 J, "
                                       "front light 50 mW for 30 min/day")
    out = dict(meta=meta("run_generalization.py", seed=None, solver="closed forms; LSODA for the RP-63 gaming TTE",
                         h_s=None, runtime_s=round(time.time() - t0, 1)), **out)
    write_json("generalization.json", out)
    make_figure(json.loads((RESULTS / "generalization.json").read_text()))
    print(f"done in {time.time() - t0:.1f} s")


# ------------------------------------------------------------------------------------------- figure
LABEL_OFF = {  # (dx, dy) in points, ha; leader line drawn when the label is moved away from its point
    "a": {"RP-63": (-12, -12, "right"), "tablet": (14, -16, "left"), "laptop": (-14, 7, "right"),
          "smartwatch": (-14, 3, "right"), "earbuds": (12, -3, "left"), "e-reader": (-6, 3, "right"),
          "drone": (6, -2, "left")},
    "b": {"RP-63": (-4, 16, "right"), "tablet": (-7, 3, "right"), "laptop": (7, -1, "left"),
          "smartwatch": (6, -15, "left"), "earbuds": (-6, -3, "right"), "e-reader": (6, 2, "left"),
          "drone": (6, -3, "left")},
}


def make_figure(o):
    apply()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(W1, 2.6), gridspec_kw=dict(wspace=0.3))
    rp = o["RP63"]
    pts = {"RP-63": dict(Pi1=dict(mid=rp["Pi1_S4"], lo=None, hi=None), Pi4=dict(mid=rp["Pi4"], lo=None, hi=None),
                         Pi4_cold_aged=dict(mid=rp["Pi4_cold_aged"]), Pi2=dict(mid=rp["Pi2"], lo=None, hi=None),
                         Pi5=dict(mid=rp["Pi5"], lo=None, hi=None))}
    pts.update(o["classes"])

    CAP = 40.0  # plotted value for Π4 = inf (JSON null): margin <= 0, the mean load alone breaches V_sys

    def cap(v):
        return CAP if v is None else min(v, CAP)

    def draw(ax, xk, yk, panel):
        for name, e in pts.items():
            is_rp = name == "RP-63"
            colr = BLUE if is_rp else INK2
            x, y = e[xk]["mid"], cap(e[yk]["mid"])
            if panel == "a":
                yc = e["Pi4_cold_aged"]["mid"]
                ax.plot([x, x], [y, cap(yc)], color=RED if is_rp else MUTED, lw=0.8, ls=(0, (2, 1.5)), zorder=2)
                ax.plot([x], [cap(yc)], marker="^" if yc is None else "o", ms=4.2, mfc="white",
                        mec=RED if is_rp else INK2, mew=0.9, zorder=4)
            ax.plot([x], [y], marker="o", ms=5.0 if is_rp else 4.2, color=colr, zorder=5)
            dx, dy, ha = LABEL_OFF[panel][name]
            far = abs(dx) > 8 or abs(dy) > 8
            ax.annotate(name, (x, y), xytext=(dx, dy), textcoords="offset points", fontsize=6.8, ha=ha,
                        va="center", color=INK, fontweight="bold" if is_rp else "normal", zorder=6,
                        arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.5, shrinkA=1, shrinkB=3) if far else None)
        ax.set_xscale("log")
        ax.set_yscale("log")

    # (a) Π1 vs Π4
    draw(ax1, "Pi1", "Pi4", "a")
    ax1.axhspan(1.0, 60, color=NEUTRAL, zorder=0)
    ax1.axhline(1.0, color=RED, lw=0.8, ls=(0, (3, 2)), zorder=1)
    ax1.text(1.1e-3, 30, "brown-out region (Π4 > 1)", fontsize=6.5, color=INK2, va="center")
    xd = pts["drone"]["Pi1"]["mid"]
    ax1.annotate("∞: hover alone\nbreaches V$_{sys}$", (xd, CAP), xytext=(6, -4), textcoords="offset points",
                 fontsize=6.0, color=INK2, va="top", ha="left")
    ax1.set_xlim(1e-3, 0.1)
    ax1.set_ylim(0.01, 60)
    ax1.set_xlabel("Π1 = P R$_\\Sigma$/Ū²  (resistive loss fraction)")
    ax1.set_ylabel("Π4 = ΔP$_b$/(P$_{cut,sys}$ − P)  (burst margin)")
    ax1.set_title("Resistance and bursts", loc="left")
    h1 = [plt.Line2D([], [], color=INK2, marker="o", lw=0, ms=4.2, label="25 °C, new, half charge"),
          plt.Line2D([], [], color=INK2, marker="o", lw=0, ms=4.2, mfc="white", mew=0.9,
                     label="−20 °C, SOH 0.8, 20 % charge")]
    ax1.legend(handles=h1, loc="lower right", fontsize=6.3, handlelength=1.0, borderaxespad=0.2)
    panel_label(ax1, "(a)")

    # (b) Π2 vs Π5
    draw(ax2, "Pi2", "Pi5", "b")
    ax2.axvline(1.0, color=RED, lw=0.8, ls=(0, (3, 2)), zorder=1)
    ax2.axvspan(1.0, 20, color=NEUTRAL, zorder=0)
    ax2.text(1.12, 1.25, "throttling\nbinds (Π2 > 1)", fontsize=6.5, color=INK2, va="top")
    ax2.set_xlim(0.03, 10)
    ax2.set_ylim(0.01, 1.5)
    ax2.set_xlabel("Π2 = P$_{dem}$/P$_{sust}$  (thermal index)")
    ax2.set_ylabel("Π5 = CV of TTE  (usage predictability)")
    ax2.set_title("Heat and usage randomness", loc="left")
    ax2.text(0.035, 0.0115, "classes other than RP-63: illustrative orders of magnitude,\n"
                            "input ranges move each group by up to ~1 decade (JSON)", fontsize=6.0, color=INK2,
             va="bottom")
    panel_label(ax2, "(b)")
    save(fig, "fig_generalization")


if __name__ == "__main__":
    main()
