"""AN0: resolved parameter table (FINAL_SPEC §4, §14).

Outputs
  results/params.json          {symbol: {value, unit, dist, lo, mode, hi, source, conf, ...}} for every Params field
                               (code names), plus "_derived" (E_n, E_u, z_e, U(1), R0/R1/tau1 at T) and
                               "_usage_factors" (§4.6), and "_checks" (code vs a transcription of the spec text)
  paper/tables/tab_params.tex  booktabs table T1 (~38 most important rows; no table wrapper)

Values, units and distributions come from the code (``phonebatt.params``), which is authoritative.  The source
keys and confidence levels are transcribed from the §4 tables.  ``SPEC4`` below is an independent transcription of
the §4 values and UQ ranges; the script compares it with ``UQ_SPEC``/``Params`` and reports every mismatch.
"""
from __future__ import annotations

import math
import re
import time

import numpy as np

from basics_common import bib_keys, f, meta, write_json, write_tex
from phonebatt import RP63, cell
from phonebatt.params import K0, UQ_SPEC, USAGE_UQ_SPEC, param_table

# ------------------------------------------------------------------------------------------------ §4 transcription
# name: (default, dist, lo, hi) with dist in U/T/LU/fixed; SI units as in phonebatt.params (J/mol, W, Ohm, s, ...)
SPEC4 = {
    # 4.1
    "Q_n": (4.70, "T", 4.56, 4.80), "dU": (0.0, "U", -0.020, 0.020), "k_S": (1.0, "U", 0.5, 1.5),
    "R0_ref": (0.060, "LU", 0.035, 0.130), "A_R": (0.80, "U", 0.5, 2.0), "z0": (0.139, "fixed", None, None),
    "Ea_R0": (18e3, "T", 12e3, 27e3), "R1_ref": (0.010, "LU", 0.003, 0.025), "tau1_ref": (30.0, "U", 10.0, 60.0),
    "Ea_R1": (30e3, "T", 25e3, 55e3), "Ea_tau": (21e3, "U", 15e3, 27e3), "V_cut": (3.40, "U", 3.30, 3.50),
    "d_sys": (0.10, "U", 0.0, 0.20), "T_shut_C": (68.0, "fixed", None, None), "eta_c": (0.99, "U", 0.98, 1.00),
    "dV_hv": (0.242, "fixed", None, None), "p_hv": (3.86, "fixed", None, None),
    # 4.2
    "C_th": (150.0, "T", 100.0, 230.0), "R_th": (3.0, "T", 2.0, 4.0), "R_skin": (5.5, "T", 4.0, 9.0),
    "T_lim_C": (45.0, "U", 41.0, 48.0), "dT_thr": (2.0, "fixed", None, None), "w_thr": (1.0, "U", 0.5, 3.0),
    "P_fl": (0.5, "U", 0.3, 1.0), "k_env": (1.0, "U", 0.8, 1.5), "dT_ch": (8.0, "U", 5.0, 12.0),
    "eta_charger": (0.90, "U", 0.85, 0.95),
    # 4.3
    "P_susp": (0.020, "U", 0.010, 0.030), "P_awake": (0.090, "LU", 0.05, 0.25), "E_wake": (0.10, "LU", 0.05, 0.30),
    "P_d0": (0.12, "U", 0.08, 0.25), "p_c": (2.567e-3, "U", 1.23e-3, 5.13e-3), "APL_ref": (0.70, "fixed", None, None),
    "gamma_b": (2.0, "U", 1.0, 2.9), "k_RR": (2.8e-3, "U", 1.0e-3, 3.3e-3), "P_SoC_max": (5.0, "U", 3.0, 7.0),
    "v0": (0.446, "U", 0.35, 0.55), "P_G": (0.27, "T", 0.14, 0.45), "P_cam": (1.15, "U", 0.9, 2.0),
    "P_camf": (0.30, "U", 0.2, 0.5), "P_aud": (0.05, "U", 0.02, 0.15), "P_BT": (0.04, "U", 0.03, 0.05),
    "m_b": (1.0, "U", 0.5, 2.0),
    # 4.4
    "wifi_P_idle": (0.004, "U", 0.002, 0.010), "lte_P_idle": (0.015, "U", 0.005, 0.035),
    "nr_P_idle": (0.030, "U", 0.015, 0.060), "wifi_P_tail": (0.12, "fixed", None, None),
    "lte_P_tail": (0.12, "U", 0.066, 0.178), "nr_P_tail": (0.255, "U", 0.249, 0.593),
    "wifi_T_tail": (0.24, "fixed", None, None), "lte_T_tail": (10.0, "U", 5.0, 11.6),
    "nr_T_tail": (10.3, "U", 10.2, 21.4), "wifi_P_act": (0.13, "U", 0.10, 0.30),
    "lte_P_act": (0.60, "LU", 0.30, 1.30), "nr_P_act": (0.90, "LU", 0.50, 1.80),
    "wifi_T_act": (0.24, "fixed", None, None), "lte_T_act": (0.10, "U", 0.08, 0.20),
    "nr_T_act": (0.10, "U", 0.08, 0.20), "wifi_E_pro": (0.0098, "fixed", None, None),
    "lte_E_pro": (0.31, "U", 0.1, 0.5), "nr_E_pro": (1.0, "U", 0.5, 2.0), "wifi_e_b": (0.015, "LU", 0.005, 0.06),
    "lte_e_b": (0.030, "LU", 0.014, 0.08), "nr_e_b": (0.030, "LU", 0.014, 0.08), "k_ul": (5.0, "U", 2.0, 8.4),
    "wifi_P_rt": (0.30, "U", 0.13, 0.8), "lte_P_rt": (1.20, "U", 0.6, 2.6), "nr_P_rt": (1.40, "U", 1.2, 3.0),
    "gamma_sig": (0.26, "LU", 0.10, 1.0), "RSRP_ref": (-85.0, "fixed", None, None), "m_max": (3.5, "fixed", None, None),
    # 4.5
    "k_c": (1.413e-3, "U", 1.05e-3, 2.09e-3), "E_cal": (58e3, "U", 20e3, 60e3), "g_hi": (1.75, "U", 1.3, 1.9),
    "k_HV": (1.0, "U", 0.0, 1.5), "k_H": (2.697e-3, "fixed", None, None), "k_L": (1.288e-3, "fixed", None, None),
    "w_L": (0.3, "U", 0.1, 0.5), "E_H": (33e3, "U", 25e3, 45e3), "E_L": (55.5e3, "U", 40e3, 70e3),
    "beta_L": (2.64, "U", 0.5, 2.64), "c_ref": (1.0, "fixed", None, None), "gamma_R": (2.5, "LU", 1.5, 7.0),
    "N_80": (1000.0, "fixed", None, None),
}
# spec ranges that are not sampled distributions in the code (discrete sets, "follows calibration", storage only)
SPEC_NOTE = {"k_H": "follows the calibration (§11.1)", "k_L": "follows the calibration (§11.1)",
             "N_80": "{500, 800, 1000} (scenario choice)", "self_discharge_per_month": "0.05-2 %/month (storage only)",
             "q_loss": "state of health input (0 = new)"}

# ------------------------------------------------------------------------------------------------ sources (§4 tables)
# name: (source keys, confidence).  Tags in brackets are FINAL_SPEC provenance tags, not bib keys.
SRC = {
    "Q_n": ("PIX9, MR24", "H"), "b0": ("[RUN] CRM06 fit to AI20", "M"), "b1": ("[RUN] CRM06 fit to AI20", "M"),
    "a0": ("[RUN] CRM06 fit to AI20", "M"), "a1": ("[RUN] CRM06 fit to AI20", "M"),
    "a2": ("[RUN] CRM06 fit to AI20", "M"), "a3": ("[RUN] CRM06 fit to AI20", "M"),
    "dV_hv": ("cell_ecm a.2 construction", "L"), "p_hv": ("cell_ecm a.2 construction", "L"),
    "dU": ("OVEJ19 (mechanism) [ASSUMPTION]", "L"), "k_S": ("[RUN] AI20 entropic fit", "M-L"),
    "R0_ref": ("PINE, BU802a, cell_ecm a.3; plausibility HE19, LEE20 [ASSUMPTION]", "L"), "A_R": ("[RUN] AI20", "M"), "z0": ("[RUN] AI20", "M"),
    "Ea_R0": ("[RUN], BARC22, SUN23, LEE20", "M"), "R1_ref": ("[RUN] AI20, cell_ecm a.3", "M/L"),
    "tau1_ref": ("[RUN]", "M"), "Ea_R1": ("[RUN], BARC22", "M"), "Ea_tau": ("[RUN] MQ19", "M"),
    "V_cut": ("QC-FG, LEE20", "M/H"), "d_sys": ("QC-SMB [ASSUMPTION]", "M/L"), "T_shut_C": ("KANG19", "M"),
    "eta_c": ("PLETT-N", "M"),
    "C_th": ("MAL99, HE19, KANG19", "M"), "R_th": ("HE19 [DERIVED]", "M/L"), "R_skin": ("PU23, KANG19", "M"),
    "T_lim_C": ("CHI15, KANG19, UL-TS", "H"), "dT_thr": ("[ASSUMPTION]", "L"), "w_thr": ("AOSP-TH [ASSUMPTION]", "L"),
    "P_fl": ("[ASSUMPTION]", "L"), "k_env": ("KANG19 [ASSUMPTION]", "L"), "dT_ch": ("KANG19 [ASSUMPTION]", "L"),
    "eta_charger": ("[ASSUMPTION]", "L"),
    "P_susp": ("CH13, CH15, HU12", "H/M"), "P_awake": ("[CAL] BONN26 (V1.1); prior CH10, AOSP", "M/L"),
    "E_wake": ("[DERIVED]", "L"), "P_d0": ("CH13, GU24", "M"), "p_c": ("GU24, CH13", "M"),
    "A_disp": ("definition: 2.455 d^2, d = 6.3 in, 19.5:9", "-"), "APL_ref": ("DZ11", "H"),
    "gamma_b": ("CH13, DH21", "M/L"), "k_RR": ("AT20", "M/L"), "P_SoC_max": ("CH13", "M"),
    "v0": ("LX [DERIVED]", "M"), "P_G": ("KW19, CH13, CH10", "H"), "P_cam": ("CX13 [DERIVED], CH13", "M"),
    "P_camf": ("GU24", "M"), "P_aud": ("CH13, HE20", "M"), "P_BT": ("CH10", "H"),
    "m_b": ("HU12, NA21, LEE20, CH13 (§3.7) [ASSUMPTION]", "L"),
    "wifi_P_idle": ("HU12, CH15", "M"), "lte_P_idle": ("HU12 [DERIVED], CH15", "M"), "nr_P_idle": ("GU24", "M"),
    "wifi_P_tail": ("HU12", "M"), "lte_P_tail": ("NA21", "M"), "nr_P_tail": ("NA21", "M"),
    "wifi_T_tail": ("HU12", "M"), "lte_T_tail": ("NA21", "M"), "nr_T_tail": ("NA21, XU20", "M"),
    "wifi_P_act": ("HU12", "M"), "lte_P_act": ("GU24, HU12 [DERIVED]", "M"), "nr_P_act": ("NA21, GU24, XU20", "M"),
    "wifi_T_act": ("HU12", "M"), "lte_T_act": ("XU20", "M"), "nr_T_act": ("XU20", "M"),
    "wifi_E_pro": ("HU12", "M"), "lte_E_pro": ("HU12", "M"), "nr_E_pro": ("NA21 [DERIVED]", "M"),
    "wifi_e_b": ("SU14, NA21, HU12", "M"), "lte_e_b": ("SU14, NA21, HU12", "M"), "nr_e_b": ("SU14, NA21, HU12", "M"),
    "k_ul": ("NA21, HU12", "M"), "wifi_P_rt": ("SU14", "M"), "lte_P_rt": ("HU12, GU24", "M"),
    "nr_P_rt": ("XU20, GU24", "M"), "gamma_sig": ("DI13 [DERIVED]", "M"), "RSRP_ref": ("definition", "-"),
    "m_max": ("DI13 (Tx/Rx current ratios)", "M"),
    "k_c": ("SCH14, KEIL16", "M"), "E_cal": ("SCH14 via CA21, SCM18", "M"), "g_hi": ("KEIL16", "M"),
    "k_HV": ("SCH14, NEL18", "L"), "k_H": ("calibrated (§11.1): AP-BATT, GOOG", "M/L"),
    "k_L": ("calibrated (§11.1): AP-BATT, GOOG", "M/L"), "w_L": ("WAL14, WU17 [ASSUMPTION]", "L"),
    "E_H": ("SCM18", "M"), "E_L": ("SCM18", "M"), "beta_L": ("SCM18", "M/L"), "c_ref": ("SCM18 (definition)", "-"),
    "gamma_R": ("SCH14, KEIL16", "L"), "N_80": ("AP-BATT, GOOG", "H"), "self_discharge_per_month": ("ZIL19", "M/L"),
    "q_loss": ("state input", "-"),
}
SECTION = {}
for _n in ("Q_n b0 b1 a0 a1 a2 a3 dV_hv p_hv dU k_S R0_ref A_R z0 Ea_R0 R1_ref tau1_ref Ea_R1 Ea_tau V_cut d_sys "
           "T_shut_C eta_c").split():
    SECTION[_n] = "4.1"
for _n in "C_th R_th R_skin T_lim_C dT_thr w_thr P_fl k_env dT_ch eta_charger".split():
    SECTION[_n] = "4.2"
for _n in ("P_susp P_awake E_wake P_d0 p_c A_disp APL_ref gamma_b k_RR P_SoC_max v0 P_G P_cam P_camf P_aud P_BT "
           "m_b").split():
    SECTION[_n] = "4.3"
for _n in "k_c E_cal g_hi k_HV k_H k_L w_L E_H E_L beta_L c_ref gamma_R N_80 self_discharge_per_month".split():
    SECTION[_n] = "4.5"

# ------------------------------------------------------------------------------------------------ table T1 rows
# (name, LaTeX symbol, description, display scale, display unit, decimals)
GROUPS = [
    ("Cell and pack (\\S4.1)", [
        ("Q_n", r"$Q_n$", "nominal capacity", 1, "Ah", 2),
        ("R0_ref", r"$R_{0,\mathrm{ref}}$", "pack $R_0$ at $z=0.5$, 25\\,\\textdegree C", 1e3, "m$\\Omega$", 0),
        ("A_R", r"$A_R$", "low-SOC rise of $R_0$", 1, "--", 2),
        ("Ea_R0", r"$E_{a,R_0}$", "activation energy of $R_0$", 1e-3, "kJ/mol", 0),
        ("R1_ref", r"$R_{1,\mathrm{ref}}$", "RC-branch resistance", 1e3, "m$\\Omega$", 0),
        ("tau1_ref", r"$\tau_{1,\mathrm{ref}}$", "RC time constant", 1, "s", 0),
        ("Ea_R1", r"$E_{a,R_1}$", "activation energy of $R_1$", 1e-3, "kJ/mol", 0),
        ("Ea_tau", r"$E_{a,\tau}$", "activation energy of $\\tau_1$", 1e-3, "kJ/mol", 0),
        ("V_cut", r"$V_{\mathrm{cut}}$", "fuel-gauge cut-off voltage", 1, "V", 2),
        ("d_sys", r"$\Delta_{\mathrm{sys}}$", "brown-out level below $V_{\\mathrm{cut}}$", 1, "V", 2),
        ("dU", r"$\delta U$", "OCV offset (hysteresis)", 1e3, "mV", 0),
        ("k_S", r"$k_S$", "multiplier on $\\partial U/\\partial T$", 1, "--", 1),
    ]),
    ("Thermal (\\S4.2)", [
        ("C_th", r"$C_{\mathrm{th}}$", "heat capacity of the phone", 1, "J/K", 0),
        ("R_th", r"$R_{\mathrm{th}}$", "thermal resistance to ambient", 1, "K/W", 1),
        ("R_skin", r"$R_{\mathrm{skin}}$", "skin hot-spot resistance", 1, "K/W", 1),
        ("T_lim_C", r"$T_{\mathrm{lim}}$", "skin temperature limit", 1, "\\textdegree C", 0),
        ("w_thr", r"$w_{\mathrm{thr}}$", "throttle logistic width", 1, "K", 1),
        ("P_fl", r"$P_{\mathrm{fl}}$", "SoC floor when throttled", 1, "W", 1),
        ("k_env", r"$k_{\mathrm{env}}$", "enclosure factor (case, pocket)", 1, "--", 1),
    ]),
    ("Device power, battery side (\\S4.3)", [
        ("P_susp", r"$P_{\mathrm{susp}}$", "suspend floor", 1e3, "mW", 0),
        ("P_awake", r"$\Delta P_{\mathrm{awake}}$", "awake overhead (calibrated)", 1e3, "mW", 0),
        ("P_d0", r"$P_{d0}$", "OLED fixed overhead", 1, "W", 2),
        ("p_c", r"$p_c$", "OLED content power per area", 1e3, "mW/cm$^2$", 2),
        ("gamma_b", r"$\gamma_b$", "brightness exponent", 1, "--", 1),
        ("k_RR", r"$\kappa_{RR}$", "refresh-rate cost above 60\\,Hz", 1e3, "mW/Hz", 1),
        ("P_SoC_max", r"$P_{\mathrm{SoC,max}}$", "maximum CPU+GPU power", 1, "W", 1),
        ("v0", r"$v_0$", "DVFS voltage intercept", 1, "--", 3),
        ("P_G", r"$P_G$", "GNSS receiver", 1, "W", 2),
        ("P_cam", r"$\Delta P_{\mathrm{cam}}$", "rear camera pipeline", 1, "W", 2),
        ("m_b", r"$m_b$", "burst-envelope multiplier", 1, "--", 1),
    ]),
    ("Radio, LTE (\\S4.4)", [
        ("lte_P_act", r"$P_{\mathrm{act}}$", "connected-active power", 1, "W", 2),
        ("lte_e_b", r"$e_b$", "energy per bit", 1e3, "mW/(Mbit/s)", 0),
        ("lte_P_rt", r"$P_{\mathrm{rt}}$", "real-time (call) power", 1, "W", 2),
        ("gamma_sig", r"$\gamma_{\mathrm{sig}}$", "weak-signal slope (per 10\\,dB)", 1, "--", 2),
    ]),
    ("Aging (\\S4.5)", [
        ("k_c", r"$k_c$", "calendar-fade rate", 1e3, "$10^{-3}\\,$d$^{-1/2}$", 3),
        ("E_cal", r"$E_{\mathrm{cal}}$", "calendar activation energy", 1e-3, "kJ/mol", 0),
        ("k_H", r"$k_H$ / $k_L$", "SEI / plating cycle-fade rates", 1e3, "$10^{-3}\\,$FCE$^{-1/2}$", 3),
        ("gamma_R", r"$\gamma_R$", "resistance growth per unit fade", 1, "--", 1),
    ]),
]
TAG_TEX = {"[RUN]": "fit", "[DERIVED]": "der.", "[ASSUMPTION]": "assum.", "[CAL]": "\\textbf{cal.}"}


def source_tex(src: str, keys: set[str]) -> str:
    """Source string -> LaTeX: bib keys become one \\cite{}, provenance tags short words, other notes dropped."""
    tags = [TAG_TEX[t] for t in TAG_TEX if t in src]
    toks = re.findall(r"[A-Za-z][A-Za-z0-9\-]*", re.sub(r"\[[A-Z]+\]", " ", src))
    cites = []
    for t in toks:
        if t in keys and t not in cites:
            cites.append(t)
    parts = []
    if "[CAL]" in src:
        parts.append(TAG_TEX["[CAL]"])
        tags = [t for t in tags if t != TAG_TEX["[CAL]"]]
    parts += tags
    if cites:
        parts.append("\\cite{" + ",".join(cites) + "}")
    return " ".join(parts) if parts else "--"


def fmt(x: float, nd: int) -> str:
    s = f"{x:.{nd}f}"
    return s.replace("-", "$-$") if s.startswith("-") else s


def uq_tex(name: str, scale: float, nd: int) -> str:
    e = UQ_SPEC.get(name)
    if e is None:
        return "calibrated" if name in ("k_H", "k_L") else "fixed"
    dist, lo, _, hi = e
    nd2 = nd if nd > 0 else 0
    lo_s, hi_s = fmt(lo * scale, nd2), fmt(hi * scale, nd2)
    if dist == "T":
        return f"T({lo_s}, {hi_s})"
    return f"{dist}({lo_s}, {hi_s})"


def main() -> None:
    t0 = time.time()
    keys = bib_keys()
    tab = param_table(RP63)
    out = {"meta": meta("run_params.py", seed=None, solver=None, h_s=None,
                        note="value/unit/dist/lo/mode/hi from phonebatt.params (code); source/conf transcribed "
                             "from FINAL_SPEC §4. Units SI as in the code (J/mol, Ohm, W, W/cm2, s, K).")}
    for name, d in tab.items():
        src, conf = SRC.get(name, ("", ""))
        e = dict(d)
        e["value"] = float(np.asarray(e["value"])) if not isinstance(e["value"], str) else e["value"]
        e.update(source=src, conf=conf, section=SECTION.get(name, "4.4" if name.split("_")[0] in ("wifi", "lte", "nr")
                                                            or name in ("k_ul", "gamma_sig", "RSRP_ref", "m_max")
                                                            else "2.10"))
        if name in SPEC_NOTE:
            e["spec_range"] = SPEC_NOTE[name]
        out[name] = e

    # ---- cross-check: code vs the transcription of §4 ------------------------------------------------------
    mism = []
    for name, (val, dist, lo, hi) in SPEC4.items():
        code_val = float(np.asarray(getattr(RP63, name)))
        if not math.isclose(code_val, val, rel_tol=2e-3, abs_tol=1e-12):
            mism.append(dict(name=name, what="default", spec=val, code=code_val))
        e = UQ_SPEC.get(name)
        if dist == "fixed":
            if e is not None:
                mism.append(dict(name=name, what="dist", spec="fixed", code=list(e)))
            continue
        if e is None:
            mism.append(dict(name=name, what="dist", spec=[dist, lo, hi], code="fixed"))
            continue
        if e[0] != dist or not math.isclose(e[1], lo, rel_tol=1e-9, abs_tol=1e-12) or \
                not math.isclose(e[3], hi, rel_tol=1e-9, abs_tol=1e-12) or \
                not math.isclose(e[2], code_val, rel_tol=1e-9, abs_tol=1e-12):
            mism.append(dict(name=name, what="UQ", spec=[dist, lo, val, hi], code=list(e)))
    missing_src = [n for n in tab if n not in SRC]
    unspecified = [n for n in UQ_SPEC if n not in SPEC4]
    out["_checks"] = dict(n_compared=len(SPEC4), mismatches=mism, fields_without_source=missing_src,
                          uq_entries_not_in_transcription=unspecified)

    # ---- derived quantities (§2.5, §2.6, §4) --------------------------------------------------------------
    der = dict(z_e=f(cell.z_empty(RP63)), U1_V=f(cell.ocv(1.0, 298.15, RP63)), mean_ocv_V=f(cell.mean_ocv(RP63)),
               E_n_Wh=f(cell.nominal_energy(RP63)), E_u_Wh=f(cell.usable_energy(RP63)),
               V_sys_V=float(RP63.V_sys), A_disp_cm2=float(RP63.A_disp), p_c_A_W=float(RP63.p_c * RP63.A_disp))
    der["R_vs_T"] = {f"{Tc}": dict(R0_mOhm=1e3 * f(cell.R0(0.5, Tc + K0, RP63)), R1_mOhm=1e3 * f(cell.R1(Tc + K0, RP63)),
                                   tau1_s=f(cell.tau1(Tc + K0, RP63))) for Tc in (40, 25, 10, 0, -10, -20)}
    der["proto"] = dict(z_e=0.0263, U1_V=4.4496, mean_ocv_V=3.8699, E_n_Wh=18.189, E_u_Wh=17.788)
    out["_derived"] = der
    out["_usage_factors"] = {k: dict(dist=v[0], lo=v[1], mode=v[2], hi=v[3]) if v[0] != "DIR" else
                             dict(dist="Dirichlet", concentration=v[1], base=list(v[2]))
                             for k, v in USAGE_UQ_SPEC.items()}
    write_json("params.json", out)

    # ---- LaTeX table T1 ---------------------------------------------------------------------------------
    lines = [
        "% T1: resolved RP-63 parameters (generated by analysis/run_params.py; do not edit by hand).",
        "% UQ: U uniform, LU log-uniform, T triangular with mode = value. Source tags: fit = our fit to open data or",
        "% PyBaMM; der. = derived from cited numbers; assum. = engineering assumption; cal. = the single calibrated",
        "% constant (V1.1). Full table with all 90+ fields: results/params.json.",
        "\\small",
        "\\setlength{\\tabcolsep}{3.5pt}",
        "\\begin{tabular}{@{}llrlll@{}}",
        "\\toprule",
        "Symbol & Meaning & Value & Unit & UQ range & Source \\\\",
        "\\midrule",
    ]
    n_rows = 0
    used_cites = set()
    for gi, (gname, rows) in enumerate(GROUPS):
        if gi:
            lines.append("\\addlinespace[2pt]")
        lines.append(f"\\multicolumn{{6}}{{@{{}}l}}{{\\textit{{{gname}}}}} \\\\")
        for name, sym, desc, scale, unit, nd in rows:
            v = float(np.asarray(getattr(RP63, name))) * scale
            vs = fmt(v, nd)
            if name == "k_H":  # one row for the two calibrated cycle-fade rates
                vs = f"{vs} / {fmt(float(RP63.k_L) * scale, nd)}"
            src = source_tex(SRC[name][0], keys)
            used_cites |= set(re.findall(r"\\cite\{([^}]*)\}", src))
            lines.append(f"{sym} & {desc} & {vs} & {unit} & {uq_tex(name, scale, nd)} & {src} \\\\")
            n_rows += 1
    lines += ["\\bottomrule", "\\end{tabular}", ""]
    write_tex("tab_params.tex", "\n".join(lines))
    cited = {k for c in used_cites for k in c.split(",")}
    bad = sorted(cited - keys)
    print(f"params.json: {len(tab)} fields; spec-vs-code mismatches: {len(mism)}; fields without source: "
          f"{missing_src}; UQ entries not transcribed: {unspecified}")
    for m in mism:
        print("  MISMATCH", m)
    print(f"tab_params.tex: {n_rows} rows; unresolved cite keys: {bad}")
    print(f"derived: z_e={der['z_e']:.4f} U(1)={der['U1_V']:.4f} meanOCV={der['mean_ocv_V']:.4f} "
          f"E_n={der['E_n_Wh']:.3f} E_u={der['E_u_Wh']:.3f}  ({time.time() - t0:.1f} s)")


if __name__ == "__main__":
    main()
