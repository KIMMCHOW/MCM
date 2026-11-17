"""AN6: pre-registered validation V1–V7 (FINAL_SPEC §8, §14).

Outputs
-------
* ``results/validation.json``: ``V1.{dev}|{case}.{model_h, obs_h, err_pct, P_W, ...}``, ``V1.additivity_residual_pct``,
  ``V1.slope_test{D1, D2}``, ``V2.{T}|{rate}.{fQ_model, fQ_obs, rmse_V, ...}``, ``V3``, ``V4.rmse_V``,
  ``V5.{case}.dtte_pct``, ``V6.{case}.{model_h, obs_h, err_pct, pass}``, ``V7`` and ``criteria`` (pass/fail of
  every pre-registered criterion).
* ``paper/figures/fig_validation_{bonn,hg2,band,pybamm}.pdf``; ``paper/tables/tab_validation.tex``;
  ``results/validation_summary.md``.

Rules (§8): the data test predictions; the only tuned constants are the declared [CAL] steps (V1.1 ΔP_awake;
V2 HG2 OCV/thermal at 25 °C; V4 Oxford OCV from low-current samples).  Everything else is the §4 default.

Usage: ``cd /home/user/MCM && PYTHONPATH=src:analysis python3 analysis/run_validation.py [--n-uq 512]
[--n-prior 4096] [--skip-pybamm]``.
"""
from __future__ import annotations

import argparse
import math
import time
import warnings
from dataclasses import replace

import numpy as np

import phonebatt
from phonebatt import UQ_SPEC, cell, run, run_batch, sample
from phonebatt import validation as V
from phonebatt.params import RP63, K0, device, iphone13pro, iphone14pro, nasa_b0005
from phonebatt.power import demand_power, steady_power
from phonebatt.usage import BONN_CASES, PRIMITIVES, Schedule, UsageState

from basics_common import meta, write_json, write_tex
from validation_cells import CellModel, rmse, run_to_event

SEED = 2026
DEVS = {"14 Pro": iphone14pro, "13 Pro": iphone13pro}
CASES = ("idle", "cam", "gps", "cxl")
CASE_LABEL = {"idle": "idle", "cam": "camera", "gps": "GPS", "cxl": "camera+GPS"}

# [PROTO] regression targets quoted in FINAL_SPEC §8 (not validation results)
PROTO_V1 = {("14 Pro", "idle"): 54.48, ("14 Pro", "cam"): 8.96, ("14 Pro", "gps"): 24.98, ("14 Pro", "cxl"): 7.48,
            ("13 Pro", "idle"): 51.85, ("13 Pro", "cam"): 8.52, ("13 Pro", "gps"): 23.77, ("13 Pro", "cxl"): 7.11}
PROTO_V1_BLIND_ERR = {("14 Pro", "idle"): -5.1, ("14 Pro", "cam"): 5.6, ("14 Pro", "gps"): -49.5,
                      ("14 Pro", "cxl"): -11.3, ("13 Pro", "idle"): -2.2, ("13 Pro", "cam"): -11.2,
                      ("13 Pro", "gps"): -21.7, ("13 Pro", "cxl"): -15.2}
PROTO_V6 = {"14P_local_video": 26.08, "14P_stream_video": 23.40, "14P_audio_BT": 55.90, "14PM_wifi_web_150": 15.18,
            "14PM_local_video_150": 33.01, "14PM_idle_minbright": 66.63, "14P_cell_web_150": 10.29,
            "Pixel8_wifi_web": 16.37, "Pixel8_local_video": 35.70, "Pixel8_idle_minbright": 70.52,
            "S24_wifi_web": 14.30, "S24_local_video": 31.20, "S24_idle_minbright": 61.64, "Pixel8a_S12_day": 34.67}
PROTO_SLOPE_OCV = (4.059, 3.756, 3.645)
WINDOWS = ((99, 50), (50, 20), (20, 1))
H_TOP = 0.065  # iOS top hold (§2.11)


def _t():
    return time.perf_counter()


# =========================================================================================== V1 Bonn
def _device_theta(th: dict, Q_dev: float) -> dict:
    """Map a §4 sample to a validation device (§4.7): relative capacity uncertainty and R·Q = const."""
    out = {k: np.array(v, float) for k, v in th.items() if np.ndim(v) == 1}
    Qd = Q_dev * out["Q_n"] / 4.70
    out["R0_ref"] = out["R0_ref"] * 4.70 / Qd
    out["R1_ref"] = out["R1_ref"] * 4.70 / Qd
    out["Q_n"] = Qd
    return out


def bonn_uq(th_base: dict, label: str) -> dict:
    """TTE quantiles of the eight Bonn cases for one parameter sample set (run_batch, RK4)."""
    out = {}
    for dev, ctor in DEVS.items():
        p = ctor()
        th = _device_theta(th_base, float(p.Q_n))
        for case in CASES:
            # h = 10 s (§6.3).  h = 30 s is NOT allowed here: τ1,ref ~ U(10, 60) s puts |hλ| above RK4's 2.785
            # stability limit for the fastest samples (found while testing: one sample ended after 0.14 h).
            with warnings.catch_warnings():  # benign overflow in the logistic of RK4 stages; see spot check
                warnings.simplefilter("ignore", RuntimeWarning)
                b = run_batch(Schedule.constant(BONN_CASES[case]), 25.0, th, h=UQ_H, p=p)
            out[(dev, case)] = b["tte_h"]
            out[(dev, case, "event")] = b["event"]
        print(f"  UQ {label}: {dev} done", flush=True)
    return out


UQ_H = 10.0


def uq_spot_check(th_base: dict, uq: dict, idx) -> dict:
    """LSODA reference runs for a few UQ samples: max relative TTE difference and event mismatches."""
    errs, mism = [], 0
    for dev, ctor in DEVS.items():
        p = ctor()
        th = _device_theta(th_base, float(p.Q_n))
        for case in CASES:
            for k in idx:
                pk = p.replace(**{n: float(v[k]) for n, v in th.items()})
                r = run(Schedule.constant(BONN_CASES[case]), 25.0, p=pk)
                errs.append(abs(uq[(dev, case)][k] / r.tte_h - 1))
                mism += int(uq[(dev, case, "event")][k] != r.event)
    return dict(n_runs=len(errs), max_rel_err=float(max(errs)), event_mismatches=mism, samples=[int(k) for k in idx])


def _display_times(t_h, frac, levels):
    """First time the iOS top-hold display (floor_1%(min(1, frac/(1-h)))) shows each level L."""
    thr = (np.asarray(levels, float) + 1.0) / 100.0 * (1.0 - H_TOP)
    # frac is decreasing in time: invert by interpolation
    return np.interp(-thr, -frac, t_h)


def _window_slopes(times: dict) -> np.ndarray:
    return np.array([(times[b] - times[a]) * 60.0 / (a - b) for a, b in WINDOWS])  # min per %


def v1_bonn(n_uq: int, n_prior: int) -> tuple[dict, dict]:
    t0 = _t()
    obs = V.bonn_observations()
    pred = V.bonn_predictions()  # P_awake = 0.090 (§4.3 default = rounded CAL)
    pred_blind = V.bonn_predictions(P_awake=0.10)
    aw_cal = V.bonn_calibrate_awake()
    pred_cal = V.bonn_predictions(P_awake=aw_cal)
    # RK4 step check for the UQ runs (nominal parameters)
    step_chk = {}
    for case, h in (("idle", UQ_H), ("cam", UQ_H)):
        b = run_batch(Schedule.constant(BONN_CASES[case]), 25.0, None, h=h, p=iphone14pro(), n=1)
        step_chk[case] = dict(h_s=h, rel_err_pct=100 * (b["tte_h"][0] / pred[("14 Pro", case)]["model_h"] - 1))

    # ---- UQ: (A) P_awake fixed at the [CAL] default, all other §4 factors sampled (primary)
    th_cal = sample(UQ_SPEC, n_uq, seed=SEED)
    th_cal["P_awake"] = np.full(n_uq, float(RP63.P_awake))
    uq_cal = bonn_uq(th_cal, f"CAL-fixed N={n_uq}")
    # ---- (B) blind prior: P_awake from its LU(0.05, 0.25) prior; (C) ABC: (B) conditioned on the 14 Pro idle
    th_pr = sample(UQ_SPEC, n_prior, seed=SEED + 1)
    uq_pr = bonn_uq(th_pr, f"prior N={n_prior}")
    o_idle = obs[("14 Pro", "idle")]["obs_h"]
    acc = np.abs(uq_pr[("14 Pro", "idle")] / o_idle - 1) <= 0.05
    spot = uq_spot_check(th_cal, uq_cal, np.linspace(0, n_uq - 1, 8).astype(int))
    print(f"  UQ spot check vs LSODA: {spot}", flush=True)
    q = lambda a: [float(x) for x in np.percentile(a[np.isfinite(a)], [5, 50, 95])]

    V1 = {}
    for (dev, case), r in pred.items():
        o = obs[(dev, case)]
        qa, qp, qb = q(uq_cal[(dev, case)]), q(uq_pr[(dev, case)]), q(uq_pr[(dev, case)][acc])
        hold = not (dev == "14 Pro" and case == "idle")
        cover = lambda qq: bool(qq[0] <= o["min_h"] and o["max_h"] <= qq[2])
        V1[f"{dev}|{case}"] = dict(
            model_h=r["model_h"], obs_h=r["obs_h"], err_pct=r["err_pct"], P_W=r["P_W"], event=r["event"],
            obs_min_h=o["min_h"], obs_max_h=o["max_h"], n_devices=o["n"], experiment=o["experiment"],
            holdout=hold, err_blind_pct=pred_blind[(dev, case)]["err_pct"],
            model_cal_exact_h=pred_cal[(dev, case)]["model_h"], err_cal_exact_pct=pred_cal[(dev, case)]["err_pct"],
            uq_cal_fixed=dict(q05_h=qa[0], q50_h=qa[1], q95_h=qa[2], covers_obs_range=cover(qa)),
            uq_prior=dict(q05_h=qp[0], q50_h=qp[1], q95_h=qp[2], covers_obs_range=cover(qp)),
            uq_abc=dict(q05_h=qb[0], q50_h=qb[1], q95_h=qb[2], covers_obs_range=cover(qb)),
            pass_a=(abs(r["err_pct"]) <= 15.0) if hold else None, pass_b=cover(qa),
            proto_model_h=PROTO_V1[(dev, case)], dev_vs_proto_pct=100 * (r["model_h"] / PROTO_V1[(dev, case)] - 1),
            proto_err_blind_pct=PROTO_V1_BLIND_ERR[(dev, case)])

    # ---- additivity (criterion c): model and data
    groups = V.load_bonn_summary()["groups"]

    def data_add(dev, cam_exp):
        g = lambda e: groups[f"{e}|{dev} (inferred)"]["runtime_h_median"]
        return 100 * ((1 / g("23-11-2023-CXL")) / (1 / g(cam_exp) + 1 / g("02-10-2023-GPS") - 1 / g("07-09-2023-Idle")) - 1)

    add = {dev: dict(model=100 * V.bonn_additivity(pred, dev), data_cam30=data_add(dev, "09-11-2023-CAM"),
                     data_cam60=data_add(dev, "14-11-2023-CAM"), pass_c=abs(100 * V.bonn_additivity(pred, dev)) <= 1.0)
           for dev in DEVS}
    add["proto_model_14Pro_pct"] = 0.41
    add["data_reported_pct"] = {"14 Pro": 0.15, "13 Pro": -0.10}

    # ---- V1.4 SOC-shape (display-map) test
    curves = V.bonn_soc_curves()
    lv = curves["level_pct"].to_numpy()
    model_tr = {}
    for dev, ctor in DEVS.items():
        p = ctor()
        for case in CASES:
            r = run(Schedule.constant(BONN_CASES[case]), 25.0, p=p, record=True)
            tr = r.traj
            s1 = np.clip(tr["s"], -1, 2)
            s2 = cell.display_soc(tr["z"], p, "D2")
            lv_need = sorted({a for w in WINDOWS for a in w})
            out = {}
            for name, frac in (("D1", s1), ("D2", s2)):
                tt = _display_times(tr["t_h"], frac, lv_need)
                out[name] = dict(zip(lv_need, tt))
            model_tr[(dev, case)] = out
    kind_of = {"Idle": "idle", "GPS": "gps", "CAM": "cam", "CXL": "cxl"}
    traces = {}
    errs = {"D1": [], "D2": []}
    for col in curves.columns[1:]:
        exp, devtag, _ = col.rsplit("_", 2)
        dev = devtag.replace("Pro", " Pro")
        kind = kind_of[exp.split("-")[-1]]
        tcol = curves[col].to_numpy()
        tobs = {L: float(tcol[lv == L][0]) for L in (99, 50, 20, 1)}
        m_obs = _window_slopes(tobs)
        r_obs = m_obs / ((tobs[1] - tobs[99]) * 60 / 98)
        rec = dict(device=dev, kind=kind, slopes_obs_min_per_pct=m_obs.tolist(), ratio_obs=r_obs.tolist())
        for name in ("D1", "D2"):
            tm = model_tr[(dev, kind)][name]
            m_mod = _window_slopes(tm)
            r_mod = m_mod / ((tm[1] - tm[99]) * 60 / 98)
            e = 100 * np.sqrt(np.mean((r_mod - r_obs) ** 2))
            rec[name] = dict(slopes_model_min_per_pct=m_mod.tolist(), ratio_model=r_mod.tolist(), rms_err_pct=e,
                             max_abs_err_pct=float(100 * np.max(np.abs(r_mod - r_obs))))
            errs[name].append(e)
        traces[col] = rec
    slope = {}
    for name in ("D1", "D2"):
        e = np.array(errs[name])
        slope[name] = dict(rms_err_pct_mean=float(e.mean()), rms_err_pct_median=float(np.median(e)),
                           rms_err_pct_max=float(e.max()), n_traces=int(e.size),
                           ratio_model_14Pro_idle=traces["07-09-2023-Idle_14Pro_h"][name]["ratio_model"])
    n_d2_better = int(sum(traces[c]["D2"]["rms_err_pct"] < traces[c]["D1"]["rms_err_pct"] for c in traces))
    slope["n_traces_D2_better"] = n_d2_better
    slope["favours"] = "D2" if slope["D2"]["rms_err_pct_mean"] < slope["D1"]["rms_err_pct_mean"] else "D1"
    # [PROTO] check: mean OCV in the D1 display windows (14 Pro, constant P)
    p14 = iphone14pro()
    ze = float(cell.z_empty(p14))
    wocv = []
    for lo, hi in ((0.99, 0.50), (0.50, 0.20), (0.20, 0.01)):
        zl, zh = ze + lo * (1 - ze), ze + hi * (1 - ze)
        wocv.append(float((cell.ocv_integral(zl, p14) - cell.ocv_integral(zh, p14)) / (zl - zh)))
    slope["proto_check"] = dict(window_mean_ocv_V=wocv, proto_V=list(PROTO_SLOPE_OCV),
                                ratio=[w / wocv[0] for w in wocv],
                                max_dev_mV=float(1000 * max(abs(a - b) for a, b in zip(wocv, PROTO_SLOPE_OCV))))
    m_idle = traces["07-09-2023-Idle_14Pro_h"]["slopes_obs_min_per_pct"]
    slope["obs_14Pro_idle_min_per_pct"] = m_idle
    slope["traces"] = traces

    V1["additivity_residual_pct"] = add
    V1["slope_test"] = slope
    V1["calibration"] = dict(P_awake_exact_W=aw_cal, P_awake_default_W=float(RP63.P_awake), proto_exact_W=0.0877,
                             P_awake_blind_prior_W=0.10, target_h=o_idle)
    V1["uq_settings"] = dict(
        primary="uq_cal_fixed: P_awake held at the [CAL] default 0.090 W; every other §4 factor (UQ_SPEC) sampled "
                "with scrambled Sobol' (seed 2026); device capacity Q_dev*Q_n/4.70 and R·Q = const (A12)",
        n_cal_fixed=n_uq, n_prior=n_prior, prior_seed=SEED + 1,
        abc="uq_abc: prior samples whose 14 Pro idle TTE is within ±5 % of the observed 55.03 h (§9 item 4)",
        n_abc_accepted=int(acc.sum()), solver="phonebatt.run_batch (RK4)", h_s=UQ_H,
        step_check_vs_lsoda=step_chk, lsoda_spot_check=spot,
        non_finite=int(sum((~np.isfinite(v)).sum() for k, v in uq_cal.items() if len(k) == 2)))
    V1["_runtime_s"] = _t() - t0
    fig = dict(pred=pred, obs=obs, uq_cal=uq_cal, uq_pr=uq_pr, acc=acc, traces=traces, model_tr=model_tr,
               curves=curves)
    return V1, fig


# =========================================================================================== V2 LG HG2
RATES = {"C/20": "C20DisCh", "0.5C": "Dis_0p5C", "2C": "Dis_2C"}
HG2_T = (40, 25, 10, 0, -10, -20)
HG2_OBS_FQ = {"0.5C": {-20: 0.73, -10: 0.86, 0: 0.91, 10: 0.93, 40: 0.97},
              "2C": {-20: 0.49, -10: 0.73, 0: 0.88, 10: 0.93, 40: 0.98}}  # §8 V2 targets (validation_data a.1)


def _hg2_run(p, variant, d, T, **kw):
    I = float(d["I_A"].mean())
    m = run_to_event(CellModel(p, variant, **kw), "I", I, float(T), t_max_s=1.5 * float(d["t_s"][-1]) + 7200)
    Vm = m.V_at(d["t_s"])
    return dict(Ah=m.Ah, Wh=m.Wh, V_mean=m.V_mean, tte_s=m.tte_s, event=m.event, T_end_C=float(m.T_C[-1]),
                rmse_V=rmse(Vm, d["V"]), n_rmse=int(np.isfinite(Vm).sum()), sim=m, I=I)


def v2_hg2() -> tuple[dict, dict]:
    t0 = _t()
    cal = V.hg2_calibration()
    p = V.hg2_params(cal)
    data = {(T, r): V.load_hg2(T, k) for T in HG2_T for r, k in RATES.items()}
    variants = {"L2": dict(variant="L2"), "X9": dict(variant="X9"), "X9_noR1": dict(variant="X9", x9_r1=False)}
    res = {}
    for vn, kw in variants.items():
        for (T, r), d in data.items():
            vv = kw["variant"]
            extra = {k: v for k, v in kw.items() if k != "variant"}
            res[(vn, T, r)] = _hg2_run(p, vv, d, T, **extra)
    out = {}
    for (T, r), d in data.items():
        Qo, Q25o = float(d["Ah"][-1]), float(data[(25, r)]["Ah"][-1])
        dt = np.diff(np.concatenate([[0.0], d["t_s"]]))
        Eo = float(np.sum(d["I_A"] * d["V"] * dt) / 3600)
        dt25 = np.diff(np.concatenate([[0.0], data[(25, r)]["t_s"]]))
        E25o = float(np.sum(data[(25, r)]["I_A"] * data[(25, r)]["V"] * dt25) / 3600)
        fQo = Qo / Q25o
        c20 = float(data[(T, "C/20")]["Ah"][-1]) / float(data[(25, "C/20")]["Ah"][-1])
        rec = dict(fQ_obs=fQo, fE_obs=Eo / E25o, Q_obs_Ah=Qo, V_mean_obs=Eo / Qo, I_A=float(d["I_A"].mean()),
                   T_rise_obs_K=float(d["T_C"].max() - d["T_C"][0]), t_end_obs_h=float(d["t_s"][-1] / 3600),
                   fQ_obs_over_C20=fQo / c20, fQ_target_spec=HG2_OBS_FQ.get(r, {}).get(T))
        for vn in variants:
            m, m25 = res[(vn, T, r)], res[(vn, 25, r)]
            sfx = "" if vn == "L2" else f"_{vn}"
            fq = m["Ah"] / m25["Ah"]
            rec.update({f"fQ_model{sfx}": fq, f"fE_model{sfx}": m["Wh"] / m25["Wh"], f"rmse_V{sfx}": m["rmse_V"],
                        f"Q_model_Ah{sfx}": m["Ah"], f"V_mean_model{sfx}": m["V_mean"],
                        f"tte_model_h{sfx}": m["tte_s"] / 3600, f"event{sfx}": m["event"],
                        f"T_rise_model_K{sfx}": m["T_end_C"] - T, f"dfQ{sfx}": fq - fQo})
        # pre-registered pass rules (0.5C and 2C only): T >= 0: |dfQ| <= 0.05 and RMSE_V <= 50 mV; -20 C: 0.10
        if r in ("0.5C", "2C"):
            for sfx in ("", "_X9"):
                dfq, rv = abs(rec[f"dfQ{sfx}"]), rec[f"rmse_V{sfx}"]
                if T >= 0:
                    ok_f, ok_v = (dfq <= 0.05) if T != 25 else None, rv <= 0.050
                elif T == -20:
                    ok_f, ok_v = dfq <= 0.10, None
                else:
                    ok_f, ok_v = None, None  # -10 C: no pre-registered criterion
                rec[f"pass_fQ{sfx}"], rec[f"pass_rmse{sfx}"] = ok_f, ok_v
        rec["in_sample"] = "OCV [CAL] from this curve" if (T == 25 and r == "C/20") else (
            "thermal [CAL] from this curve" if (T == 25 and r == "2C") else None)
        out[f"{T}|{r}"] = rec
    # apparent resistance: data vs model (same definition: (V̄_C/20 − V̄_0.5C)/ΔI)
    s = V.load_hg2_summary()["arrhenius_apparent_resistance"]
    app = {}
    for T in HG2_T:
        rC, r5 = out[f"{T}|C/20"], out[f"{T}|0.5C"]
        dI = r5["I_A"] - rC["I_A"]
        app[str(T)] = dict(R_app_obs_ohm=float(s["R_apparent_ohm"][str(T)]),
                           R_app_model_L2_ohm=(rC["V_mean_model"] - r5["V_mean_model"]) / dI,
                           R_app_model_X9_ohm=(rC["V_mean_model_X9"] - r5["V_mean_model_X9"]) / dI,
                           R_sigma_model_ohm=float(cell.R0(0.5, T + K0, p) + cell.R1(T + K0, p)))

    def ea_fit(key):
        Ts = [T for T in HG2_T if T <= 25]
        x = np.array([1 / (T + K0) - 1 / (25 + K0) for T in Ts])
        y = np.log([app[str(T)][key] / app["25"][key] for T in Ts])
        return float(np.sum(x * y) / np.sum(x * x) * 8.314 / 1e3)

    out["apparent_resistance"] = app
    out["Ea_apparent_kJmol"] = dict(obs=ea_fit("R_app_obs_ohm"), model_L2=ea_fit("R_app_model_L2_ohm"),
                                    model_X9=ea_fit("R_app_model_X9_ohm"), R_sigma_model=ea_fit("R_sigma_model_ohm"),
                                    obs_summary=float(s["Ea_kJ_per_mol"]), spec_expected_R_sigma=20.5)
    out["R_sigma_m20C_ohm"] = dict(model=app["-20"]["R_sigma_model_ohm"], spec_expected=0.165,
                                   obs_apparent=app["-20"]["R_app_obs_ohm"])
    # X9 parameter sensitivity at -20 C (kappa_K 33-47 h^-1, E_a,K 25-55 kJ/mol; §10.4 ranges)
    sens = {}
    for r in ("0.5C", "2C"):
        d, d25 = data[(-20, r)], data[(25, r)]
        for lab, x9 in (("kappa33", dict(kappa_h=33.0)), ("kappa47", dict(kappa_h=47.0)), ("Ea25", dict(Ea=25e3)),
                        ("Ea55", dict(Ea=55e3))):
            a = _hg2_run(p, "X9", d, -20, x9=x9)["Ah"]
            b = _hg2_run(p, "X9", d25, 25, x9=x9)["Ah"]
            sens[f"-20|{r}|{lab}"] = a / b
    out["X9_sensitivity_fQ_m20C"] = sens
    out["calibration"] = dict(Q_Ah=cal["Q_Ah"], C_th_JK=cal["C_th"], R_th_KW=cal["R_th"], rmse_T_K=cal["rmse_T_K"],
                              R_sigma_ohm=0.038, split="6:1",
                              note="25 °C only: OCV = V + I·R_Σ on the C/20 curve; C_th, R_th fitted to the 2C "
                                   "temperature rise (validation.hg2_calibration).  All activation energies RP-63.")
    out["dynamic_UDDS"] = dict(status="not_run", reason="no UDDS (or other drive-cycle) file in "
                               "docs/research/data/validation/lg_hg2_raw: only the C/20, 0.5C and 2C discharges were "
                               "extracted; downloading was out of scope")
    out["solver"] = "analysis/validation_cells.py (core cell functions, constant current, LSODA rtol 1e-7)"
    out["_runtime_s"] = _t() - t0
    return out, dict(data=data, res=res)


# =========================================================================================== V3 NASA
def v3_nasa() -> dict:
    t0 = _t()
    idx = {c["file"]: c for c in V.nasa_index()["curves"]}

    def obs(fn, vcut):
        d = V.load_nasa(fn)
        on = d["I_A"] > 0.5 * d["I_A"].max()
        I = float(np.median(d["I_A"][on]))
        t_on = d["t_s"][on]
        below = np.nonzero(on & (d["V"] <= vcut))[0]
        t_cut = float(d["t_s"][below[0]] - t_on[0]) if below.size else math.nan
        return d, I, t_cut, float(t_on[0])

    out = {}
    cases = dict(fresh=dict(), aged=dict(aged_cell=True), Q_only=dict(aged_cell=False),
                 R_only=dict(aged_cell=True, q_loss=0.0))
    q_aged = 1.0 - phonebatt.params.NASA_B0005["Q_aged"] / phonebatt.params.NASA_B0005["Q_fresh"]
    cases["Q_only"]["q_loss"] = q_aged
    cases["Q_only"]["r_growth"] = 0.0
    df, If, tcf, tonf = obs("nasa_B0005_24C_2A_CC_fresh.csv", 2.7)
    da, Ia, tca, tona = obs("nasa_B0005_24C_2A_CC_aged.csv", 2.7)
    res = {}
    for k, kw in cases.items():
        ov = {kk: vv for kk, vv in kw.items() if kk != "aged_cell"}
        p = nasa_b0005(kw.get("aged_cell", False), **ov)
        I = Ia if k in ("aged", "R_only") else If
        m = run_to_event(CellModel(p, "L2"), "I", I, 24.0)
        res[k] = m
        out[f"model_{k}"] = dict(tte_min=m.tte_s / 60, event=m.event, Ah=m.Ah, T_rise_K=float(m.T_C.max() - 24.0),
                                 R0_ohm=float(p.R0_ref), R1_ohm=float(p.R1_ref), Q_max_Ah=float(p.Q_max))
    obs_f = idx["nasa_B0005_24C_2A_CC_fresh.csv"]["t_load_end_s"] / 60
    obs_a = idx["nasa_B0005_24C_2A_CC_aged.csv"]["t_load_end_s"] / 60
    out["obs"] = dict(tte_fresh_min=obs_f, tte_aged_min=obs_a, ratio_pct=100 * (obs_a / obs_f - 1),
                      t_to_2p7V_fresh_min=tcf / 60, t_to_2p7V_aged_min=tca / 60,
                      Q_fresh_Ah=idx["nasa_B0005_24C_2A_CC_fresh.csv"]["Q_Ah"],
                      Q_aged_Ah=idx["nasa_B0005_24C_2A_CC_aged.csv"]["Q_Ah"],
                      T_rise_fresh_K=idx["nasa_B0005_24C_2A_CC_fresh.csv"]["Tmax"] - idx["nasa_B0005_24C_2A_CC_fresh.csv"]["T0"],
                      I_fresh_A=If, I_aged_A=Ia,
                      note="observed TTE = load duration t_load_end_s (FINAL_SPEC §8 V3: 55.8 -> 39.7 min)")
    tf = res["fresh"].tte_s / 60
    out["model_ratio_pct"] = 100 * (res["aged"].tte_s / res["fresh"].tte_s - 1)
    out["model_ratio_Q_only_pct"] = 100 * (res["Q_only"].tte_s / res["fresh"].tte_s - 1)
    out["model_ratio_R_only_pct"] = 100 * (res["R_only"].tte_s / res["fresh"].tte_s - 1)
    out["err_ratio_pp"] = out["model_ratio_pct"] - out["obs"]["ratio_pct"]
    out["err_fresh_pct"] = 100 * (tf / obs_f - 1)
    out["err_aged_pct"] = 100 * (res["aged"].tte_s / 60 / obs_a - 1)
    out["spec_capacity_only_pct"] = -28.5
    # V(t) RMSE on the load phase
    for k, d, ton in (("fresh", df, tonf), ("aged", da, tona)):
        on = d["I_A"] > 0.5 * d["I_A"].max()
        tt = d["t_s"][on] - ton
        out[f"rmse_V_{k}"] = rmse(res[k].V_at(tt), d["V"][on])
    # temperature sanity checks with the fresh B0005 parameters (other cells of the same make)
    san = {}
    for fn, Ta, vcut in (("nasa_B0055_4C_2A_CC_fresh.csv", 4.0, 2.5), ("nasa_B0029_43C_4A_CC_fresh.csv", 43.0, 2.0),
                         ("nasa_B0047_4C_1A_CC_fresh.csv", 4.0, 2.5)):
        d, I, tc, _ = obs(fn, vcut)
        p = nasa_b0005(False, V_cut=vcut)
        m = run_to_event(CellModel(p, "L2"), "I", I, Ta)
        Qo = idx[fn]["Q_Ah"]
        san[fn.split("_")[1]] = dict(T_amb_C=Ta, I_A=I, V_cut=vcut, Q_obs_Ah=Qo, Q_model_Ah=m.Ah,
                                     ratio_obs_to_B0005_24C=Qo / idx["nasa_B0005_24C_2A_CC_fresh.csv"]["Q_Ah"],
                                     ratio_model_to_B0005_24C=m.Ah / res["fresh"].Ah, event=m.event,
                                     T_rise_model_K=float(m.T_C.max() - Ta),
                                     T_rise_obs_K=idx[fn]["Tmax"] - idx[fn]["T0"])
    san["spec_B0055_ratio"] = 0.71
    out["temperature_sanity"] = san
    out["note"] = ("B0005: R0 <- EIS R_e, R1 <- R_ct, Q_max <- integrated 2 A capacity (params.nasa_b0005); Ai2020 "
                   "4.2 V OCV (U(0) = 3.018 V), V_cut 2.7 V, RP-63 activation energies, cell-only heat balance")
    out["_runtime_s"] = _t() - t0
    return out


# =========================================================================================== V4 Oxford
def v4_oxford() -> tuple[dict, dict]:
    t0 = _t()
    d = V.load_oxford()
    t, Vd, I, Tc = d["t_s"], d["V"], d["I_A"], d["T_C"]
    Q = 0.74
    dt = np.diff(np.concatenate([[0.0], t]))
    Ah = np.cumsum(I * dt) / 3600.0
    z = 1.0 - (Ah - I * dt / 3600.0) / Q  # charge removed before each sample
    low = np.abs(I) < 0.02
    coef = np.polyfit(z[low], Vd[low], 3)  # [CAL] OCV from the low-current samples
    zg = np.linspace(0, 1, 201)
    zr = (z[low].min(), z[low].max())
    Ug = np.polyval(coef, np.clip(zg, *zr))
    cal_rmse = float(np.sqrt(np.mean((np.polyval(coef, z[low]) - Vd[low]) ** 2)))
    # cell ECM scaled to 0.74 Ah by R·Q = const (A12), cell-level values (no PCM): R0·Q = 19 mΩ × 4.70 Ah
    # ([RUN] cell, §4.1), R1·Q = 10 mΩ × 4.70 Ah (RP-63 R1_ref), τ1 = 30 s, RP-63 activation energies.
    base = dict(Q_n=Q, ocv_table=(zg, Ug), k_S=0.0, V_cut=2.5, d_sys=0.0, dV_hv=0.0, P_susp=0.0, P_awake=0.0,
                E_wake=0.0, P_d0=0.0, p_c=0.0, k_RR=0.0, P_SoC_max=0.0, m_b=0.0, burst=False, throttling=False,
                load_heat=False)
    pars = {
        "L2_cellRQ": (replace(RP63, **base, R0_ref=0.019 * 4.70 / Q, R1_ref=0.010 * 4.70 / Q), "L2"),
        "L1_cellRQ": (replace(RP63, **base, R0_ref=0.019 * 4.70 / Q, R1_ref=0.010 * 4.70 / Q), "L1"),
        "L2_HG2RQ": (replace(RP63, **base, R0_ref=0.038 * 6 / 7 * 2.781 / Q, R1_ref=0.038 / 7 * 2.781 / Q, A_R=0.0),
                     "L2"),
        "OCV_only": (replace(RP63, **base, R0_ref=1e-9, R1_ref=1e-9), "L1"),
    }
    out, traces = {}, {}
    for k, (p, var) in pars.items():
        r = CellModel(p, var).simulate_profile(t, I, 40.0, z0=1.0, T_data_C=Tc, h=0.25)
        e = r["V"] - Vd
        hi = np.abs(I) > 2.0
        out[k] = dict(rmse_V=float(np.sqrt(np.mean(e ** 2))), bias_V=float(np.mean(e)),
                      rmse_V_high_current=float(np.sqrt(np.mean(e[hi] ** 2))),
                      R0_40C_ohm=float(cell.R0(0.6, 40 + K0, p)), R1_40C_ohm=float(cell.R1(40 + K0, p)),
                      tau1_40C_s=float(cell.tau1(40 + K0, p)))
        traces[k] = r["V"]
    dI, dV = np.diff(I), np.diff(Vd)
    big = np.abs(dI) > 1.0
    Rstep = -dV[big] / dI[big]
    res = dict(rmse_V=out["L2_cellRQ"]["rmse_V"], variants=out,
               obs=dict(step_resistance_1s_median_ohm=float(np.median(Rstep)),
                        step_resistance_1s_p10_p90_ohm=[float(x) for x in np.percentile(Rstep, [10, 90])],
                        n_steps=int(big.sum()), I_peak_A=float(I.max()), Ah_net=float(Ah[-1]),
                        V_range=[float(Vd.min()), float(Vd.max())], T_range_C=[float(Tc.min()), float(Tc.max())]),
               calibration=dict(kind="OCV = cubic in coulomb-counted z through the |I| < 0.02 A samples",
                                n_low=int(low.sum()), z_range=list(map(float, zr)), fit_rmse_V=cal_rmse,
                                coef=coef.tolist()),
               settings="cell temperature prescribed from the log; RK4 h = 0.25 s, zero-order hold of the 1 Hz "
                        "current; primary = L2_cellRQ (R0·Q = 89.3 mΩ·Ah, R1·Q = 47 mΩ·Ah, τ1 = 30 s)",
               _runtime_s=_t() - t0)
    return res, dict(t=t, V=Vd, I=I, traces=traces)


# =========================================================================================== V5 PyBaMM
AI_RUN = dict(  # [RUN] research fits to the Ai2020 DFN (docs/research/data/cell_ecm_pybamm/pybamm_rc_fit2.json)
    r1=dict(T25=(41.53e-3, 14.36e-3, 34.77), T0=(49.02e-3, 39.97e-3, 98.99)),
    r2=dict(T25=(40.74e-3, 2.51e-3, 6.74, 14.08e-3, 47.43), T0=(47.7e-3, 3.81e-3, 11.95, 39.93e-3, 124.05)))
P1C = 2.28 * 3.820  # W: 1C-equivalent power (nominal 2.28 Ah × Ai2020 mean OCV 3.820 V)


def _ea(v25, v0):
    return 8.314 * math.log(v0 / v25) / (1 / 273.15 - 1 / 298.15)


def v5_pybamm() -> tuple[dict, dict]:
    t0 = _t()
    warnings.filterwarnings("ignore")
    import pybamm

    pv0 = pybamm.ParameterValues("Ai2020")
    x0, x100, _, _ = pybamm.lithium_ion.get_min_max_stoichiometries(pv0)
    Qe = float(pv0.evaluate(pybamm.LithiumIonParameters().n.Q_init))
    Qw = (x100 - x0) * Qe  # capacity of the z window (x0 -> x100), = C/50 capacity
    Vmin = float(pv0["Lower voltage cut-off [V]"])
    a, b = AI_RUN["r1"]["T25"], AI_RUN["r1"]["T0"]
    Ea0, Ea1, Eat = _ea(a[0], b[0]), _ea(a[1], b[1]), _ea(a[2], b[2])
    a2, b2 = AI_RUN["r2"]["T25"], AI_RUN["r2"]["T0"]
    x5 = dict(R0=a2[0], R1=a2[1], tau1=a2[2], R2=a2[3], tau2=a2[4], Ea_R1=_ea(a2[1], b2[1]), Ea_tau1=_ea(a2[2], b2[2]),
              Ea_R2=_ea(a2[3], b2[3]), Ea_tau2=_ea(a2[4], b2[4]))
    base = dict(Q_n=Qw, dV_hv=0.0, dU=0.0, V_cut=Vmin, d_sys=0.0, R0_ref=a[0], R1_ref=a[1], tau1_ref=a[2],
                Ea_R0=Ea0, Ea_R1=Ea1, Ea_tau=Eat, A_R=0.80, z0=0.139, P_susp=0.0, P_awake=0.0, E_wake=0.0, P_d0=0.0,
                p_c=0.0, k_RR=0.0, P_SoC_max=0.0, m_b=0.0, burst=False, throttling=False, load_heat=False,
                C_th=1e12)
    pa = replace(RP63, **base)
    pa_rp = replace(pa, Ea_R0=RP63.Ea_R0, Ea_R1=RP63.Ea_R1, Ea_tau=RP63.Ea_tau)
    # X5 Ea for R0: 2RC fit R0 25/0 C
    x5["Ea_R0"] = _ea(a2[0], b2[0])
    pa5 = replace(pa, Ea_R0=x5["Ea_R0"])
    models = {"L1": CellModel(pa, "L1", isothermal=True), "L2": CellModel(pa, "L2", isothermal=True),
              "X5": CellModel(pa5, "X5", isothermal=True, x5=x5), "X9": CellModel(pa, "X9", isothermal=True),
              "L2_rp63Ea": CellModel(pa_rp, "L2", isothermal=True)}

    def dfn(steps, T, soc0):
        pv = pybamm.ParameterValues("Ai2020")
        pv["Ambient temperature [K]"] = 273.15 + T
        pv["Initial temperature [K]"] = 273.15 + T
        exp = pybamm.Experiment(steps, period="2 seconds")
        try:
            sol = pybamm.Simulation(pybamm.lithium_ion.DFN(), parameter_values=pv, experiment=exp).solve(
                initial_soc=soc0)
        except Exception as ex:  # solver failure is reported, not hidden
            return None, str(ex)[:120]
        return dict(t=sol["Time [s]"].entries, V=sol["Voltage [V]"].entries,
                    Q=sol["Discharge capacity [A.h]"].entries), None

    out = {}
    curves = {}
    # ---- constant-power discharges
    for T in (25, 0, -10):
        for c in (0.1, 0.5, 1.0):
            P = c * P1C
            key = f"cp_{c:g}C_{T}C"
            tt = _t()
            ref, err = dfn([f"Discharge at {P:.4f} W until {Vmin} V"], T, 1.0)
            rec = dict(P_W=P, T_C=T, c_rate_eq=c)
            if ref is None:
                rec.update(status="dfn_failed", error=err)
                out[key] = rec
                continue
            tte_ref = float(ref["t"][-1])
            rec.update(tte_dfn_h=tte_ref / 3600, Q_dfn_Ah=float(ref["Q"][-1]), dfn_s=_t() - tt)
            dtte, rm, tte_m = {}, {}, {}
            sims = {}
            for vn, m in models.items():
                s = run_to_event(m, "P", P, float(T))
                sims[vn] = s
                tte_m[vn] = s.tte_s / 3600
                dtte[vn] = 100 * (s.tte_s / tte_ref - 1)
                Vm = s.V_at(ref["t"])
                rm[vn] = rmse(Vm, ref["V"])
            rec.update(dtte_pct=dtte, rmse_V=rm, tte_model_h=tte_m,
                       pass_3pct={k: (abs(v) <= 3.0) if T >= 0 else None for k, v in dtte.items()})
            out[key] = rec
            curves[key] = dict(ref=ref, sims=sims)
            print(f"  V5 {key}: DFN {tte_ref / 3600:.3f} h; " + ", ".join(f"{k} {v:+.2f}%" for k, v in dtte.items()),
                  flush=True)
    # ---- 10 s pulses (1C-equivalent power) from rest at z0 in {0.8, 0.5, 0.2}
    for T in (25, 0, -10):
        for z0 in (0.8, 0.5, 0.2):
            key = f"pulse_10s_1C_z{z0:g}_{T}C"
            ref, err = dfn(["Rest for 60 seconds", f"Discharge at {P1C:.4f} W for 10 seconds or until {Vmin} V",
                            "Rest for 300 seconds"], T, z0)
            rec = dict(P_W=P1C, T_C=T, z0=z0)
            if ref is None:
                rec.update(status="dfn_failed", error=err)
                out[key] = rec
                continue
            segs = [(60.0, "I", 0.0), (10.0, "P", P1C), (300.0, "I", 0.0)]
            i_end = np.searchsorted(ref["t"], 69.99)
            dV_ref = float(ref["V"][0] - ref["V"][i_end])
            rm, dvp, sims = {}, {}, {}
            for vn, m in models.items():
                s = m.simulate(segs, float(T), z0=z0, max_step=1.0)
                sims[vn] = s
                Vm = s.V_at(ref["t"])
                rm[vn] = rmse(Vm, ref["V"])
                dvp[vn] = float(Vm[0] - s.V_at([69.99])[0])
            rec.update(rmse_V=rm, dV_pulse_end_V=dvp, dV_pulse_end_dfn_V=dV_ref)
            out[key] = rec
            curves[key] = dict(ref=ref, sims=sims)
    # ---- check that the local L2 equals the core LSODA (constant power, isothermal via C_th = 1e12)
    s = run_to_event(models["L2"], "P", 0.5 * P1C, 25.0)
    r = run(UsageState("cp", P_extra=0.5 * P1C), 25.0, p=pa)
    core_chk = dict(local_h=s.tte_s / 3600, core_h=r.tte_h, rel_diff=s.tte_s / 3600 / r.tte_h - 1)
    summ = {}
    for vn in models:
        for band, cond in (("ge0C", lambda T: T >= 0), ("m10C", lambda T: T < 0)):
            vals = [abs(v["dtte_pct"][vn]) for v in out.values() if "dtte_pct" in v and cond(v["T_C"])]
            summ[f"{vn}|{band}"] = dict(max_abs_dtte_pct=max(vals) if vals else None, n=len(vals))
    res = dict(cases=out, summary=summ, core_L2_check=core_chk,
               settings=dict(pybamm=pybamm.__version__, parameter_set="Ai2020", model="DFN (isothermal)",
                             Q_window_Ah=Qw, V_cut=Vmin, P_1C_W=P1C,
                             ecm_25C=dict(R0=a[0], R1=a[1], tau1=a[2], A_R=0.8, z0=0.139),
                             ecm_Ea_J_mol=dict(R0=Ea0, R1=Ea1, tau1=Eat), x5=x5,
                             x9=dict(c=0.166, kappa_h=40.0, Ea=35e3),
                             note=("ECM parameters are the research-step [RUN] 1RC/2RC fits to this DFN's 1C pulses at "
                                   "25 and 0 °C (two-point Arrhenius); the 25/0 °C pulse comparisons are therefore "
                                   "near in-sample for R0/R1/tau1, the constant-power TTEs and all -10 °C cases are "
                                   "out-of-sample.  All ECMs isothermal like the DFN; same OCV (CRM06 fit of the "
                                   "Ai2020 OCV) and the same z window.  L2_rp63Ea uses RP-63 activation energies.")),
               _runtime_s=_t() - t0)
    return res, curves


# =========================================================================================== V6 band
def v6_band() -> dict:
    t0 = _t()
    pred = V.v6_predictions()
    out = {}
    for k, v in pred.items():
        rec = dict(model_h=v["model_h"], obs_h=v["obs_h"], err_pct=v["err_pct"], pass_=v["pass_"], src=v["src"],
                   P_mean_W=v["P_mean_W"], proto_model_h=PROTO_V6[k],
                   dev_vs_proto_pct=100 * (v["model_h"] / PROTO_V6[k] - 1))
        rec["pass"] = rec.pop("pass_")
        out[k] = rec
    # X14 (audio offloaded to a DSP, a = 0) on the marginal audio case
    au = replace(PRIMITIVES["AU"], awake=0)
    r = run(au, 25.0, p=iphone14pro())
    out["14P_audio_BT"]["X14_model_h"] = r.tte_h
    out["14P_audio_BT"]["X14_err_pct"] = 100 * (r.tte_h / 75.0 - 1)
    out["14P_audio_BT"]["X14_P_W"] = float(demand_power(au, iphone14pro()))
    # power bands
    gm = float(steady_power("GM", RP63, 25.0))
    out["max_load_gaming"] = dict(model_W=gm, obs_band_W=[2.4, 5.8], proto_W=3.12, pass_=None, src="NBC22-24")
    out["max_load_gaming"]["pass"] = bool(2.4 <= gm <= 5.8)
    out["max_load_gaming"].pop("pass_")
    s1 = run("SB-W", 25.0)
    rate = 100.0 / s1.tte_h
    out["standby_drain"] = dict(model_pct_per_h=rate, obs_band_pct_per_h=[0.15, 0.5], proto_pct_per_h=0.23,
                                pass_=None, src="device_power BG9", S1_tte_h=s1.tte_h)
    out["standby_drain"]["pass"] = bool(0.15 <= rate <= 0.5)
    out["standby_drain"].pop("pass_")
    n_pass = sum(bool(v["pass"]) for v in out.values() if isinstance(v, dict) and "pass" in v)
    n_tot = sum(1 for v in out.values() if isinstance(v, dict) and "pass" in v)
    out["n_pass"], out["n_cases"] = n_pass, n_tot
    out["proto_n_pass"] = 11
    out["_runtime_s"] = _t() - t0
    return out


# =========================================================================================== V7 cold phenomenology
def v7_cold() -> dict:
    t0 = _t()
    T20 = 20.0 + K0
    f0 = cell.a_R(T20, RP63.Ea_R0)
    f1 = cell.a_R(T20, RP63.Ea_R1)
    R1 = 0.25 / (6 * f0 + f1)  # R0 + R1 = 0.25 Ω at 20 °C, z = 0.5, RP-63 6:1 split
    p = device(2.7, 5.2, R0_ref=6 * float(R1), R1_ref=float(R1))
    Eu = float(cell.usable_energy(p))
    P_bar = Eu / (159.0 / 60.0)  # workload intensity from LEE20's 159 min at 25 °C (energy bucket, an input)
    ui0 = UsageState("UI-2015", s_scr=1, awake=1, b=0.5, APL=0.5, f_r=60, x_c=0.0, w=1, c="lte", rho="wifi",
                     lam=2.0, d_dl=1.0, dPb=3.0)
    ui = replace(ui0, P_extra=float(P_bar - demand_power(ui0, p)))
    out = dict(config=dict(Q_Ah=2.7, diag_in=5.2, R0_ref_ohm=float(p.R0_ref), R1_ref_ohm=float(p.R1_ref),
                           R_sigma_20C_ohm=float(cell.R0(0.5, T20, p) + cell.R1(T20, p)), E_u_Wh=Eu,
                           P_mean_W=P_bar, dPb_W=3.0,
                           I_peak_A_at_3p8V=(P_bar + 3.0) / 3.8,
                           note="workload P̄ = E_u/159 min (LEE20 25 °C operating time, energy bucket): an input "
                                "intensity, not a fitted model parameter"))
    runs = {}
    for T in (25, 10, 0, -5, -10, -15):
        r = run(ui, float(T), p=p)
        runs[str(T)] = dict(tte_min=60 * r.tte_h, event=r.event, s_off_pct=100 * r.s_off, T_end_C=r.T_end_C)
    out["runs"] = runs
    loss0 = runs["0"]["tte_min"] / runs["25"]["tte_min"] - 1
    out["tte_loss_0C_pct"] = 100 * loss0
    out["tte_loss_0C_obs_pct"] = -15.1
    out["pass_tte_loss_same_order"] = bool(-30.2 <= 100 * loss0 <= -7.55)  # within a factor 2 of -15.1 %
    s_cold = {T: runs[T]["s_off_pct"] for T in ("0", "-5", "-10", "-15")}
    out["pass_brownout_soc_15_35_at_m5C"] = bool(runs["-5"]["event"] == "brownout" and 15 <= s_cold["-5"] <= 35)
    out["brownout_soc_at_le0C_pct"] = s_cold
    # stranded charge comes back after warming (C12): shut off at -5 °C, warm to 25 °C (z kept, v1 relaxed)
    r_c = run(ui, -5.0, p=p)
    r_w = run(ui, 25.0, s0=r_c.s_off, T0_C=25.0, p=p)
    Qrec = (r_c.z_end - r_w.z_end) * float(p.Q_max) * 1000.0
    out["recovery_after_warming"] = dict(s_off_cold_pct=100 * r_c.s_off, extra_min=60 * r_w.tte_h, extra_mAh=Qrec,
                                         obs_mAh=330.0, obs_min=18.0, event_warm=r_w.event,
                                         pass_same_order=bool(165.0 <= Qrec <= 660.0))
    # sensitivities of the (read-off / [unverified]) configuration -- not fitted, reported for interpretation
    def variant(pp, uu):
        rr = {str(T): run(uu, float(T), p=pp) for T in (25, 0, -5, -10)}
        rc = rr["-5"]
        rw = run(uu, 25.0, s0=rc.s_off, T0_C=25.0, p=pp) if rc.event != "none" else None
        return dict(tte_loss_0C_pct=100 * (rr["0"].tte_h / rr["25"].tte_h - 1),
                    s_off_pct={T: 100 * r.s_off for T, r in rr.items()}, events={T: r.event for T, r in rr.items()},
                    recovered_mAh=(rc.z_end - rw.z_end) * float(pp.Q_max) * 1000.0 if rw else None)

    ui_h = replace(ui0, P_extra=float(0.5 * P_bar - demand_power(ui0, p)))
    R1b = 0.13 / (6 * f0 + f1)
    p_lo = device(2.7, 5.2, R0_ref=6 * float(R1b), R1_ref=float(R1b))
    out["sensitivity"] = {
        "half_mean_load": variant(p, ui_h),
        "R_sigma_0.13ohm_LEE20_C5": variant(p_lo, replace(ui0, P_extra=float(P_bar - demand_power(ui0, p_lo)))),
        "X7_no_burst_envelope": variant(replace(p, burst=False), ui),
    }
    out["_runtime_s"] = _t() - t0
    return out


# =========================================================================================== figures
def fig_bonn(V1, F):
    import matplotlib.pyplot as plt
    from style import apply, save, BLUE, ORANGE, AQUA, YELLOW, INK, INK2, MUTED, GRID, W1, panel_label
    apply()
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(W1, 2.7), gridspec_kw=dict(width_ratios=[1.25, 1.0], wspace=0.32))
    col = {"14 Pro": BLUE, "13 Pro": ORANGE}
    mk = {"idle": "o", "cam": "s", "gps": "^", "cxl": "D"}
    lo, hi = 5.5, 80
    xx = np.array([lo, hi])
    ax.fill_between(xx, xx * 0.85, xx * 1.15, color=GRID, alpha=0.7, lw=0, label="±15 %")
    ax.plot(xx, xx, color=MUTED, lw=0.8)
    for dev in DEVS:
        for case in CASES:
            r = V1[f"{dev}|{case}"]
            u = r["uq_cal_fixed"]
            ax.errorbar(r["obs_h"], r["model_h"], yerr=[[r["model_h"] - u["q05_h"]], [u["q95_h"] - r["model_h"]]],
                        xerr=[[r["obs_h"] - r["obs_min_h"]], [r["obs_max_h"] - r["obs_h"]]], fmt=mk[case],
                        color=col[dev], mfc="white" if (dev, case) == ("14 Pro", "idle") else col[dev],
                        mec=col[dev], ms=4.5,
                        elinewidth=0.9, capsize=0, zorder=3)
    ax.set_xscale("log"), ax.set_yscale("log")
    ax.set_xlim(lo, hi), ax.set_ylim(lo, hi)
    ticks = [6, 10, 20, 30, 50]
    ax.set_xticks(ticks), ax.set_yticks(ticks)
    ax.set_xticklabels([str(t) for t in ticks]), ax.set_yticklabels([str(t) for t in ticks])
    ax.minorticks_off()
    ax.set_xlabel("observed TTE (h), median and device range")
    ax.set_ylabel("model TTE (h), 90 % parameter UQ")
    # direct labels per case cluster
    lab = {"cam": (8.6, 13.5, "camera"), "cxl": (12.5, 6.4, "camera+GPS"), "idle": (38, 66, "idle (14 Pro = CAL)"),
           "gps": (31, 16.0, "GPS")}
    for case, (x, y, s) in lab.items():
        ax.text(x, y, s, fontsize=7, color=INK2, ha="center")
    ax.annotate("14 Pro GPS: −48 %\n(pre-registered fail)", xy=(48.5, 25.0), xytext=(52, 11),
                fontsize=6.8, color=INK2, ha="center", arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6))
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker="o", ls="", color=col[d], label=f"iPhone {d}") for d in DEVS]
    h.append(Line2D([], [], color=GRID, lw=5, label="±15 % band"))
    ax.legend(handles=h, loc="upper left", fontsize=7)
    panel_label(ax, "(a)")
    # (b) slope test
    tr = F["traces"]
    xs = np.arange(3)
    obs_r = np.array([v["ratio_obs"] for v in tr.values()])
    for i, rr in enumerate(obs_r):
        bx.plot(xs + (i % 7 - 3) * 0.025, rr, ls="", marker="o", ms=2.2, color=MUTED, alpha=0.7, zorder=2)
    bx.plot(xs, np.median(obs_r, axis=0), color=INK, lw=1.4, marker="o", ms=3.5, label="observed (28 traces, median)",
            zorder=4)
    d1 = np.median([v["D1"]["ratio_model"] for v in tr.values()], axis=0)
    d2 = np.median([v["D2"]["ratio_model"] for v in tr.values()], axis=0)
    e1, e2 = V1["slope_test"]["D1"]["rms_err_pct_mean"], V1["slope_test"]["D2"]["rms_err_pct_mean"]
    bx.plot(xs, d1, color=AQUA, lw=1.8, marker="s", ms=3.5, label=f"model, coulombic D1 (err {e1:.1f} %)", zorder=3)
    bx.plot(xs, d2, color=YELLOW, lw=1.8, marker="^", ms=3.5, label=f"model, energy D2 (err {e2:.1f} %)", zorder=3)
    bx.set_xticks(xs)
    bx.set_xticklabels(["99→50 %", "50→20 %", "20→1 %"])
    bx.set_xlabel("displayed-SOC window")
    bx.set_ylabel("drain time per % / run mean")
    bx.set_ylim(0.83, 1.13)
    bx.legend(loc="lower left", fontsize=6.8)
    panel_label(bx, "(b)")
    save(fig, "fig_validation_bonn")


def fig_hg2(V2, F):
    import matplotlib.pyplot as plt
    from style import apply, save, BLUE, ORANGE, INK, INK2, MUTED, GRID, W1, panel_label, SEQ_BLUE
    apply()
    fig, axs = plt.subplots(1, 3, figsize=(W1, 2.45), gridspec_kw=dict(wspace=0.34))
    Ts = sorted(HG2_T)
    for ax, r, lab in ((axs[0], "0.5C", "(a)"), (axs[1], "2C", "(b)")):
        fo = [V2[f"{T}|{r}"]["fQ_obs"] for T in Ts]
        tol = [0.10 if T == -20 else (0.05 if T >= 0 else np.nan) for T in Ts]
        for T, f, tl in zip(Ts, fo, tol):
            if np.isfinite(tl) and T != 25:
                ax.plot([T, T], [f - tl, f + tl], color=GRID, lw=5, solid_capstyle="butt", zorder=1)
        ax.plot(Ts, fo, ls="", marker="o", color=INK, ms=4, label="data (HG2)", zorder=4)
        ax.plot(Ts, [V2[f"{T}|{r}"]["fQ_model"] for T in Ts], color=BLUE, marker="s", ms=3, label="L2 (1RC)", zorder=3)
        ax.plot(Ts, [V2[f"{T}|{r}"]["fQ_model_X9"] for T in Ts], color=ORANGE, marker="^", ms=3, label="X9 (KiBaM)",
                zorder=3)
        ax.set_xlabel("ambient temperature (°C)")
        ax.set_ylim(0.4, 1.05)
        ax.set_xticks([-20, -10, 0, 10, 25, 40])
        ax.text(0.5, 0.04 if r == "0.5C" else 0.11, f"{r} constant current", transform=ax.transAxes, fontsize=7,
                color=INK2, ha="center")
        panel_label(ax, lab)
    axs[0].set_ylabel("Q(T) / Q(25 °C)")
    axs[0].legend(loc="center right", fontsize=6.8)
    axs[1].text(0.5, 0.04, "gray bars: pass tolerance", transform=axs[1].transAxes, fontsize=6.5, color=INK2,
                ha="center")
    # (c) V vs delivered charge at 2C
    ax = axs[2]
    data, res = F["data"], F["res"]
    from matplotlib.lines import Line2D
    hh = []
    for T, mkr in ((25, "o"), (0, "s"), (-20, "^")):
        d = data[(T, "2C")]
        ax.plot(d["Ah"], d["V"], ls="", marker=mkr, ms=2.4, color=INK, zorder=4)
        hh.append(Line2D([], [], ls="", marker=mkr, ms=3, color=INK, label=f"{T} °C"))
        for vn, c in (("L2", BLUE), ("X9", ORANGE)):
            s = res[(vn, T, "2C")]["sim"]
            q = np.concatenate([[0.0], np.cumsum(0.5 * (s.I[1:] + s.I[:-1]) * np.diff(s.t)) / 3600])
            ax.plot(q, s.V, color=c, lw=1.3, zorder=3)
    ax.set_xlabel("delivered charge (A h)")
    ax.set_ylabel("terminal voltage (V)")
    ax.set_xlim(0, 3.75)
    ax.set_xticks([0, 1, 2, 3])
    ax.set_ylim(2.75, 4.2)
    hh += [Line2D([], [], color=BLUE, label="L2"), Line2D([], [], color=ORANGE, label="X9")]
    ax.legend(handles=hh, loc="upper right", fontsize=6.2, title="2C data", title_fontsize=6.5, handletextpad=0.3,
              borderaxespad=0.1)
    panel_label(ax, "(c)")
    save(fig, "fig_validation_hg2")


def fig_band(V6):
    import matplotlib.pyplot as plt
    from style import apply, save, BLUE, ORANGE, AQUA, INK, INK2, MUTED, GRID, W1, panel_label
    apply()
    order = ["14P_local_video", "14P_stream_video", "14P_audio_BT", "14P_cell_web_150", "14PM_wifi_web_150",
             "14PM_local_video_150", "14PM_idle_minbright", "Pixel8_wifi_web", "Pixel8_local_video",
             "Pixel8_idle_minbright", "Pixel8a_S12_day", "S24_wifi_web", "S24_local_video", "S24_idle_minbright"]
    label = {"14P_local_video": "iPhone 14 Pro · local video (Apple)", "14P_stream_video": "iPhone 14 Pro · streamed video (Apple)",
             "14P_audio_BT": "iPhone 14 Pro · BT audio (Apple)", "14P_cell_web_150": "iPhone 14 Pro · cellular web, 150 nit (LMAG)",
             "14PM_wifi_web_150": "iPhone 14 Pro Max · Wi-Fi web, 150 nit (NBC)",
             "14PM_local_video_150": "iPhone 14 Pro Max · local video (NBC)",
             "14PM_idle_minbright": "iPhone 14 Pro Max · idle, min. brightness (NBC)",
             "Pixel8_wifi_web": "Pixel 8 · Wi-Fi web (NBC)", "Pixel8_local_video": "Pixel 8 · local video (NBC)",
             "Pixel8_idle_minbright": "Pixel 8 · idle, min. brightness (NBC)",
             "Pixel8a_S12_day": "Pixel 8a · typical day vs ≥ 31 h claim (Google)",
             "S24_wifi_web": "Galaxy S24 · Wi-Fi web (NBC)", "S24_local_video": "Galaxy S24 · local video (NBC)",
             "S24_idle_minbright": "Galaxy S24 · idle, min. brightness (NBC)"}
    vendor = lambda k: "Apple" if k.startswith("14P") else ("Google" if k.startswith("Pixel") else "Samsung")
    col = {"Apple": BLUE, "Google": ORANGE, "Samsung": AQUA}
    fig, ax = plt.subplots(figsize=(W1, 3.05))
    y = np.arange(len(order))[::-1]
    ax.axvspan(-25, 25, color=GRID, alpha=0.6, lw=0)
    ax.axvline(0, color=MUTED, lw=0.8)
    for yi, k in zip(y, order):
        r = V6[k]
        e = r["err_pct"]
        c = col[vendor(k)]
        m = ">" if k == "Pixel8a_S12_day" else "o"
        ax.plot(e, yi, marker=m, ls="", ms=5, color=c, mfc=c if r["pass"] else "white", mec=c, mew=1.2, zorder=3)
        if k == "14P_audio_BT":
            ex = r["X14_err_pct"]
            ax.annotate("", xy=(ex, yi), xytext=(e, yi), arrowprops=dict(arrowstyle="->", color=MUTED, lw=0.8))
            ax.text(ex + 3, yi, "X14 DSP audio", va="center", fontsize=6.5, color=INK2)
        ax.text(min(e, 190) + (6 if e < 150 else -6), yi + 0.32, f"{e:+.0f} %", fontsize=6.3, color=INK2,
                ha="left" if e < 150 else "right", va="center")
    ax.set_yticks(y)
    ax.set_yticklabels([label[k] for k in order], fontsize=6.8)
    ax.set_xlabel("model error ε = TTE_model / TTE_published − 1 (%)")
    ax.set_xlim(-45, 200)
    ax.set_ylim(-0.7, len(order) - 0.3)
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker="o", ls="", color=col[v], label=v) for v in ("Apple", "Google", "Samsung")]
    h += [Line2D([], [], marker="o", ls="", color=MUTED, mfc="white", label="outside ±25 %"),
          Line2D([], [], color=GRID, lw=6, label="±25 % band")]
    ax.legend(handles=h, loc="upper right", bbox_to_anchor=(1.0, 0.86), fontsize=6.8)
    g, s = V6["max_load_gaming"], V6["standby_drain"]
    ax.text(198, len(order) - 1.2,
            f"power-band checks: capped gaming {g['model_W']:.2f} W in 2.4–5.8 W ✓;\n"
            f"standby drain {s['model_pct_per_h']:.2f} %/h in 0.15–0.5 %/h ✓   ({V6['n_pass']}/{V6['n_cases']} pass)",
            fontsize=6.6, color=INK2, ha="right", va="top")
    save(fig, "fig_validation_band")


def fig_pybamm(V5, C):
    import matplotlib.pyplot as plt
    from style import apply, save, BLUE, ORANGE, AQUA, YELLOW, INK, INK2, MUTED, GRID, W1, panel_label
    apply()
    vc = {"L2": BLUE, "X9": ORANGE, "L1": AQUA, "X5": YELLOW}
    vl = {"L2": "L2 (1RC)", "X9": "X9 (KiBaM)", "L1": "L1 (quasi-static)", "X5": "X5 (2RC)"}
    fig, axs = plt.subplots(1, 3, figsize=(W1, 2.5), gridspec_kw=dict(wspace=0.36, width_ratios=[1, 1.15, 1]))
    ax = axs[0]
    for T in (25, 0, -10):
        key = f"cp_0.5C_{T}C"
        if key not in C:
            continue
        ref = C[key]["ref"]
        ax.plot(ref["Q"], ref["V"], color=INK, lw=1.6, zorder=4)
        for vn in ("L2", "X9"):
            s = C[key]["sims"][vn]
            q = np.concatenate([[0.0], np.cumsum(0.5 * (s.I[1:] + s.I[:-1]) * np.diff(s.t)) / 3600])
            ax.plot(q, s.V, color=vc[vn], lw=1.2, ls="--" if vn == "X9" else "-", zorder=3)
        qe = float(ref["Q"][-1])
        if T == -10:
            ax.text(qe - 0.08, 3.08, f"{T} °C", fontsize=6.8, color=INK2, ha="right")
        else:
            ax.annotate(f"{T} °C", xy=(qe - 0.02, 3.18 if T == 0 else 3.45), xytext=(1.45, 3.12 if T == 0 else 3.32),
                        fontsize=6.8, color=INK2, arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6))
    ax.set_xlabel("delivered charge (A h), 0.5C-eq. power")
    ax.set_ylabel("terminal voltage (V)")
    ax.set_ylim(2.95, 4.25)
    ax.set_xlim(0, 2.6)
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], color=INK, label="DFN (PyBaMM)"), Line2D([], [], color=BLUE, label="L2"),
                       Line2D([], [], color=ORANGE, ls="--", label="X9")], loc="upper right", fontsize=6.5)
    panel_label(ax, "(a)")
    ax = axs[1]
    keys = [k for k in V5["cases"] if k.startswith("cp_") and "dtte_pct" in V5["cases"][k]]
    xs = np.arange(len(keys))
    ax.axhspan(-3, 3, color=GRID, alpha=0.7, lw=0)
    ax.axhline(0, color=MUTED, lw=0.8)
    for j, vn in enumerate(("L1", "L2", "X5", "X9")):
        yv = [V5["cases"][k]["dtte_pct"][vn] for k in keys]
        ax.plot(xs + (j - 1.5) * 0.14, yv, ls="", marker="o", ms=3.2, color=vc[vn], label=vl[vn], zorder=3)
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{V5['cases'][k]['c_rate_eq']:g}" for k in keys], fontsize=6.5)
    temps = [V5["cases"][k]["T_C"] for k in keys]
    for T in sorted(set(temps), reverse=True):
        ii = [i for i, tt in enumerate(temps) if tt == T]
        ax.text(np.mean(ii), -0.095, f"{T} °C", transform=ax.get_xaxis_transform(), ha="center", va="top",
                fontsize=6.8, color=INK2)
        if ii[0] > 0:
            ax.axvline(ii[0] - 0.5, color=GRID, lw=0.8)
    ax.set_ylabel("TTE difference vs DFN (%)")
    ax.set_yscale("symlog", linthresh=5)
    ax.set_yticks([-50, -10, -3, 0, 3, 10, 50, 300])
    ax.set_yticklabels(["−50", "−10", "−3", "0", "3", "10", "50", "300"])
    ax.set_ylim(-60, 700)
    ax.minorticks_off()
    ax.set_xlabel("C-rate equivalent", labelpad=13)
    ax.legend(loc="upper left", fontsize=6.0, handletextpad=0.2, borderaxespad=0.2)
    panel_label(ax, "(b)")
    ax = axs[2]
    key = "pulse_10s_1C_z0.5_0C"
    if key in C:
        ref = C[key]["ref"]
        ax.plot(ref["t"] - 60, ref["V"], color=INK, lw=1.8, label="DFN", zorder=4)
        for vn in ("L1", "L2", "X5", "X9"):
            s = C[key]["sims"][vn]
            tt = np.linspace(0, ref["t"][-1], 1500)
            ax.plot(tt - 60, s.V_at(tt), color=vc[vn], lw=1.1, label=vl[vn].split(" ")[0], zorder=3)
        ax.set_xlim(-10, 200)
        ax.set_xlabel("time from pulse start (s)")
        ax.set_ylabel("terminal voltage (V)")
        ax.text(0.96, 0.62, "10 s 1C-eq. pulse\n0 °C, z = 0.5", transform=ax.transAxes, fontsize=6.8,
                color=INK2, ha="right")
        ax.legend(loc="lower right", fontsize=6.3, ncol=2)
    panel_label(ax, "(c)")
    save(fig, "fig_validation_pybamm")


# =========================================================================================== table + summary
def _pf(x):
    return "--" if x is None else ("pass" if x else "fail")


def criteria(R) -> dict:
    V1, V2, V3, V4, V5, V6, V7 = (R.get(k, {}) for k in ("V1", "V2", "V3", "V4", "V5", "V6", "V7"))
    C = {}
    for dev in DEVS:
        for case in CASES:
            r = V1[f"{dev}|{case}"]
            if r["holdout"]:
                C[f"V1a|{dev}|{case}"] = dict(metric=f"eps = {r['err_pct']:+.1f} %", rule="|eps| <= 15 %",
                                              pass_=r["pass_a"])
            C[f"V1b|{dev}|{case}"] = dict(metric=f"90% UQ [{r['uq_cal_fixed']['q05_h']:.1f}, "
                                                 f"{r['uq_cal_fixed']['q95_h']:.1f}] h vs obs "
                                                 f"[{r['obs_min_h']:.2f}, {r['obs_max_h']:.2f}] h",
                                          rule="obs range inside 90 % UQ", pass_=r["pass_b"])
    for dev in DEVS:
        a = V1["additivity_residual_pct"][dev]
        C[f"V1c|{dev}"] = dict(metric=f"model {a['model']:+.2f} %, data {a['data_cam60']:+.2f} %",
                               rule="|model residual| <= 1 %", pass_=a["pass_c"])
    st = V1["slope_test"]
    C["V1.4"] = dict(metric=f"RMS shape error D1 {st['D1']['rms_err_pct_mean']:.1f} %, D2 "
                            f"{st['D2']['rms_err_pct_mean']:.1f} %", rule="expected: data favour D2",
                     pass_=st["favours"] == "D2")
    for r in ("0.5C", "2C"):
        for T in (40, 10, 0, 25, -20):
            k = f"{T}|{r}"
            v = V2[k]
            for sfx, nm in (("", "L2"), ("_X9", "X9")):
                if v.get(f"pass_fQ{sfx}") is not None:
                    C[f"V2|{k}|fQ|{nm}"] = dict(metric=f"dfQ = {v[f'dfQ{sfx}']:+.3f}",
                                                rule="|dfQ| <= 0.05 (T >= 0), 0.10 (-20 C)", pass_=v[f"pass_fQ{sfx}"])
                if v.get(f"pass_rmse{sfx}") is not None:
                    C[f"V2|{k}|rmseV|{nm}"] = dict(metric=f"RMSE_V = {1000 * v[f'rmse_V{sfx}']:.0f} mV",
                                                   rule="<= 50 mV (T >= 0)", pass_=v[f"pass_rmse{sfx}"])
    if V3:
        C["V3"] = dict(metric=f"aged/fresh TTE model {V3['model_ratio_pct']:+.1f} % vs obs {V3['obs']['ratio_pct']:+.1f} %",
                       rule="no numeric criterion (structural Q/R test)", pass_=None)
    if V4:
        C["V4"] = dict(metric=f"RMSE_V = {1000 * V4['rmse_V']:.0f} mV", rule="report only", pass_=None)
    if V5:
        for k, v in V5["cases"].items():
            if "dtte_pct" in v and v["T_C"] >= 0:
                for vn in ("L1", "L2", "X5", "X9"):
                    C[f"V5|{k}|{vn}"] = dict(metric=f"dTTE = {v['dtte_pct'][vn]:+.2f} %", rule="|dTTE| <= 3 % (T >= 0)",
                                             pass_=v["pass_3pct"][vn])
    if V6:
        for k, v in V6.items():
            if isinstance(v, dict) and "pass" in v:
                met = (f"eps = {v['err_pct']:+.1f} %" if "err_pct" in v else
                       (f"{v.get('model_W', v.get('model_pct_per_h')):.3f}"))
                C[f"V6|{k}"] = dict(metric=met, rule="|eps| <= 25 % (or inside band / >= bound)", pass_=v["pass"])
    if V7:
        C["V7|tte_loss_0C"] = dict(metric=f"{V7['tte_loss_0C_pct']:+.1f} % vs -15.1 %", rule="same order (factor 2)",
                                   pass_=V7["pass_tte_loss_same_order"])
        C["V7|brownout_soc"] = dict(metric=f"s_off(-5 C) = {V7['brownout_soc_at_le0C_pct']['-5']:.1f} %",
                                    rule="brown-out at 15-35 % displayed", pass_=V7["pass_brownout_soc_15_35_at_m5C"])
        C["V7|recovery"] = dict(metric=f"{V7['recovery_after_warming']['extra_mAh']:.0f} mAh vs 330 mAh",
                                rule="charge comes back (same order)", pass_=V7["recovery_after_warming"]["pass_same_order"])
    for v in C.values():
        v["pass"] = v.pop("pass_")
    return C


def _m(x, fmt="{:+.1f}"):
    """Signed number for LaTeX with a true minus sign."""
    t = fmt.format(x)
    return t.replace("-", "$-$", 1) if t.startswith("-") else t


def _v1b_fails(C) -> str:
    """Compact list of failed V1b cases: a case failing on both phones is listed once."""
    bad = [(k.split("|")[1], k.split("|")[2]) for k in C if k.startswith("V1b") and C[k]["pass"] is False]
    if not bad:
        return "none"
    out = []
    for case in CASES:
        devs = [d for d, c in bad if c == case]
        if devs:
            out.append(f"{CASE_LABEL[case]} ({'both' if len(devs) == 2 else devs[0]})")
    return "fails: " + ", ".join(out)


def write_table(R):
    """paper/tables/tab_validation.tex: one row per pre-registered test (booktabs, no float wrapper)."""
    V1, V2, V3, V4, V5, V6, V7 = (R.get(k) for k in ("V1", "V2", "V3", "V4", "V5", "V6", "V7"))
    C = R["criteria"]

    def cnt(pre, suffix=""):
        ks = [k for k in C if k.startswith(pre) and k.endswith(suffix) and C[k]["pass"] is not None]
        return f"{sum(C[k]['pass'] for k in ks)}/{len(ks)}"

    e = lambda d, c: V1[f"{d}|{c}"]["err_pct"]
    ad = V1["additivity_residual_pct"]
    st = V1["slope_test"]
    rows = [
        ("V1a Bonn", r"$|\varepsilon|\le 15\%$, hold-outs", cnt("V1a"),
         f"GPS fails: {_m(e('14 Pro', 'gps'), '{:+.0f}')}\\% (14 Pro), {_m(e('13 Pro', 'gps'), '{:+.0f}')}\\% (13 Pro)"),
        ("V1b Bonn", r"obs.\ range in 90\% UQ", cnt("V1b"), _v1b_fails(C)),
        ("V1c Bonn", r"additivity $\le 1\%$", cnt("V1c"),
         f"model {_m(ad['14 Pro']['model'], '{:+.2f}')}\\%, data {_m(ad['14 Pro']['data_cam60'], '{:+.2f}')}\\% (14 Pro)"),
        ("V1.4 Bonn", "SOC shape: D1 vs D2", st["favours"],
         f"slope error D1 {st['D1']['rms_err_pct_mean']:.1f}\\%, D2 {st['D2']['rms_err_pct_mean']:.1f}\\%"),
        ("V2 HG2, L2", r"$|\Delta f_Q|$; RMSE$_V\le$ 50 mV", f"{cnt('V2|', '|fQ|L2')}; {cnt('V2|', '|rmseV|L2')}",
         f"$f_Q$($-$20\\,\\textdegree C, 2C) {V2['-20|2C']['fQ_model']:.2f} vs {V2['-20|2C']['fQ_obs']:.2f}"),
        ("V2 HG2, X9", "same", f"{cnt('V2|', '|fQ|X9')}; {cnt('V2|', '|rmseV|X9')}",
         f"$f_Q$($-$20\\,\\textdegree C, 2C) {V2['-20|2C']['fQ_model_X9']:.2f} vs {V2['-20|2C']['fQ_obs']:.2f}"),
        ("V3 NASA", "aged/fresh TTE", "--",
         f"{_m(V3['model_ratio_pct'])}\\% vs {_m(V3['obs']['ratio_pct'])}\\% (charge-limited)"),
        ("V4 Oxford", r"RMSE$_V$ (report)", "--", f"{1000 * V4['rmse_V']:.0f} mV; OCV only "
                                                    f"{1000 * V4['variants']['OCV_only']['rmse_V']:.0f} mV"),
    ]
    if V5:
        mx25 = max(abs(v["dtte_pct"]["L2"]) for v in V5["cases"].values() if "dtte_pct" in v and v["T_C"] == 25)
        rows.append(("V5 DFN, L2; X9", r"$|\Delta\mathrm{TTE}|\le 3\%$, $\ge 0$\,\textdegree C",
                     f"{cnt('V5', '|L2')}; {cnt('V5', '|X9')}",
                     f"L2 $\\le$ {mx25:.1f}\\% at 25\\,\\textdegree C; "
                     f"{_m(V5['cases']['cp_0.5C_0C']['dtte_pct']['L2'])}\\% at 0\\,\\textdegree C, 0.5C"))
    rows += [
        ("V6 band", r"$|\varepsilon|\le 25\%$", f"{V6['n_pass']}/{V6['n_cases']}", "fails: Android floors, BT audio"),
        ("V7 cold", "LEE20 phenomenology", cnt("V7"),
         f"0\\,\\textdegree C loss {_m(V7['tte_loss_0C_pct'], '{:+.0f}')}\\% vs $-$15\\%"),
    ]
    lines = [r"\begin{tabular}{@{}llll@{}}", r"\toprule", r"Test & Criterion & Passed & Key result \\", r"\midrule"]
    lines += [" & ".join(r) + r" \\" for r in rows]
    lines += [r"\bottomrule", r"\end{tabular}", ""]
    write_tex("tab_validation.tex", "\n".join(lines))


def summary_md(R) -> str:
    """results/validation_summary.md: digest for the paper writer (<= 80 lines)."""
    V1, V2, V3, V4, V5, V6, V7 = (R.get(k) for k in ("V1", "V2", "V3", "V4", "V5", "V6", "V7"))
    C = R["criteria"]
    main = [v["pass"] for k, v in C.items() if v["pass"] is not None and not k.endswith(("|X9", "|L1", "|X5"))]
    var = [v["pass"] for k, v in C.items() if v["pass"] is not None and k.endswith(("|X9", "|L1", "|X5"))]
    L = ["# Validation (AN6) summary", "",
         "`analysis/run_validation.py` -> `results/validation.json` (FINAL_SPEC §8, pre-registered). Only the declared "
         "[CAL] steps were fitted (V1.1 ΔP_awake; V2 HG2 OCV/thermal at 25 °C; V4 Oxford OCV). [PROTO] numbers are "
         f"regression targets, not validation. **Main model (L2): {sum(main)} of {len(main)} criteria pass**; "
         f"structural variants (X9, L1, X5): {sum(var)} of {len(var)}.", ""]
    st, ad = V1["slope_test"], V1["additivity_residual_pct"]
    L.append(f"## V1 Bonn iPhones (ΔP_awake 0.090 W; exact CAL {V1['calibration']['P_awake_exact_W']:.4f} W, "
             "[PROTO] 0.0877)")
    L.append("| case | model h | obs h (range) | ε | ε blind 0.10 W | 90 % UQ, CAL fixed | (a) | (b) |")
    L.append("|---|---|---|---|---|---|---|---|")
    for dev in DEVS:
        for case in CASES:
            r = V1[f"{dev}|{case}"]
            u = r["uq_cal_fixed"]
            L.append(f"| {dev} {CASE_LABEL[case]} | {r['model_h']:.2f} | {r['obs_h']:.2f} ({r['obs_min_h']:.1f}–"
                     f"{r['obs_max_h']:.1f}) | {r['err_pct']:+.1f} % | {r['err_blind_pct']:+.1f} % | {u['q05_h']:.1f}–"
                     f"{u['q95_h']:.1f} h | {_pf(r['pass_a'])} | {_pf(r['pass_b'])} |")
    mx = max(abs(V1[f"{d}|{c}"]["dev_vs_proto_pct"]) for d in DEVS for c in CASES)
    nb = sum(V1[f"{d}|{c}"]["uq_abc"]["covers_obs_range"] for d in DEVS for c in CASES)
    L += [f"* Regression: model TTEs within {mx:.2f} % of [PROTO]; additivity (c) model {ad['14 Pro']['model']:+.2f} / "
          f"{ad['13 Pro']['model']:+.2f} % ([PROTO] +0.41), data {ad['14 Pro']['data_cam60']:+.2f} / "
          f"{ad['13 Pro']['data_cam60']:+.2f} % (14 / 13 Pro, camera 60 fps as in the combined run).",
          "* Reading: idle and camera pass (a) on both phones; GPS fails on both (pre-registered for the 14 Pro: the "
          "data imply ≈ +0.03 W, the KW19 prior is 0.27 W; the 13 Pro also fails at −20 %). Camera+GPS passes (a) on "
          "both phones, narrowly on the 13 Pro (−14.7 %); on the 14 Pro the GPS excess is partly offset by the "
          "camera's low power.",
          f"* (b) is strict: the §4 ranges of P_d0 and ΔP_cam lie mostly above their defaults, so the CAL-fixed 90 % "
          f"band sits below the point prediction. It misses the 13 Pro camera range and both camera+GPS ranges, which lie "
          f"above the point prediction, besides GPS. The ABC "
          f"posterior ({V1['uq_settings']['n_abc_accepted']}/{V1['uq_settings']['n_prior']} prior samples within "
          f"±5 % of the 14 Pro idle) covers {nb}/8.",
          f"* V1.4 SOC-shape test, 28 traces: RMS window-slope error D1 {st['D1']['rms_err_pct_mean']:.1f} %, D2 "
          f"{st['D2']['rms_err_pct_mean']:.1f} %; D2 better on {st['n_traces_D2_better']}/28. The data favour the "
          f"energy-fraction SOC (D2), as expected. 14 Pro idle observed "
          f"{'/'.join(f'{x:.1f}' for x in st['obs_14Pro_idle_min_per_pct'])} min per % (spec 31.1/30.8/32.4). D1 window "
          f"OCVs match [PROTO] to {st['proto_check']['max_dev_mV']:.1f} mV.", ""]
    L.append("## V2 LG HG2 (CAL at 25 °C only; f_Q = Q(T)/Q(25 °C); criteria |Δf_Q| ≤ 0.05 (T ≥ 0), 0.10 (−20))")
    L.append("| T, rate | f_Q obs | L2 | X9 | RMSE_V L2 / X9 (mV) |")
    L.append("|---|---|---|---|---|")
    for r in ("0.5C", "2C"):
        for T in (40, 10, 0, -10, -20):
            v = V2[f"{T}|{r}"]
            L.append(f"| {T} °C {r} | {v['fQ_obs']:.3f} | {v['fQ_model']:.3f} | {v['fQ_model_X9']:.3f} | "
                     f"{1000 * v['rmse_V']:.0f} / {1000 * v['rmse_V_X9']:.0f} |")
    ea = V2["Ea_apparent_kJmol"]
    rs = V2["R_sigma_m20C_ohm"]
    ea2 = 8.314 * math.log(rs["model"] / 0.038) / (1 / 253.15 - 1 / 298.15) / 1e3
    L += [f"* As pre-registered, L2 under-predicts the cold loss: R_Σ(−20 °C) = {rs['model']:.3f} Ω (spec 0.165) vs "
          f"apparent {rs['obs_apparent']:.3f} Ω; effective E_a {ea2:.1f} kJ/mol (spec 20.5) vs data {ea['obs']:.1f}. "
          "X9 (KiBaM: c 0.166, κ_K 40 h⁻¹, E_a,K 35 kJ/mol, R0+R1 kept) closes most of the 2C gap. At 0.5C it "
          "is still 0.04–0.12 too high (−20 °C: 0.85 vs 0.73).",
          "* The data also contain a non-kinetic loss: the C/20 ratio is 0.97/0.95/0.93/0.92/0.90 at "
          "40/10/0/−10/−20 °C (sequential-test aging and thermodynamic capacity). The model's C/20 ratio is 1.00. "
          f"Dividing it out (fQ_obs_over_C20) leaves L2 {V2['0|0.5C']['fQ_model'] - V2['0|0.5C']['fQ_obs_over_C20']:+.3f} "
          "at 0 °C 0.5C (diagnostic only).",
          f"* 25 °C RMSE_V (CAL temperature): {1000 * V2['25|0.5C']['rmse_V']:.0f} mV (0.5C), "
          f"{1000 * V2['25|2C']['rmse_V']:.0f} mV (2C). UDDS dynamic check: not run (no drive-cycle file in lg_hg2_raw).", ""]
    L.append("## V3 NASA B0005 (aging -> TTE), V4 Oxford (pulsed load)")
    L += [f"* V3: observed 2 A TTE {V3['obs']['tte_fresh_min']:.1f} -> {V3['obs']['tte_aged_min']:.1f} min "
          f"({V3['obs']['ratio_pct']:+.1f} %); model {V3['model_fresh']['tte_min']:.1f} -> {V3['model_aged']['tte_min']:.1f} "
          f"min ({V3['model_ratio_pct']:+.1f} %; Q only {V3['model_ratio_Q_only_pct']:+.1f} %, spec −28.5 %; R only "
          f"{V3['model_ratio_R_only_pct']:+.2f} %). Both runs end on the `charge` event (z = 0 before 2.7 V; "
          "core_issues C2), so V3 confirms Q↓ only, not R↑ through the cut-off. RMSE_V "
          f"{1000 * V3['rmse_V_fresh']:.0f}/{1000 * V3['rmse_V_aged']:.0f} mV. 4 °C check (B0055): model Q ratio "
          f"{V3['temperature_sanity']['B0055']['ratio_model_to_B0005_24C']:.2f} vs {V3['temperature_sanity']['B0055']['ratio_obs_to_B0005_24C']:.2f} "
          "(no criterion; clearly off, same cause).",
          f"* V4: RMSE_V {1000 * V4['rmse_V']:.0f} mV for the R·Q-scaled cell ECM (R0·Q 89 mΩ·Ah). L1 gives "
          f"{1000 * V4['variants']['L1_cellRQ']['rmse_V']:.0f} mV, so the RC branch helps. But OCV-only gives "
          f"{1000 * V4['variants']['OCV_only']['rmse_V']:.0f} mV: the Kokam power cell's measured 1-s step resistance "
          f"is {1000 * V4['obs']['step_resistance_1s_median_ohm']:.0f} mΩ vs the model's "
          f"{1000 * V4['variants']['L2_cellRQ']['R0_40C_ohm']:.0f} mΩ. R·Q scaling (A12) does not transfer from energy "
          "cells to a power cell.", ""]
    if V5:
        L.append("## V5 PyBaMM DFN (Ai2020, isothermal; ΔTTE % vs DFN; criterion ≤ 3 % at ≥ 0 °C)")
        L.append("| case | DFN h | L1 | L2 | X5 | X9 |")
        L.append("|---|---|---|---|---|---|")
        for k, v in V5["cases"].items():
            if "dtte_pct" in v:
                d = v["dtte_pct"]
                L.append(f"| {k[3:]} | {v['tte_dfn_h']:.3f} | {d['L1']:+.1f} | {d['L2']:+.1f} | {d['X5']:+.1f} | "
                         f"{d['X9']:+.1f} |")
        pr = [v for k, v in V5["cases"].items() if k.startswith("pulse")]
        L += [f"* At 25 °C every structure is within 2.4 % (X9 within 0.7 %). At 0 °C, 0.5C only X9 passes. At 1C "
              "(0 °C) and ≥ 0.5C (−10 °C) the DFN delivers only 1.67 / 1.21 of 2.47 Ah, a transport-limited collapse "
              "(mechanism not isolated here) that no ECM represents: structural error +30 % to +360 %. Sustained phone "
              "loads (≤ 0.21C, gaming) lie between the passing 0.1C and failing 0.5C cases at 0 °C.",
              f"* 10 s pulses (9 cases): RMSE_V L2 {min(v['rmse_V']['L2'] for v in pr) * 1000:.0f}–"
              f"{max(v['rmse_V']['L2'] for v in pr) * 1000:.0f} mV, X5 marginally lower, L1/X9 "
              f"{min(v['rmse_V']['L1'] for v in pr) * 1000:.0f}–{max(v['rmse_V']['L1'] for v in pr) * 1000:.0f} mV. "
              "R1 seen instantly overshoots the pulse drop; the 25/0 °C pulses are near in-sample for the [RUN] ECM fits. "
              f"Local L2 = core LSODA to {abs(V5['core_L2_check']['rel_diff']):.0e}.", ""]
    L.append("## V6 published band (±25 %), V7 cold phenomenology")
    L += [f"* V6: {V6['n_pass']}/{V6['n_cases']} pass ([PROTO] 11/16; model hours within 0.03 % of [PROTO]). Fails: " +
          ", ".join(f"{k} {V6[k]['err_pct']:+.0f} %" for k in V6 if isinstance(V6[k], dict) and V6[k].get("pass") is False
                    and "err_pct" in V6[k]) + f". Every iPhone case passes except BT audio; X14 (DSP audio) gives "
          f"{V6['14P_audio_BT']['X14_err_pct']:+.1f} %. Android screen-on floors are 1.6–2.9× ours (not tuned).",
          f"* V7 (pre-registered config: Q 2.7 Ah, R_Σ(20 °C) 0.25 Ω, ΔP_b 3 W, P̄ = E_u/159 min = "
          f"{V7['config']['P_mean_W']:.2f} W): 0 °C TTE loss {V7['tte_loss_0C_pct']:+.0f} % (LEE20 −15 %), s_off "
          + "/".join(f"{V7['brownout_soc_at_le0C_pct'][T]:.0f}" for T in ("0", "-5", "-10", "-15"))
          + f" % at 0/−5/−10/−15 °C (want 15–35 %), recovered {V7['recovery_after_warming']['extra_mAh']:.0f} mAh "
          "(LEE20 330): all three fail. Even at 25 °C the phone browns out at "
          f"{V7['runs']['25']['s_off_pct']:.0f} %.",
          "* V7 sensitivity (not fitted): with LEE20's own shut-off resistance 0.13 Ω (C5) the 0 °C loss is "
          f"{V7['sensitivity']['R_sigma_0.13ohm_LEE20_C5']['tte_loss_0C_pct']:+.0f} %, s_off(−5 °C) "
          f"{V7['sensitivity']['R_sigma_0.13ohm_LEE20_C5']['s_off_pct']['-5']:.0f} %, recovery "
          f"{V7['sensitivity']['R_sigma_0.13ohm_LEE20_C5']['recovered_mAh']:.0f} mAh, all inside the bands. The "
          "outcome hinges on the read-off phone resistance; the mechanisms (brown-out, recovery on warming) are reproduced.",
          ]
    return "\n".join(L) + "\n"


# =========================================================================================== main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-uq", type=int, default=512)
    ap.add_argument("--n-prior", type=int, default=4096)
    ap.add_argument("--skip-pybamm", action="store_true")
    a = ap.parse_args()
    T0 = _t()
    R, F = {}, {}
    print("V1 Bonn ...", flush=True)
    R["V1"], F["V1"] = v1_bonn(a.n_uq, a.n_prior)
    print("V2 HG2 ...", flush=True)
    R["V2"], F["V2"] = v2_hg2()
    print("V3 NASA ...", flush=True)
    R["V3"] = v3_nasa()
    print("V4 Oxford ...", flush=True)
    R["V4"], F["V4"] = v4_oxford()
    if not a.skip_pybamm:
        print("V5 PyBaMM ...", flush=True)
        R["V5"], F["V5"] = v5_pybamm()
    print("V6 band ...", flush=True)
    R["V6"] = v6_band()
    print("V7 cold ...", flush=True)
    R["V7"] = v7_cold()
    R["criteria"] = criteria(R)
    out = dict(meta=meta("run_validation.py", seed=SEED, solver="LSODA (phonebatt.run, rtol 1e-7) for point "
                         "predictions; phonebatt.run_batch RK4 for V1 UQ; local cell drivers (LSODA) for V2-V5",
                         h_s=dict(V1_uq=UQ_H, V4_rk4=0.25), n_uq=a.n_uq,
                         n_prior=a.n_prior, runtime_s=None), **R)
    if "V5" in out:  # §14 layout: V5.{case}.dtte_pct (cases at the top level of V5)
        v5 = out["V5"]
        out["V5"] = dict(**v5["cases"], summary=v5["summary"], core_L2_check=v5["core_L2_check"],
                         settings=v5["settings"], _runtime_s=v5["_runtime_s"])
    print("figures ...", flush=True)
    fig_bonn(R["V1"], F["V1"])
    fig_hg2(R["V2"], F["V2"])
    fig_band(R["V6"])
    if "V5" in R:
        fig_pybamm(R["V5"], F["V5"])
    write_table(R)
    out["meta"]["runtime_s"] = _t() - T0
    write_json("validation.json", out)
    from basics_common import RESULTS
    (RESULTS / "validation_summary.md").write_text(summary_md(R))
    print(f"done in {_t() - T0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
