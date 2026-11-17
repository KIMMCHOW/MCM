"""Validation data loaders and pre-registered validation cases (FINAL_SPEC §8, §15.1).

Data are read from ``docs/research/data/validation/`` (licences in its README.md).  Nothing here tunes model
parameters except the single [CAL] step of V1.1 (``bonn_calibrate_awake``) and the cell-only 25 °C CAL of V2
(``hg2_calibration``), both declared in the spec.
"""
from __future__ import annotations

import csv
import json
import math
from dataclasses import replace
from functools import lru_cache
from pathlib import Path

import numpy as np

from .params import (RP63, Params, galaxy_s24, hg2_cell, iphone13pro, iphone14pro, iphone14promax, nits_to_b,
                     pixel8, pixel8a)
from .usage import BONN_CASES, PRIMITIVES, Schedule, UsageState, s12_schedule

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "docs" / "research" / "data" / "validation"


# ------------------------------------------------------------------------------------------- V1 Bonn iPhones
BONN_EXPERIMENTS = {"idle": "07-09-2023-Idle", "cam": "09-11-2023-CAM", "gps": "02-10-2023-GPS",
                    "cxl": "23-11-2023-CXL"}  # the runs quoted in FINAL_SPEC §8 V1 (validation_data a.2)
BONN_DEVICES = {"14 Pro": iphone14pro, "13 Pro": iphone13pro}


@lru_cache(maxsize=1)
def load_bonn_summary(path: str | None = None) -> dict:
    with open(Path(path) if path else DATA_DIR / "bonn_iphone_summary.json") as f:
        return json.load(f)


def bonn_observations() -> dict:
    """Observed runtimes to the last logged 1 % for the V1 cases: {(device, case): {obs_h, min_h, max_h, n,
    all_medians_h}} (device '14 Pro' | '13 Pro'; case idle | cam | gps | cxl)."""
    d = load_bonn_summary()["groups"]
    sensor = {"idle": "Idle", "cam": "CAM", "gps": "GPS", "cxl": "CXL"}
    out = {}
    for dev in BONN_DEVICES:
        for case, exp in BONN_EXPERIMENTS.items():
            g = d[f"{exp}|{dev} (inferred)"]
            others = [v["runtime_h_median"] for k, v in d.items()
                      if k.endswith(f"|{dev} (inferred)") and k.split("|")[0].endswith(sensor[case])]
            out[(dev, case)] = dict(obs_h=g["runtime_h_median"], min_h=g["runtime_h_min"], max_h=g["runtime_h_max"],
                                    n=g["n_devices"], experiment=exp, all_medians_h=sorted(others))
    return out


def bonn_soc_curves():
    """Bonn displayed-level curves (hours to reach each displayed %), as a pandas DataFrame."""
    import pandas as pd
    return pd.read_csv(DATA_DIR / "bonn_iphone_summary_soc_curves.csv")


def bonn_predictions(P_awake: float | None = None, solver: str = "lsoda", **run_kw) -> dict:
    """Model predictions for the eight V1 cases (25 °C, airplane mode, black OLED, no bursts)."""
    from .model import run
    from .power import demand_power
    obs = bonn_observations()
    out = {}
    for dev, ctor in BONN_DEVICES.items():
        p = ctor() if P_awake is None else ctor(P_awake=P_awake)
        for case, u in BONN_CASES.items():
            r = run(Schedule.constant(u), 25.0, p=p, solver=solver, **run_kw)
            o = obs[(dev, case)]["obs_h"]
            out[(dev, case)] = dict(model_h=r.tte_h, obs_h=o, err_pct=100 * (r.tte_h / o - 1), P_W=float(
                demand_power(u, p)), event=r.event)
    return out


def bonn_additivity(pred: dict, dev: str = "14 Pro") -> float:
    """Model additivity residual: (1/TTE_cxl)/(1/TTE_cam + 1/TTE_gps − 1/TTE_idle) − 1 (V1 criterion c)."""
    t = {c: pred[(dev, c)]["model_h"] for c in BONN_CASES}
    return (1 / t["cxl"]) / (1 / t["cam"] + 1 / t["gps"] - 1 / t["idle"]) - 1


def bonn_calibrate_awake(target_h: float | None = None) -> float:
    """V1.1 [CAL]: ΔP_awake such that the iPhone 14 Pro black-screen idle TTE equals the observed median."""
    from scipy.optimize import brentq
    from .model import run
    tgt = bonn_observations()[("14 Pro", "idle")]["obs_h"] if target_h is None else target_h
    sch = Schedule.constant(BONN_CASES["idle"])
    return brentq(lambda aw: run(sch, 25.0, p=iphone14pro(P_awake=aw)).tte_h - tgt, 0.03, 0.2, xtol=1e-6)


# ------------------------------------------------------------------------------------------- V6 published band
def _v6_cases() -> list[dict]:
    so, vw = PRIMITIVES["SO"], PRIMITIVES["VW"]
    b150 = float(nits_to_b(150.0))
    local = replace(vw, name="VW-local", lam=0.0, d_dl=0.0)
    local150 = replace(local, b=b150, w=0, rho="none")
    idle_min = UsageState("idle-minbright", s_scr=1, b=0.0, APL=0.65, f_r=60, w=0, c="lte", rho="none")
    return [
        dict(id="14P_local_video", src="APL14P", p=iphone14pro(), u=replace(local, bt=1, au=0), obs_h=23.0),
        dict(id="14P_stream_video", src="APL14P", p=iphone14pro(), u=replace(vw, bt=1, au=0), obs_h=20.0),
        dict(id="14P_audio_BT", src="APL14P", p=iphone14pro(), u=PRIMITIVES["AU"], obs_h=75.0),
        dict(id="14PM_wifi_web_150", src="NBC22", p=iphone14promax(), u=replace(so, b=b150), obs_h=15.8),
        dict(id="14PM_local_video_150", src="NBC22", p=iphone14promax(), u=local150, obs_h=28.13),
        dict(id="14PM_idle_minbright", src="NBC22", p=iphone14promax(), u=idle_min, obs_h=57.9),
        dict(id="14P_cell_web_150", src="LMAG22", p=iphone14pro(), u=replace(PRIMITIVES["SO-L"], b=b150),
             obs_h=10.37),
        dict(id="Pixel8_wifi_web", src="NBC23", p=pixel8(), u=replace(so, b=b150), obs_h=13.57),
        dict(id="Pixel8_local_video", src="NBC23", p=pixel8(), u=local150, obs_h=22.72),
        dict(id="Pixel8_idle_minbright", src="NBC23", p=pixel8(), u=idle_min, obs_h=24.88),
        dict(id="S24_wifi_web", src="NBC24", p=galaxy_s24(), u=replace(so, b=b150), obs_h=14.67),
        dict(id="S24_local_video", src="NBC24", p=galaxy_s24(), u=local150, obs_h=20.45),
        dict(id="S24_idle_minbright", src="NBC24", p=galaxy_s24(), u=idle_min, obs_h=39.97),
        dict(id="Pixel8a_S12_day", src="GPX", p=pixel8a(), u=s12_schedule(), obs_h=31.0, bound="lower"),
    ]


V6_TOL = 0.25


def v6_cases() -> list[dict]:
    """V6 published-specification cases (§8): dicts with id, src, p (device Params), u (input), obs_h."""
    return _v6_cases()


def v6_predictions(**run_kw) -> dict:
    """Model TTE for each V6 case with ε = TTE_model/TTE_obs − 1 and the ±25 % pass flag."""
    from .model import run
    out = {}
    for c in v6_cases():
        r = run(c["u"], 25.0, p=c["p"], **run_kw)
        e = r.tte_h / c["obs_h"] - 1
        ok = (r.tte_h >= c["obs_h"]) if c.get("bound") == "lower" else abs(e) <= V6_TOL
        out[c["id"]] = dict(model_h=r.tte_h, obs_h=c["obs_h"], err_pct=100 * e, pass_=bool(ok), src=c["src"],
                            P_mean_W=r.P_mean_W)
    return out


# ------------------------------------------------------------------------------------------- V2 LG HG2 cell
HG2_TESTS = {40: "555", 25: "549", 10: "575", 0: "585", -10: "593", -20: "607"}
HG2_KINDS = ("C20DisCh", "Dis_0p5C", "Dis_2C")


def _hms(s: str) -> float:
    h, m, sec = s.strip().split(":")
    return 3600.0 * float(h) + 60.0 * float(m) + float(sec)


def load_hg2(T_C: int, kind: str = "Dis_0p5C") -> dict:
    """Longest contiguous discharge (DCH) segment of one LG HG2 raw CSV (KOL20, CC BY 4.0).

    Returns arrays: t_s (from the step start; the logger's first row is ≈ 60 s into the step), V, I_A (> 0
    discharge), T_C, Ah (integrated from t = 0 with the first logged current held over the first interval).
    """
    f = DATA_DIR / "lg_hg2_raw" / f"{HG2_TESTS[int(T_C)]}_{kind}.csv"
    rows = []
    with open(f, newline="") as fh:
        it = csv.reader(fh)
        for line in it:
            if line and line[0].startswith("Time Stamp"):
                header = [c.strip() for c in line]
                next(it)  # units row
                break
        iS, iSt, iT, iV, iI, iTe = (header.index(k) for k in ("Step", "Status", "Step Time", "Voltage", "Current",
                                                               "Temperature"))
        for line in it:
            if len(line) <= iTe or not line[iV].strip():
                continue
            rows.append((int(line[iS]), line[iSt].strip(), _hms(line[iT]), float(line[iV]), float(line[iI]),
                         float(line[iTe])))
    segs, cur = [], []
    for r in rows:
        if r[1] == "DCH" and (not cur or r[0] == cur[-1][0]):
            cur.append(r)
        else:
            if cur:
                segs.append(cur)
            cur = [r] if r[1] == "DCH" else []
    if cur:
        segs.append(cur)
    seg = max(segs, key=len)
    a = np.array([r[2:] for r in seg], dtype=float)
    t, V, I, T = a[:, 0], a[:, 1], -a[:, 2], a[:, 3]
    dt = np.diff(np.concatenate([[0.0], t]))
    Ah = np.cumsum(I * dt) / 3600.0
    return dict(t_s=t, V=V, I_A=I, T_C=T, Ah=Ah, file=f.name)


def hg2_calibration(R_sigma: float = 0.038, n_grid: int = 201) -> dict:
    """V2 [CAL] at 25 °C only (§4.7): NMC OCV from the C/20 discharge (V + I·R_Σ), Q from its capacity, and the
    cell C_th, R_th fitted to the 25 °C 2C temperature rise (heat = I·(U − V), measured overpotential).

    Returns dict(ocv_table=(z, U), Q_Ah, C_th, R_th, rmse_T_K) for ``params.hg2_cell``.
    """
    from scipy.optimize import least_squares
    c20 = load_hg2(25, "C20DisCh")
    Q = float(c20["Ah"][-1])
    z = 1.0 - c20["Ah"] / Q
    U = c20["V"] + c20["I_A"] * R_sigma
    zg = np.linspace(0.0, 1.0, n_grid)
    o = np.argsort(z)
    Ug = np.interp(zg, np.concatenate([[0.0], z[o]]), np.concatenate([[U[o][0]], U[o]]))
    table = (zg, Ug)
    d2 = load_hg2(25, "Dis_2C")
    zt = 1.0 - d2["Ah"] / Q
    q = d2["I_A"] * (np.interp(np.clip(zt, 0, 1), zg, Ug) - d2["V"])
    t = d2["t_s"]
    Tm = d2["T_C"]
    Ta = float(Tm[0])

    def sim(par):
        C, R = par
        T = np.empty_like(Tm)
        T[0] = Tm[0]
        for i in range(1, len(t)):
            dt = t[i] - t[i - 1]
            T[i] = T[i - 1] + dt * (q[i - 1] - (T[i - 1] - Ta) / R) / C
        return T

    fit = least_squares(lambda par: sim(par) - Tm, x0=[50.0, 5.0], bounds=([5.0, 0.5], [500.0, 50.0]))
    rmse = float(np.sqrt(np.mean(fit.fun ** 2)))
    return dict(ocv_table=table, Q_Ah=Q, C_th=float(fit.x[0]), R_th=float(fit.x[1]), rmse_T_K=rmse)


def hg2_params(cal: dict | None = None) -> Params:
    """HG2 cell parameters with the V2 25 °C calibration applied."""
    cal = hg2_calibration() if cal is None else cal
    return hg2_cell(ocv_table=cal["ocv_table"], C_th=cal["C_th"], R_th=cal["R_th"], Q_Ah=cal["Q_Ah"])


def load_hg2_summary() -> dict:
    with open(DATA_DIR / "lg_hg2_summary.json") as f:
        return json.load(f)


# ------------------------------------------------------------------------------------------- V3 NASA, V4 Oxford
def _read_csv(path: Path) -> dict:
    with open(path, newline="") as fh:
        r = csv.reader(fh)
        head = next(r)
        a = np.array([[float(x) for x in row] for row in r if row], dtype=float)
    return {h: a[:, i] for i, h in enumerate(head)}


def nasa_index() -> dict:
    with open(DATA_DIR / "nasa_curves_index.json") as f:
        return json.load(f)


def load_nasa(name: str = "nasa_B0005_24C_2A_CC_fresh.csv") -> dict:
    """NASA PCoE discharge curve (t_s, V, I_A > 0 discharge, T_cell_C), decimated to 10 s."""
    return _read_csv(DATA_DIR / name)


def nasa_summary() -> list:
    with open(DATA_DIR / "nasa_pcoe_summary.json") as f:
        return json.load(f)


def load_oxford() -> dict:
    """Oxford ExampleDC cell 1 Artemis discharge at 40 °C (t_s, V, I_A > 0 discharge, T_C), 1 s logging."""
    return _read_csv(DATA_DIR / "oxford_ExampleDC_C1_artemis_discharge.csv")


def load_device_benchmarks() -> dict:
    with open(DATA_DIR / "device_benchmarks.json") as f:
        return json.load(f)


__all__ = ["DATA_DIR", "BONN_EXPERIMENTS", "BONN_DEVICES", "load_bonn_summary", "bonn_observations",
           "bonn_soc_curves", "bonn_predictions", "bonn_additivity", "bonn_calibrate_awake", "v6_cases",
           "v6_predictions", "V6_TOL", "HG2_TESTS", "load_hg2", "hg2_calibration", "hg2_params",
           "load_hg2_summary", "nasa_index", "load_nasa", "nasa_summary", "load_oxford", "load_device_benchmarks"]
