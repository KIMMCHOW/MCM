"""Summarise the NASA Ames PCoE 'Battery Data Set' (Saha & Goebel 2007) for model validation.

Usage:
    python3 -I validation_nasa_summary.py <dir_with_extracted_B00xx.mat_tree> <out.json>

Data source (not redistributed here; download it yourself):
    https://phm-datasets.s3.amazonaws.com/NASA/5.+Battery+Data+Set.zip   (209,708,670 bytes)
For each cell this reports the ambient temperature, the discharge-cycle count, and for the first and last
discharge cycles: the delivered capacity (A h, integrated from Current_measured down to the run's cutoff), the
duration, the mean current, the initial and final terminal voltage, the cell temperature rise, and the mean
energy-weighted voltage. It also gives the first and last fitted EIS resistances (Re, Rct).
"""
import sys, os, json, glob
import numpy as np
from scipy.io import loadmat


def _f(x):
    return np.asarray(x, dtype=float).ravel()


def summarise_cell(path):
    name = os.path.splitext(os.path.basename(path))[0]
    m = loadmat(path, squeeze_me=False, struct_as_record=True)
    cyc = m[name][0, 0]["cycle"][0]
    dis, imp = [], []
    for c in cyc:
        typ = str(c["type"][0])
        amb = float(np.asarray(c["ambient_temperature"]).ravel()[0])
        d = c["data"][0, 0]
        if typ == "discharge":
            t = _f(d["Time"]); v = _f(d["Voltage_measured"]); i = _f(d["Current_measured"])
            T = _f(d["Temperature_measured"])
            if t.size < 5:
                continue
            # discharge current is reported negative; integrate |I|
            q_int = float(np.trapezoid(np.abs(i), t) / 3600.0)
            e_int = float(np.trapezoid(np.abs(i) * v, t) / 3600.0)
            cap = float(_f(d["Capacity"])[0]) if "Capacity" in d.dtype.names and _f(d["Capacity"]).size else float("nan")
            on = np.nonzero(np.abs(i) > 0.1)[0]
            t_load = float(t[on[-1]] - t[0]) if on.size else float("nan")  # time until load removed (= time to cutoff)
            dis.append(dict(amb_C=amb, capacity_Ah_reported=cap, q_int_Ah=q_int, energy_Wh=e_int,
                            duration_s=float(t[-1] - t[0]), t_load_s=t_load, mean_abs_I_A=float(np.mean(np.abs(i))),
                            V_start=float(v[0]), V_end_min=float(np.min(v)),
                            V_mean_energy=float(e_int / q_int) if q_int > 0 else float("nan"),
                            T_start_C=float(T[0]), T_max_C=float(np.max(T))))
        elif typ == "impedance":
            try:
                re_ = float(_f(d["Re"])[0]); rct = float(_f(d["Rct"])[0])
                imp.append(dict(amb_C=amb, Re_ohm=re_, Rct_ohm=rct))
            except Exception:
                pass
    # drop obviously aborted runs (< 5 min) for first/last statistics
    good = [x for x in dis if x["duration_s"] > 300]
    # per-ambient statistics over cycles that started from a full charge (V_start >= 4.1 V)
    by_amb = {}
    for x in good:
        if x["V_start"] < 4.1:
            continue
        by_amb.setdefault(str(int(round(x["amb_C"]))), []).append(x)
    amb_stats = {}
    for k, xs in by_amb.items():
        first5 = xs[:5]
        amb_stats[k] = dict(n_full=len(xs),
                            Q_first5_median_Ah=float(np.median([x["q_int_Ah"] for x in first5])),
                            E_first5_median_Wh=float(np.median([x["energy_Wh"] for x in first5])),
                            t_load_first5_median_s=float(np.median([x["t_load_s"] for x in first5])),
                            I_first5_median_A=float(np.median([x["mean_abs_I_A"] for x in first5])),
                            Vbar_first5_median=float(np.median([x["V_mean_energy"] for x in first5])),
                            dT_first5_median_C=float(np.median([x["T_max_C"] - x["T_start_C"] for x in first5])),
                            Q_max_Ah=float(max(x["q_int_Ah"] for x in xs)),
                            Q_last_Ah=float(xs[-1]["q_int_Ah"]))
    out = dict(cell=name, n_discharge=len(dis), n_discharge_gt5min=len(good),
               ambient_C=sorted({round(x["amb_C"]) for x in dis}),
               first=good[0] if good else None, last=good[-1] if good else None,
               q_int_Ah_max=max((x["q_int_Ah"] for x in good), default=None),
               impedance_first=imp[0] if imp else None, impedance_last=imp[-1] if imp else None,
               n_impedance=len(imp), by_ambient_full_charge=amb_stats)
    return out


def main():
    root, out = sys.argv[1], sys.argv[2]
    res = []
    for p in sorted(glob.glob(os.path.join(root, "**", "B00*.mat"), recursive=True)):
        try:
            res.append(summarise_cell(p))
        except Exception as e:  # keep going on odd files
            res.append(dict(cell=os.path.basename(p), error=repr(e)))
    # de-duplicate (some cells appear in two inner zips)
    seen, uniq = set(), []
    for r in res:
        if r["cell"] in seen:
            continue
        seen.add(r["cell"]); uniq.append(r)
    with open(out, "w") as g:
        json.dump(uniq, g, indent=1)
    for r in uniq:
        if "error" in r:
            print(r["cell"], "ERROR", r["error"]); continue
        for k, a in r.get("by_ambient_full_charge", {}).items():
            print(f"{r['cell']} T_amb={k:>2}C nFull={a['n_full']:3d} I~{a['I_first5_median_A']:.2f}A "
                  f"Q5={a['Q_first5_median_Ah']:.3f}Ah E5={a['E_first5_median_Wh']:.2f}Wh "
                  f"t_cut5={a['t_load_first5_median_s']/60:.1f}min Vbar={a['Vbar_first5_median']:.3f}V "
                  f"dT={a['dT_first5_median_C']:.1f}C Qmax={a['Q_max_Ah']:.3f} Qlast={a['Q_last_Ah']:.3f}")
        fi, li = r["impedance_first"], r["impedance_last"]
        if fi and li:
            print(f"     EIS n={r['n_impedance']} Re {fi['Re_ohm']:.4f}->{li['Re_ohm']:.4f} ohm, Rct {fi['Rct_ohm']:.4f}->{li['Rct_ohm']:.4f} ohm")


if __name__ == "__main__":
    main()
