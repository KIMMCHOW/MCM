"""Summarise LG HG2 (3 Ah NMC 18650) constant-current discharges at six ambient temperatures.

Usage:
    python3 -I validation_lg_hg2_summary.py <dir_with_csvs> <out.json>

Source: P. Kollmeyer, C. Vidal, M. Naguib, M. Skells (2020), "LG 18650HG2 Li-ion Battery Data and Example Deep
Neural Network xEV SOC Estimator Script", Mendeley Data v3, doi:10.17632/cp3473x7xv.3, CC BY 4.0.
Files used: <test>_C20DisCh.csv, <test>_Dis_0p5C.csv, <test>_Dis_2C.csv (Digatron export).
Test-number -> ambient temperature mapping is taken from the .mat file names in the same dataset.
Capacity and energy are integrated from the logged current/voltage over all 'DCH' rows (trapezoid on Prog Time);
the first, unlogged minute of each step is filled with the first logged current and voltage (Step Time column).
The first logged voltage is therefore ~1 min into the discharge (V_start below is that value, not the rest OCV).
"""
import sys, os, glob, json, csv
import numpy as np

TEMP = {"549": 25, "555": 40, "575": 10, "585": 0, "593": -10, "607": -20}


def hms(s):
    h, m, sec = s.split(":")
    return int(h) * 3600 + int(m) * 60 + float(sec)


def load(path):
    rows = []
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        r = csv.reader(f)
        header = None
        for line in r:
            if header is None:
                if line and line[0].strip() == "Time Stamp":
                    header = [c.strip() for c in line]
                continue
            if not line or not line[0] or line[0].startswith(","):
                continue
            if len(line) < len(header) - 1:
                continue
            try:
                d = dict(zip(header, line))
                rows.append((d["Status"].strip(), hms(d["Prog Time"]), float(d["Voltage"]), float(d["Current"]),
                             float(d["Temperature"]), hms(d["Step Time"])))
            except (ValueError, KeyError):
                continue
    return rows


def summarise(rows):
    dch = [x for x in rows if "DCH" in x[0]]
    t = np.array([x[1] for x in dch]); v = np.array([x[2] for x in dch]); i = -np.array([x[3] for x in dch])
    T = np.array([x[4] for x in dch])
    # split into contiguous segments (gap > 120 s) and keep the longest one
    gaps = np.nonzero(np.diff(t) > 120)[0]
    segs = np.split(np.arange(t.size), gaps + 1)
    s = max(segs, key=lambda idx: t[idx[-1]] - t[idx[0]])
    t, v, i, T = t[s], v[s], i[s], T[s]
    st0 = [x[5] for x in dch][s[0]]  # step time of first logged sample (logger writes the first row after ~60 s)
    # prepend the step start (t - st0) with the first logged current and voltage, so the first minute is counted
    t = np.concatenate([[t[0] - st0], t]); v = np.concatenate([[v[0]], v]); i = np.concatenate([[i[0]], i])
    T = np.concatenate([[T[0]], T])
    q = float(np.trapezoid(i, t) / 3600); e = float(np.trapezoid(i * v, t) / 3600)
    return dict(n=int(t.size), duration_min=float((t[-1] - t[0]) / 60), I_mean_A=float(np.mean(i)), Q_Ah=q, E_Wh=e,
                V_start=float(v[0]), V_end=float(v[-1]), V_mean=float(e / q), T_start_C=float(T[0]),
                T_max_C=float(T.max()), dt_median_s=float(np.median(np.diff(t))))


def main():
    root, out = sys.argv[1], sys.argv[2]
    res = []
    for p in sorted(glob.glob(os.path.join(root, "*.csv"))):
        base = os.path.basename(p)
        test, kind = base.split("_", 1)
        kind = kind[:-4]
        if test not in TEMP:
            continue
        s = summarise(load(p))
        s.update(file=base, test=test, T_amb_C=TEMP[test], kind=kind)
        res.append(s)
    res.sort(key=lambda r: (r["kind"], -r["T_amb_C"]))
    # apparent DC resistance from the mean-voltage drop between the C/20 and 0.5C discharges at each temperature,
    # and an Arrhenius fit R(T) = R(25C) * exp[(Ea/Rg) (1/T - 1/298.15 K)] through the 25 C point (T <= 25 C)
    X = {(r["kind"], r["T_amb_C"]): r for r in res}
    rapp = {}
    for T in sorted({r["T_amb_C"] for r in res}):
        if ("C20DisCh", T) not in X or ("Dis_0p5C", T) not in X:
            continue
        a, b = X[("C20DisCh", T)], X[("Dis_0p5C", T)]
        rapp[T] = (a["V_mean"] - b["V_mean"]) / (b["I_mean_A"] - a["I_mean_A"])
    xs = [1 / (T + 273.15) - 1 / 298.15 for T in rapp if T <= 25]
    ys = [np.log(rapp[T] / rapp[25]) for T in rapp if T <= 25]
    ea_r = float(np.dot(xs, ys) / np.dot(xs, xs))
    arr = dict(R_apparent_ohm={str(k): float(v) for k, v in rapp.items()}, EaR_K=ea_r, Ea_kJ_per_mol=ea_r * 8.314 / 1000,
               note="apparent R includes ohmic + polarisation + OCV-window effects; DERIVED, not an EIS value")
    with open(out, "w") as g:
        json.dump(dict(source="Kollmeyer et al. 2020, Mendeley Data doi:10.17632/cp3473x7xv.3 (CC BY 4.0)",
                       note="Q, E integrated over the longest contiguous DCH segment", results=res,
                       arrhenius_apparent_resistance=arr), g, indent=1)
    print("apparent R (ohm):", {k: round(v, 4) for k, v in rapp.items()}, f"Ea/Rg = {ea_r:.0f} K, Ea = {ea_r*8.314/1000:.1f} kJ/mol")
    for r in res:
        print(f"{r['kind']:9s} T={r['T_amb_C']:4d}C I={r['I_mean_A']:.3f}A Q={r['Q_Ah']:.3f}Ah E={r['E_Wh']:.2f}Wh "
              f"t={r['duration_min']:.1f}min V {r['V_start']:.3f}->{r['V_end']:.3f} Vbar={r['V_mean']:.3f} "
              f"Tcell {r['T_start_C']:.1f}->{r['T_max_C']:.1f}C dt={r['dt_median_s']:.1f}s n={r['n']}")


if __name__ == "__main__":
    main()
