"""Extract a few decimated discharge curves from the NASA PCoE Battery Data Set for model validation.

Usage:
    python3 -I validation_nasa_extract_curves.py <dir_with_extracted_B00xx.mat_tree> <out_dir>

Writes one CSV per selected (cell, cycle): columns t_s, V, I_A (discharge positive), T_cell_C, plus a JSON index.
Cycle choice: the k-th discharge cycle (0-based, counting only runs > 5 min that start at V >= 4.1 V) at the
requested ambient temperature.
Source: B. Saha and K. Goebel (2007), "Battery Data Set", NASA Prognostics Data Repository, NASA Ames.
"""
import sys, os, glob, json
import numpy as np
from scipy.io import loadmat

SELECT = [  # (cell, ambient_C, k-th full-charge discharge, tag)
    ("B0005", 24, 0, "24C_2A_CC_fresh"),
    ("B0005", 24, -1, "24C_2A_CC_aged"),
    ("B0029", 43, 0, "43C_4A_CC_fresh"),
    ("B0047", 4, 0, "4C_1A_CC_fresh"),
    ("B0055", 4, 0, "4C_2A_CC_fresh"),
    ("B0025", 24, 0, "24C_4A_square_0p05Hz_fresh"),
]
DT = 10.0  # s, decimation step (square-wave curve kept at native resolution)


def f(x):
    return np.real(np.asarray(x)).astype(float).ravel()


def main():
    root, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    paths = {os.path.splitext(os.path.basename(p))[0]: p
             for p in sorted(glob.glob(os.path.join(root, "**", "B00*.mat"), recursive=True))}
    index = []
    for cell, amb, k, tag in SELECT:
        m = loadmat(paths[cell], squeeze_me=False, struct_as_record=True)
        cyc = m[cell][0, 0]["cycle"][0]
        runs = []
        for c in cyc:
            if str(c["type"][0]) != "discharge":
                continue
            if round(float(np.asarray(c["ambient_temperature"]).ravel()[0])) != amb:
                continue
            d = c["data"][0, 0]
            t, v, i, T = f(d["Time"]), f(d["Voltage_measured"]), f(d["Current_measured"]), f(d["Temperature_measured"])
            if t.size < 5 or t[-1] - t[0] < 300 or v[0] < 4.1:
                continue
            runs.append((t - t[0], v, -i, T))
        t, v, i, T = runs[k]
        if "square" not in tag:
            grid = np.arange(0, t[-1], DT)
            keep = np.unique(np.searchsorted(t, grid))
            keep = keep[keep < t.size]
            t, v, i, T = t[keep], v[keep], i[keep], T[keep]
        fn = f"nasa_{cell}_{tag}.csv"
        with open(os.path.join(out, fn), "w") as g:
            g.write("t_s,V,I_A,T_cell_C\n")
            for a, b, c_, e in zip(t, v, i, T):
                g.write(f"{a:.1f},{b:.4f},{c_:.4f},{e:.2f}\n")
        on = np.nonzero(i > 0.1)[0]
        q = float(np.trapezoid(np.clip(i, 0, None), t) / 3600)
        index.append(dict(file=fn, cell=cell, ambient_C=amb, cycle_k=k, n_full_runs=len(runs),
                          t_load_end_s=float(t[on[-1]]), Q_Ah=q, V0=float(v[0]), Vmin=float(v.min()),
                          T0=float(T[0]), Tmax=float(T.max())))
        print(index[-1])
    with open(os.path.join(out, "nasa_curves_index.json"), "w") as g:
        json.dump(dict(source="Saha & Goebel 2007, NASA PCoE Battery Data Set; "
                              "https://phm-datasets.s3.amazonaws.com/NASA/5.+Battery+Data+Set.zip",
                       sign="I_A > 0 = discharge", decimation_s=DT, curves=index), g, indent=1)


if __name__ == "__main__":
    main()
