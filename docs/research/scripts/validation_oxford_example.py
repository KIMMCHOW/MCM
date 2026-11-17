"""Summarise and export the Oxford Battery Degradation Dataset 1 example drive-cycle file (ExampleDC_C1.mat).

Usage:
    python3 -I validation_oxford_example.py <ExampleDC_C1.mat> <out_dir>

Source: C. Birkl, D. Howey (2017), "Oxford Battery Degradation Dataset 1", University of Oxford,
doi:10.5287/bodleian:KO2kdmYGg. Licence: ODC Open Database License (ODbL) v1.0; contents under DbCL v1.0.
The exported CSV is a derived database and is therefore also offered under ODbL v1.0 with attribution.
Cell: Kokam SLPB533459H4, 740 mAh pouch, tested at 40 degC chamber temperature, 1 s logging.
Sign convention in the export: I_A > 0 = discharge (the raw file uses negative mA for discharge).
"""
import sys, os, json
import numpy as np
from scipy.io import loadmat


def main():
    src, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    m = loadmat(src, squeeze_me=True, struct_as_record=False)["ExampleDC_C1"]
    res = {}
    for part in ("ch", "dc"):
        p = getattr(m, part)
        t = np.asarray(p.t, float); t = t - t[0]
        v = np.asarray(p.v, float); i = -np.asarray(p.i, float) / 1000.0; T = np.asarray(p.T, float)
        q = np.asarray(p.q, float)
        dis = np.clip(i, 0, None); reg = np.clip(i, None, 0)
        res[part] = dict(n=int(t.size), duration_s=float(t[-1]), dt_median_s=float(np.median(np.diff(t))),
                         V_start=float(v[0]), V_min=float(v.min()), V_max=float(v.max()), V_end=float(v[-1]),
                         I_mean_A=float(np.mean(i)), I_rms_A=float(np.sqrt(np.mean(i ** 2))),
                         I_max_dis_A=float(i.max()), I_max_regen_A=float(-i.min()),
                         Ah_out=float(np.trapezoid(dis, t) / 3600), Ah_in_regen=float(-np.trapezoid(reg, t) / 3600),
                         net_Ah=float(np.trapezoid(i, t) / 3600), q_logged_end_mAh=float(q[-1]),
                         E_net_Wh=float(np.trapezoid(i * v, t) / 3600),
                         T_min_C=float(T.min()), T_max_C=float(T.max()))
        if part == "dc":
            with open(os.path.join(out, "oxford_ExampleDC_C1_artemis_discharge.csv"), "w") as g:
                g.write("t_s,V,I_A,T_C\n")
                for a, b, c, d in zip(t, v, i, T):
                    g.write(f"{a:.0f},{b:.5f},{c:.5f},{d:.3f}\n")
    with open(os.path.join(out, "oxford_example_summary.json"), "w") as g:
        json.dump(dict(source="Birkl & Howey 2017, Oxford Battery Degradation Dataset 1, doi:10.5287/bodleian:KO2kdmYGg",
                       licence="ODbL v1.0 / DbCL v1.0", cell="Kokam SLPB533459H4 740 mAh pouch, 40 degC chamber",
                       sign="I_A > 0 = discharge", **res), g, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
