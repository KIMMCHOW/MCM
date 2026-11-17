"""Summarise iPhone 13 Pro / 14 Pro battery-level traces from the Bonn 'Battery Impact of Individual Hardware
Sensors in Modern iPhone Devices' dataset (Schweizer et al. 2026, bonndata doi:10.60507/FK2/LSJGWW, CC BY 4.0).

Usage:
    python3 -I validation_bonn_iphone_summary.py <dir_with_extracted_experiment_folders> <out.json>

Each device file holds iOS-reported battery level (1 % steps) vs Unix time; the app shows a black full-screen view
(OLED, auto-lock off, airplane mode on, Wi-Fi/Bluetooth off, auto-brightness off; see the paper's Table 1).
Reported per device: total runtime (startDate -> last logged sample), level at the last sample, the time spent
between the first sample at <= 95 % and the first at <= 15 % (an 80-point window), and the implied mean drain
rate in %/h over that window. Device model is not stored in the files; we infer the 13 Pro / 14 Pro split from the
idle runtimes (3 shortest = 13 Pro, matching the paper's 3 x 13 Pro + 6 x 14 Pro) and flag it as inferred.
"""
import sys, os, glob, json
import numpy as np


def first_time_at_or_below(t, lvl, x):
    idx = np.nonzero(lvl <= x + 1e-6)[0]
    return float(t[idx[0]]) if idx.size else float("nan")


def summarise_file(path):
    d = json.load(open(path))
    rows = sorted(((r["timestamp"], r["batteryLevel"]) for r in d["data"]), key=lambda r: r[0])
    t = np.array([r[0] for r in rows]); lvl = np.array([r[1] for r in rows])
    t0 = float(d.get("startDate", t[0]))
    t95, t15 = first_time_at_or_below(t, lvl, 0.95), first_time_at_or_below(t, lvl, 0.15)
    w = (t15 - t95) / 3600
    levels = {}
    for k in range(100, 0, -1):  # hours from startDate until the level first reads <= k %
        tk = first_time_at_or_below(t, lvl, k / 100)
        levels[k] = (tk - t0) / 3600 if tk == tk else None
    return dict(device=str(d.get("deviceName")), sensors=d.get("sensors"), refresh_s=d.get("refreshRate"),
                hours_to_level=levels,
                n=int(t.size), level_first=float(lvl[0]), level_last=float(lvl[-1]),
                runtime_h=float((t[-1] - t0) / 3600), t95_to_t15_h=float(w),
                drain_pct_per_h_95_15=float(80 / w) if w > 0 else float("nan"))


def main():
    root, out = sys.argv[1], sys.argv[2]
    res = []
    for exp in sorted(glob.glob(os.path.join(root, "*"))):
        if not os.path.isdir(exp):
            continue
        cfgs = glob.glob(os.path.join(exp, "*configuration.json"))
        cfg = json.load(open(cfgs[0])) if cfgs else {}
        for f in sorted(glob.glob(os.path.join(exp, "*.json"))):
            if f.endswith("configuration.json"):
                continue
            s = summarise_file(f); s.update(experiment=os.path.basename(exp), config=cfg); res.append(s)
    # infer model groups from the 9-s idle experiment(s)
    idle = [r for r in res if r["experiment"].endswith("Idle") and str(r["config"].get("samplingRate")) == "9"]
    mean_rt = {}
    for r in idle:
        mean_rt.setdefault(r["device"], []).append(r["runtime_h"])
    order = sorted(mean_rt, key=lambda k: np.mean(mean_rt[k]))
    group = {k: ("13 Pro (inferred)" if i < 3 else "14 Pro (inferred)") for i, k in enumerate(order)}
    for r in res:
        r["model_inferred"] = group.get(r["device"], "unknown")
    # group statistics
    stats = {}
    for r in res:
        key = f"{r['experiment']}|{r['model_inferred']}"
        stats.setdefault(key, []).append(r)
    summary = {}
    for k, rs in sorted(stats.items()):
        rt = np.array([x["runtime_h"] for x in rs]); dr = np.array([x["drain_pct_per_h_95_15"] for x in rs])
        summary[k] = dict(n_devices=len(rs), config=rs[0]["config"], runtime_h_median=float(np.median(rt)),
                          runtime_h_min=float(rt.min()), runtime_h_max=float(rt.max()),
                          drain_pct_per_h_median=float(np.nanmedian(dr)),
                          level_last_max=float(max(x["level_last"] for x in rs)))
    # median SOC(t) curve per experiment and model group -> CSV next to the JSON
    keys = sorted(stats)
    with open(os.path.splitext(out)[0] + "_soc_curves.csv", "w") as g:
        g.write("level_pct," + ",".join(k.replace(" (inferred)", "").replace("|", "_").replace(" ", "") + "_h" for k in keys) + "\n")
        for lev in range(100, 0, -1):
            vals = []
            for k in keys:
                v = [x["hours_to_level"][lev] for x in stats[k] if x["hours_to_level"][lev] is not None]
                vals.append(f"{np.median(v):.3f}" if v else "")
            g.write(f"{lev}," + ",".join(vals) + "\n")
    for r in res:
        r.pop("hours_to_level", None)
    with open(out, "w") as g:
        json.dump(dict(source="Schweizer et al. 2026, bonndata doi:10.60507/FK2/LSJGWW (CC BY 4.0); "
                              "paper doi:10.3390/s26102923", model_mapping_inferred=group,
                       groups=summary, per_device=res), g, indent=1)
    print("inferred model map:", group)
    for k, s in summary.items():
        print(f"{k:42s} n={s['n_devices']} cfg={s['config']} runtime median {s['runtime_h_median']:.2f} h "
              f"[{s['runtime_h_min']:.2f}, {s['runtime_h_max']:.2f}] drain(95->15%) {s['drain_pct_per_h_median']:.3f} %/h "
              f"last-level max {s['level_last_max']:.2f}")


if __name__ == "__main__":
    main()
