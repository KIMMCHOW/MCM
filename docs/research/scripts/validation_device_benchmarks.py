"""Derived average battery-side power for published smartphone battery-life figures (validation targets).

Usage:  python3 -I validation_device_benchmarks.py <bonn_iphone_summary.json> <out.json>

P_avg = E_batt / TTE, with E_batt the rated pack energy. For Android phones only mAh is published, so
E = Q * V_nom with V_nom = 3.87 V (range 3.85-3.89 V), the nominal voltage implied by Apple's own Wh/mAh pairs
(12.38 Wh / 3200 mAh = 3.869 V; 11.97 Wh / 3095 mAh = 3.867 V; 16.68 Wh / 4323 mAh = 3.858 V; 9to5Mac 2022).
Every runtime below is copied from the source named in `src` (keys resolve to validation_data.md, section e).
"""
import sys, json

VNOM, VLO, VHI = 3.87, 3.85, 3.89
DEV = {  # name: (E_Wh or None, Q_mAh or None, energy source key)
    "iPhone 13 Pro": (11.97, 3095, "N9M22"),
    "iPhone 14 Pro": (12.38, 3200, "N9M22"),
    "iPhone 14 Pro Max": (16.68, 4323, "NBC22"),
    "Pixel 8": (None, 4575, "GPX"),
    "Pixel 8a": (None, 4492, "GPX"),
    "Galaxy S24 (Exynos 2400)": (None, 4000, "NBC24"),
}
ROWS = [  # device, test, conditions, hours, source
    ("iPhone 14 Pro", "video playback (local)", "Apple lab: looped 2h23m movie, Wi-Fi associated, BT headphones, Auto-Brightness/True Tone off", 23, "APL14P"),
    ("iPhone 14 Pro", "video playback (streamed)", "Apple lab: looped 3h01m HDR movie streamed", 20, "APL14P"),
    ("iPhone 14 Pro", "audio playback", "Apple lab: 358-track AAC playlist, BT headphones, display not visible (AOD on)", 75, "APL14P"),
    ("iPhone 13 Pro", "video playback (local)", "Apple lab (Aug 2021), same protocol family", 22, "APL13P"),
    ("iPhone 13 Pro", "video playback (streamed)", "Apple lab", 20, "APL13P"),
    ("iPhone 13 Pro", "audio playback", "Apple lab", 75, "APL13P"),
    ("iPhone 14 Pro", "web browsing (cellular)", "Laptop Mag: continuous web surfing, 150 nits, T-Mobile eSIM (mobile network), mean of 3 runs", 10 + 22 / 60, "LMAG22"),
    ("iPhone 13 Pro", "web browsing (cellular)", "Laptop Mag: same test", 12.0, "LMAG22"),
    ("iPhone 14 Pro Max", "idle screen-on, min brightness, no WLAN", "Notebookcheck Reader/Idle", 57 + 53 / 60, "NBC22"),
    ("iPhone 14 Pro Max", "web browsing (Wi-Fi)", "Notebookcheck WiFi websurfing, fixed 150 cd/m2", 15 + 48 / 60, "NBC22"),
    ("iPhone 14 Pro Max", "video playback (local H.264 1080p)", "Notebookcheck Big Buck Bunny loop", 28 + 8 / 60, "NBC22"),
    ("iPhone 14 Pro Max", "max load (3D benchmark), max brightness", "Notebookcheck Load", 6 + 56 / 60, "NBC22"),
    ("Pixel 8", "idle screen-on, min brightness, no WLAN", "Notebookcheck Reader/Idle", 24 + 53 / 60, "NBC23"),
    ("Pixel 8", "web browsing (Wi-Fi)", "Notebookcheck WiFi v1.3 (150 cd/m2 script)", 13 + 34 / 60, "NBC23"),
    ("Pixel 8", "video playback (local H.264 1080p)", "Notebookcheck Big Buck Bunny loop", 22 + 43 / 60, "NBC23"),
    ("Pixel 8", "max load (Burnout benchmark), max brightness", "Notebookcheck Load", 3 + 25 / 60, "NBC23"),
    ("Galaxy S24 (Exynos 2400)", "idle screen-on, min brightness, no WLAN", "Notebookcheck Reader/Idle", 39 + 58 / 60, "NBC24"),
    ("Galaxy S24 (Exynos 2400)", "web browsing (Wi-Fi)", "Notebookcheck WiFi v1.3, Chrome 122, 150 cd/m2", 14 + 40 / 60, "NBC24"),
    ("Galaxy S24 (Exynos 2400)", "video playback (local H.264 1080p)", "Notebookcheck Big Buck Bunny loop", 20 + 27 / 60, "NBC24"),
    ("Galaxy S24 (Exynos 2400)", "max load (Burnout benchmark), max brightness", "Notebookcheck Load", 2 + 39 / 60, "NBC24"),
    ("Pixel 8a", "typical mixed use (median-user profile)", "Google: talk/data/standby mix, US carrier, late 2023; 'average ... at least 31 hours'", 31, "GPX"),
]


def main():
    bonn = json.load(open(sys.argv[1]))["groups"]
    rows = list(ROWS)
    pick = {"07-09-2023-Idle": "idle, black full-screen OLED, airplane mode (9 s log)",
            "02-10-2023-GPS": "GPS 'best' accuracy, black screen, airplane mode",
            "09-11-2023-CAM": "rear camera 30 fps (frames not rendered), black screen",
            "23-11-2023-CXL": "camera 60 fps + GPS best, black screen"}
    for exp, label in pick.items():
        for model, dev in (("13 Pro (inferred)", "iPhone 13 Pro"), ("14 Pro (inferred)", "iPhone 14 Pro")):
            g = bonn[f"{exp}|{model}"]
            rows.append((dev, label, f"Bonn dataset, median of {g['n_devices']} devices, 100% -> 1%", g["runtime_h_median"], "BONN26"))
    out = []
    for dev, test, cond, h, src in rows:
        E, Q, esrc = DEV[dev]
        if E is None:
            e, elo, ehi, ekind = Q / 1000 * VNOM, Q / 1000 * VLO, Q / 1000 * VHI, f"{Q} mAh x 3.87 V (assumed)"
        else:
            e = elo = ehi = E; ekind = f"{E} Wh rated"
        out.append(dict(device=dev, test=test, conditions=cond, hours=round(h, 3), source=src, E_Wh=round(e, 2),
                        E_basis=ekind, E_source=esrc, P_avg_W=round(e / h, 3),
                        P_range_W=[round(elo / h, 3), round(ehi / h, 3)]))
    # energy-free additivity check on the Bonn runs: 1/TTE(cam60+GPS) =? 1/TTE(cam60) + 1/TTE(GPS) - 1/TTE(idle)
    add = {}
    for model in ("13 Pro (inferred)", "14 Pro (inferred)"):
        T = lambda e: bonn[f"{e}|{model}"]["runtime_h_median"]
        pred = 1 / (1 / T("14-11-2023-CAM") + 1 / T("02-10-2023-GPS") - 1 / T("07-09-2023-Idle"))
        add[model] = dict(TTE_pred_h=round(pred, 3), TTE_obs_h=round(T("23-11-2023-CXL"), 3),
                          rel_err_pct=round(100 * (pred / T("23-11-2023-CXL") - 1), 2))
        print("additivity", model, add[model])
    json.dump(dict(note=__doc__, rows=out, additivity_check_bonn=add), open(sys.argv[2], "w"), indent=1)
    for r in out:
        print(f"{r['device']:26s} {r['test'][:46]:46s} {r['hours']:7.2f} h  E={r['E_Wh']:5.2f} Wh  P={r['P_avg_W']:.3f} W  [{r['source']}]")


if __name__ == "__main__":
    main()
