"""AN1: primitive powers and component breakdown (FINAL_SPEC §3, §5.1, §8 V6, §14).

Outputs
  results/power.json                 primitives.{ID}.{P_W, P_capped_W, components{...}, dPb_W, ...},
                                     v6.{case}.{model_P_W, model_tte_h, obs_h, err_pct, pass, P_implied_W, ...},
                                     bands, proto_check
  paper/figures/fig_power_breakdown.pdf
      (a) RP-63 primitives at 25 °C: stacked component bars (demand before throttling), F5 capped game power,
          published max-load band;  (b) V6 device configurations: model components vs the power implied by the
          published runtime (F2 inverted on the device's own battery), filled = within the ±25 % band.
  paper/tables/tab_scenarios.tex     condensed primitives table (T2)

``P_W`` is the battery-side demand before throttling (the 'P [PROTO]' column of §5.1); ``P_capped_W`` is the F5
quasi-steady power at T_amb = 25 °C (differs only for the game states).  Component powers are the E1 terms; the
radio entry includes the idle power of every associated radio.
"""
from __future__ import annotations

import time

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.optimize import brentq

from basics_common import f, load_proto, meta, write_json, write_tex
from phonebatt import PRIMITIVES, RP63, analytic, cell, power
from phonebatt import validation as V
from phonebatt.usage import Schedule
from style import (AQUA, BLUE, GREEN, GRID, INK, INK2, MAGENTA, ORANGE, VIOLET, W1, YELLOW, apply,
                   panel_label, save)

ORDER = ["SB-W", "SB-L", "SB-5", "RD", "SO", "SO-L", "VW", "VL", "V5", "NV", "GM", "GM-L", "VC", "VC-W", "CM", "AU",
         "AU-L", "TRK", "BU"]
PROTO_P = {"SB-W": 0.040, "SB-L": 0.047, "SB-5": 0.102, "RD": 0.456, "SO": 1.060, "SO-L": 1.218, "VW": 0.548,
           "VL": 0.694, "V5": 0.837, "NV": 1.486, "GM": 3.906, "GM-L": 4.302, "VC": 2.518, "VC-W": 1.622,
           "CM": 1.687, "AU": 0.224, "AU-L": 0.261, "TRK": 0.464, "BU": 4.095}  # §5.1 'P [PROTO]' column
LABEL = {"SB-W": "standby, Wi-Fi", "SB-L": "standby, LTE", "SB-5": "standby, 5G push", "RD": "reading",
         "SO": "social feed", "SO-L": "social, LTE", "VW": "video stream", "VL": "video, LTE", "V5": "video, 5G",
         "NV": "navigation", "GM": "3D game", "GM-L": "3D game, LTE", "VC": "video call, LTE",
         "VC-W": "video call", "CM": "camera", "AU": "audio, screen off", "AU-L": "audio, LTE",
         "TRK": "GNSS tracking", "BU": "photo burst (S9)"}
USED_IN = {"SB-W": "S1, S12, C", "SB-L": "S2a, C", "SB-5": "S2b", "RD": "S3, S12", "SO": "S4, S10a, S12, C",
           "SO-L": "S11a, C", "VW": "S5a, S12, C", "VL": "S5b, C", "V5": "S5c", "NV": "S6, S12, C",
           "GM": "S7, S12, C", "GM-L": "C", "VC": "S8a, S11c, S12, C", "VC-W": "S8b, C", "CM": "CM, C",
           "AU": "S16, C", "AU-L": "C", "TRK": "S9", "BU": "S9"}

# component groups of the figure: fixed categorical order (entity -> slot), also the stacking order
COMP = [("SoC (CPU+GPU)", ("soc",), BLUE), ("display", ("disp",), ORANGE), ("radio", ("net",), AQUA),
        ("floor (suspend, awake, wake-ups)", ("susp", "awake", "wake"), YELLOW), ("camera", ("cam",), MAGENTA),
        ("GNSS + sensors", ("gnss", "sens"), GREEN), ("audio / BT", ("aud",), VIOLET)]

# V6 case id -> (short label, device_benchmarks.json (device, test) for the rated-energy implied power)
V6_LABEL = {
    "14P_local_video": ("14 Pro: local video", ("iPhone 14 Pro", "video playback (local)")),
    "14P_stream_video": ("14 Pro: streamed video", ("iPhone 14 Pro", "video playback (streamed)")),
    "14P_audio_BT": ("14 Pro: audio, BT", ("iPhone 14 Pro", "audio playback")),
    "14PM_wifi_web_150": ("14 Pro Max: web, 150 nit", ("iPhone 14 Pro Max", "web browsing (Wi-Fi)")),
    "14PM_local_video_150": ("14 Pro Max: local video", ("iPhone 14 Pro Max", "video playback (local H.264 1080p)")),
    "14PM_idle_minbright": ("14 Pro Max: idle, min. bright.",
                            ("iPhone 14 Pro Max", "idle screen-on, min brightness, no WLAN")),
    "14P_cell_web_150": ("14 Pro: LTE web, 150 nit", ("iPhone 14 Pro", "web browsing (cellular)")),
    "Pixel8_wifi_web": ("Pixel 8: web", ("Pixel 8", "web browsing (Wi-Fi)")),
    "Pixel8_local_video": ("Pixel 8: local video", ("Pixel 8", "video playback (local H.264 1080p)")),
    "Pixel8_idle_minbright": ("Pixel 8: idle, min. bright.", ("Pixel 8", "idle screen-on, min brightness, no WLAN")),
    "S24_wifi_web": ("Galaxy S24: web", ("Galaxy S24 (Exynos 2400)", "web browsing (Wi-Fi)")),
    "S24_local_video": ("Galaxy S24: local video", ("Galaxy S24 (Exynos 2400)", "video playback (local H.264 1080p)")),
    "S24_idle_minbright": ("Galaxy S24: idle, min. bright.",
                           ("Galaxy S24 (Exynos 2400)", "idle screen-on, min brightness, no WLAN")),
    "Pixel8a_S12_day": ("Pixel 8a: S12 day (bound)", ("Pixel 8a", "typical mixed use (median-user profile)")),
}
V6_PROTO_H = {"14P_local_video": 26.08, "14P_stream_video": 23.40, "14P_audio_BT": 55.90, "14PM_wifi_web_150": 15.18,
              "14PM_local_video_150": 33.01, "14PM_idle_minbright": 66.63, "14P_cell_web_150": 10.29,
              "Pixel8_wifi_web": 16.37, "Pixel8_local_video": 35.70, "Pixel8_idle_minbright": 70.52,
              "S24_wifi_web": 14.30, "S24_local_video": 31.20, "S24_idle_minbright": 61.64, "Pixel8a_S12_day": 34.67}


def components(u, p, Tamb_C=25.0) -> dict:
    """E1 terms of one constant input [W]: demand (soc_dem) and the F5-capped SoC power (soc_capped)."""
    Pns, dem, dPb, br = power.power_parts(u, p)
    out = {k: f(v) for k, v in br.items() if k != "soc_dem"}
    out["soc_dem"] = f(dem)
    out["soc_capped"] = f(power.steady_power(u, p, Tamb_C)) - f(Pns)
    return out, f(dPb)


def schedule_components(sch: Schedule, p, Tamb_C=25.0) -> dict:
    """Time-averaged E1 terms of a cyclic schedule with quasi-steady (F5) SoC power per segment [W]."""
    acc: dict[str, float] = {}
    for d, u in sch:
        c, _ = components(u, p, Tamb_C)
        for k, v in c.items():
            acc[k] = acc.get(k, 0.0) + d * v
    return {k: v / sch.period for k, v in acc.items()}


def group(comp: dict, capped: bool) -> list[float]:
    vals = []
    for _, keys, _ in COMP:
        v = 0.0
        for k in keys:
            v += comp["soc_capped" if (k == "soc" and capped) else ("soc_dem" if k == "soc" else k)]
        vals.append(v)
    return vals


def implied_power(obs_h: float, p) -> float:
    """Constant battery power whose exact quasi-static TTE (F2, 25 °C, from 100 %) equals obs_h [W]."""
    return brentq(lambda P: f(analytic.tte_L1(P, p, 25.0)) - obs_h, 1e-3, 20.0, xtol=1e-9)


def main() -> None:
    t0 = time.time()
    proto = load_proto()
    # ---------------------------------------------------------------------------------------- primitives
    prim = {}
    dev_P, dev_comp = [], []
    for n in ORDER:
        u = PRIMITIVES[n]
        comp, dPb = components(u, RP63)
        P = f(power.demand_power(u, RP63))
        Pc = f(power.steady_power(u, RP63, 25.0))
        Pp, brp = proto.Pmean(n, proto.P0)
        dcomp = max(abs(comp.get("soc_dem" if k == "soc_dem" else k, 0.0) - v) for k, v in brp.items())
        prim[n] = dict(label=LABEL[n], P_W=P, P_capped_W=Pc, dPb_W=dPb, P_peak_W=Pc + dPb, components=comp,
                       groups=dict(zip([c[0] for c in COMP], group(comp, capped=False))),
                       proto_P_W=PROTO_P[n], dev_vs_spec_mW=1e3 * (P - PROTO_P[n]),
                       dev_vs_proto_module_W=P - Pp, max_component_dev_vs_proto_W=dcomp, used_in=USED_IN[n])
        dev_P.append(abs(P - PROTO_P[n]))
        dev_comp.append(dcomp)
    gm = analytic.thermal_cap("GM", RP63, 25.0)

    # ---------------------------------------------------------------------------------------- V6 band
    v6p = V.v6_predictions()
    bench = {(r["device"], r["test"]): r for r in V.load_device_benchmarks()["rows"]}
    v6 = {}
    for c in V.v6_cases():
        cid = c["id"]
        pr = v6p[cid]
        if isinstance(c["u"], Schedule):
            comp = schedule_components(c["u"], c["p"])
            P_dem = sum(group(comp, capped=True))
        else:
            comp, _ = components(c["u"], c["p"])
            P_dem = f(power.demand_power(c["u"], c["p"]))
        Pimp = implied_power(c["obs_h"], c["p"])
        b = bench.get(V6_LABEL[cid][1])
        v6[cid] = dict(label=V6_LABEL[cid][0], src=c["src"], model_P_W=P_dem, model_P_mean_L2_W=pr["P_mean_W"],
                       model_tte_h=pr["model_h"], obs_h=c["obs_h"], err_pct=pr["err_pct"], pass_=pr["pass_"],
                       bound=c.get("bound"), P_implied_W=Pimp,
                       P_implied_rated_W=None if b is None else b["P_avg_W"],
                       E_u_device_Wh=f(cell.usable_energy(c["p"])), components=comp,
                       groups=dict(zip([g[0] for g in COMP], group(comp, capped=True))),
                       proto_tte_h=V6_PROTO_H[cid], dev_vs_proto_pct=100 * (pr["model_h"] / V6_PROTO_H[cid] - 1))
    Eu = f(cell.usable_energy(RP63))
    bands = dict(
        max_load_W=dict(lo=min(bench[(d, t)]["P_avg_W"] for d, t in bench if t.startswith("max load")),
                        hi=max(bench[(d, t)]["P_avg_W"] for d, t in bench if t.startswith("max load")),
                        src="NBC22, NBC23, NBC24 (rated-energy basis, device_benchmarks.json)",
                        model_GM_capped_W=f(gm["P_star"])),
        standby_drain=dict(lo_pct_per_h=0.15, hi_pct_per_h=0.5, src="device_power BG9",
                           lo_W=0.0015 * Eu, hi_W=0.005 * Eu, model_SBW_W=prim["SB-W"]["P_W"],
                           model_S1_pct_per_h=None))
    n_pass = sum(v["pass_"] for v in v6.values())

    out = dict(meta=meta("run_power.py", solver="LSODA (V6 TTE); algebraic (powers)", h_s=None,
                         note="powers at T_amb = 25 °C, RP-63 defaults; P_capped_W = F5 fixed point"),
               primitives=prim,
               game_cap_25C=dict(P_star_W=f(gm["P_star"]), rho_perf=f(gm["rho_perf"]), P_nonsoc_W=f(gm["P_nonsoc"]),
                                 sharp_cap_W=f(gm["sharp_cap_W"])),
               v6=v6, v6_n_pass=n_pass, v6_n_cases=len(v6), bands=bands,
               proto_check=dict(max_abs_dev_P_vs_spec_W=max(dev_P), tol_W=5e-4,
                                max_abs_dev_components_vs_proto_module_W=max(dev_comp),
                                max_abs_dev_v6_tte_vs_spec_pct=max(abs(v["dev_vs_proto_pct"]) for v in v6.values())))
    # standby drain of S1 in %/h from the L2 run is in AN3; here the L0 equivalent
    out["bands"]["standby_drain"]["model_S1_pct_per_h"] = 100 * prim["SB-W"]["P_W"] / Eu
    write_json("power.json", out)

    figure(prim, v6, bands, gm)
    table(prim)
    print(f"primitives: max |P - PROTO| = {1e3 * max(dev_P):.3f} mW (tol 0.5 mW); max component dev vs proto "
          f"module = {max(dev_comp):.2e} W")
    print(f"GM capped 25 C: {f(gm['P_star']):.4f} W (PROTO 3.116); V6 pass {n_pass}/{len(v6)}; max V6 TTE dev vs "
          f"PROTO {out['proto_check']['max_abs_dev_v6_tte_vs_spec_pct']:.3f} %")
    for k, v in v6.items():
        print(f"  {k:24s} P_model={v['model_P_W']:.3f} P_impl={v['P_implied_W']:.3f} (rated {v['P_implied_rated_W']})"
              f" TTE {v['model_tte_h']:.2f}/{v['obs_h']:.2f} h err {v['err_pct']:+.1f} % pass={v['pass_']}")
    print(f"done in {time.time() - t0:.1f} s")


def figure(prim, v6, bands, gm) -> None:
    apply()
    fig = plt.figure(figsize=(W1, 3.15))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.78, left=0.14, right=0.93, top=0.83, bottom=0.12)
    axa, axb = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    h = 0.72
    # (a) primitives -----------------------------------------------------------------------------------
    names = ORDER[::-1]
    y = np.arange(len(names))
    for yi, n in zip(y, names):
        vals = group(prim[n]["components"], capped=False)
        left = 0.0
        for (lab, _, col), v in zip(COMP, vals):
            if v > 0:
                axa.barh(yi, v, left=left, height=h, color=col, edgecolor="white", linewidth=0.3)
            left += v
        axa.text(left + 0.06, yi, f"{prim[n]['P_W']:.2f}", va="center", ha="left", fontsize=6.3, color=INK2)
    lo, hi = bands["max_load_W"]["lo"], bands["max_load_W"]["hi"]
    yg = names.index("GM")
    axa.barh(yg, hi - lo, left=lo, height=0.98, color=GRID, edgecolor="none", zorder=0)
    for n in ("GM", "GM-L"):
        yi = names.index(n)
        Pc = prim[n]["P_capped_W"]
        axa.plot([Pc, Pc], [yi - 0.48, yi + 0.48], color=INK, lw=1.3, solid_capstyle="butt", zorder=4)
    leg = [Line2D([], [], color=INK, lw=1.3, label="F5 thermal cap, 25 °C"),
           Patch(facecolor=GRID, edgecolor="none", label=f"published max load {lo:.1f}–{hi:.1f} W")]
    axa.legend(handles=leg, loc="lower right", bbox_to_anchor=(1.0, 0.075), fontsize=6.3, handlelength=1.0,
               handletextpad=0.4, borderaxespad=0.2)
    axa.set_yticks(y)
    axa.set_yticklabels([f"{n}  {prim[n]['label']}" for n in names], fontsize=6.6)
    axa.set_xlim(0, 6.0)
    axa.set_ylim(-0.6, len(names) - 0.4)
    axa.set_xlabel("battery-side power (W), 25 °C, before throttling")
    axa.grid(axis="y", visible=False)
    axa.tick_params(axis="y", length=0)
    panel_label(axa, "(a)")
    axa.set_title("RP-63 activity primitives", loc="left", fontsize=8.5, pad=4)

    # (b) V6 device configurations ----------------------------------------------------------------------
    ids = list(v6)[::-1]
    yv = np.arange(len(ids))
    for yi, cid in zip(yv, ids):
        vals = group(v6[cid]["components"], capped=True)
        left = 0.0
        for (lab, _, col), v in zip(COMP, vals):
            if v > 0:
                axb.barh(yi, v, left=left, height=h, color=col, edgecolor="white", linewidth=0.3)
            left += v
        Pi = v6[cid]["P_implied_W"]
        ok = v6[cid]["pass_"]
        mk = "<" if v6[cid]["bound"] == "lower" else "D"
        axb.plot([Pi], [yi], marker=mk, ms=4.2, mfc=INK if ok else "white", mec=INK, mew=0.9, ls="none", zorder=5)
        e = v6[cid]["err_pct"]
        axb.text(1.02, yi, f"{e:+.0f}%", transform=axb.get_yaxis_transform(), va="center", ha="left", fontsize=6.3,
                 color=INK2)
    axb.text(1.02, len(ids) - 0.3, "ε$_\\mathrm{TTE}$", transform=axb.get_yaxis_transform(), va="bottom",
             ha="left", fontsize=6.5, color=INK2)
    axb.set_yticks(yv)
    axb.set_yticklabels([v6[c]["label"] for c in ids], fontsize=6.6)
    xmax = max(max(v["P_implied_W"], sum(v["groups"].values())) for v in v6.values())
    axb.set_xlim(0, 1.08 * xmax)
    axb.set_ylim(-0.6, len(ids) - 0.4)
    axb.set_xlabel("battery-side power (W)")
    axb.grid(axis="y", visible=False)
    axb.tick_params(axis="y", length=0)
    panel_label(axb, "(b)")
    axb.set_title("Published runtimes (V6)", loc="left", fontsize=8.5, pad=4)

    handles = [Patch(facecolor=c, edgecolor="none", label=l) for l, _, c in COMP]
    handles += [Line2D([], [], marker="D", ms=4, mfc=INK, mec=INK, ls="none", label="implied by runtime, |ε| ≤ 25 %"),
                Line2D([], [], marker="D", ms=4, mfc="white", mec=INK, ls="none", label="implied, outside band")]
    fig.legend(handles=handles, loc="upper center", ncol=5, bbox_to_anchor=(0.53, 1.0), fontsize=6.6,
               handlelength=1.2, columnspacing=1.0, handletextpad=0.4)
    save(fig, "fig_power_breakdown")


def table(prim) -> None:
    def num(x, nd=2):
        return f"{x:.{nd}f}"

    rows = []
    for n in ORDER:
        u = PRIMITIVES[n]
        if u.s_scr:
            disp = f"{num(u.b, 1)}/{num(u.APL, 2)}/{int(u.f_r)}"
        else:
            disp = "off"
        if u.rho == "none":
            data = "--"
        else:
            path = "Wi-Fi" if u.rho == "wifi" else {"lte": "LTE", "nr": "5G"}[u.c]
            if u.kappa == "rt":
                data = f"{path}, real-time"
            else:
                lam = u.lam
                lam_s = "0" if lam == 0 else (f"1/{round(1 / lam)}" if lam < 1 else f"{lam:g}")
                dd = f", {u.d_dl:g}" if u.d_dl else ""
                du = f", UL {u.d_ul:g}" if u.d_ul else ""
                data = f"{path}, {lam_s}{dd}{du}"
        other = []
        if u.g:
            other.append("GNSS")
        if u.c_r:
            other.append("cam$_r$")
        if u.c_f:
            other.append("cam$_f$")
        if u.au:
            other.append("aud")
        if u.bt:
            other.append("BT")
        if (u.awake and not u.s_scr):
            other.append("WL")
        P = prim[n]["P_W"]
        capped = abs(prim[n]["P_capped_W"] - P) > 5e-3 and n != "BU"  # BU lasts 30 s: no steady cap
        Ps = f"{num(P, 3)} ({num(prim[n]['P_capped_W'], 2)})" if capped else num(P, 3)
        xc = "0" if u.x_c == 0 else num(u.x_c, 2)
        rows.append(f"{n} & {prim[n]['label']} & {disp} & {xc} & {data} & {', '.join(other) or '--'} & "
                    f"{u.dPb:g} & {Ps} & {prim[n]['used_in']} \\\\")
    lines = [
        "% T2: activity primitives (generated by analysis/run_power.py; do not edit by hand).",
        "% b/APL/f_r: slider, average picture level, refresh rate (Hz); off = screen off. Data: path, burst rate (1/s),",
        "% mean downlink (Mbit/s). P: battery-side demand at 25 degC before throttling; in parentheses the F5",
        "% thermally capped power. Delta P_b (W): sub-second burst envelope amplitude (x m_b). Other: cam_r/cam_f rear/front",
        "% camera, aud audio, BT Bluetooth audio, WL wakelock (screen off, processor awake). C = usage CTMC state.",
        "\\footnotesize",
        "\\setlength{\\tabcolsep}{3pt}",
        "\\begin{tabular}{@{}llllllrrl@{}}",
        "\\toprule",
        "ID & Activity & $b$/APL/$f_r$ & $x_c$ & Data, $\\lambda$ (s$^{-1}$), $d_{\\mathrm{dl}}$ & Other & "
        "$\\Delta P_b$ & $P$ (W) & Used in \\\\",
        "\\midrule",
        *rows,
        "\\bottomrule",
        "\\end{tabular}",
        "",
    ]
    write_tex("tab_scenarios.tex", "\n".join(lines))


if __name__ == "__main__":
    main()
