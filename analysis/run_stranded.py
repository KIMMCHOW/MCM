"""AN12: cold / aged stranded charge and the O1 peak-power governor (FINAL_SPEC §5.4, §7 F4, §12.2 O1).

Writes
* results/stranded.json  - map."{Tamb}|{SOH}".{s_off_pct, tte_h, event, ...} (S10a social feed with the 3 W
  burst envelope), map_governed (same with the online O1 governor), governor."{m_b}".{s_off_pct, tte_h, event}
  (aged pack at -10 °C), governor_map."{m_b}"."{Tamb}|{SOH}", S9 grid, proto_comparison, meta
* paper/figures/fig_stranded_map.pdf, paper/tables/tab_stranded.tex

Run:  cd /home/user/MCM && PYTHONPATH=src:analysis python3 analysis/run_stranded.py

O1 governor (§12.2): the sub-second burst amplitude is capped online at ΔP_b,eff = min(ΔP_b, max(0, 0.9·ΔP_b,safe))
with ΔP_b,safe = V_sys(Ũ − V_sys)/R0 − P (F4, evaluated on the model state).  In the L2 model the burst envelope
only enters the brown-out event g3 (it carries no mean energy, §2.9), and with the cap P + ΔP_b,eff < P_cut,sys
whenever ΔP_b,safe > 0, i.e. V_pk > V_sys; when ΔP_b,safe ≤ 0 the mean-load voltage is already ≤ V_sys < V_cut,
so the cut-off has fired first.  The governed trajectory is therefore *exactly* the trajectory with the envelope
switched off (p.burst = False); the script checks V_pk,gov > V_sys along every governed trajectory and reports
the performance cost (fraction of the requested burst amplitude delivered) from the recorded state.
"""
from __future__ import annotations

import json
import math
import subprocess
import time
from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize

from phonebatt import RP63, scenario
from phonebatt.cell import R0 as R0f
from phonebatt.cell import peak_voltage
from style import INK, INK2, MUTED, SEQ_BLUE, W1, apply, cmap_seq, panel_label, save

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
TAB = ROOT / "paper" / "tables"
SPEC = "FINAL_SPEC 2026-11-17"

TGRID = (25, 10, 0, -10, -20)
SOHS = (1.0, 0.9, 0.8)
MB = (1.0, 0.75, 0.5, 0.25, 0.0)
GOV_FRAC = 0.9  # O1 cap: 0.9·ΔP_b,safe

# §5.4 [PROTO] stranded-charge map: (s_off %, TTE h); brownout except 25 °C / SOH 1.0 (cutoff)
PROTO_MAP = {
    (25, 1.0): (0.4, 16.64), (25, 0.9): (0.5, 14.94), (25, 0.8): (0.9, 13.22),
    (10, 1.0): (0.7, 16.57), (10, 0.9): (1.3, 14.80), (10, 0.8): (2.1, 13.05),
    (0, 1.0): (1.4, 16.43), (0, 0.9): (2.4, 14.62), (0, 0.8): (3.8, 12.80),
    (-10, 1.0): (2.8, 16.18), (-10, 0.9): (4.9, 14.22), (-10, 0.8): (8.5, 12.15),
    (-20, 1.0): (6.4, 15.54), (-20, 0.9): (13.1, 12.99), (-20, 0.8): (23.3, 10.23)}
# §12.2 O1 [PROTO]: aged pack (SOH 0.8) at -10 °C, social: m_b -> (s_off %, TTE h)
PROTO_GOV = {1.0: (8.5, 12.15), 0.5: (2.6, 12.88), 0.0: (1.5, 13.00)}
# §5.4 [PROTO] S9: (TTE h, s_off %), all cutoff during a BU burst
PROTO_S9 = {(25, 1.0): (26.99, 1.7), (-10, 1.0): (25.33, 7.2), (-20, 1.0): (22.67, 16.9),
            (25, 0.8): (21.32, 2.7), (-10, 0.8): (17.00, 22.2), (-20, 0.8): (12.33, 44.2)}


def clean(x):
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, np.ndarray):
        return clean(x.tolist())
    if isinstance(x, (np.floating, float)):
        x = float(x)
        return x if math.isfinite(x) else None
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def git_meta():
    def sh(*a):
        try:
            return subprocess.run(a, cwd=ROOT, capture_output=True, text=True, timeout=10).stdout.strip()
        except Exception:  # noqa: BLE001
            return ""
    return dict(git_sha=sh("git", "rev-parse", "--short", "HEAD"),
                core_uncommitted=bool(sh("git", "status", "--porcelain", "src/phonebatt")))


def mk(T, soh):
    return f"{T:g}|{soh:.1f}"


def rec(r) -> dict:
    return dict(s_off_pct=100 * r.s_off, tte_h=r.tte_h, event=r.event, soe_off_pct=100 * r.soe_off,
                T_end_C=r.T_end_C, Ts_end_C=r.Ts_end_C, V_end=r.V_end)


def run_s4(T, soh, p=RP63, record=False):
    sc = scenario("S4", Tamb_C=float(T), soh=float(soh), p=p)
    return sc, sc.run(record=record)


def governor_cost(sc, r) -> dict:
    """Performance cost of the online O1 cap along a governed (envelope-off) trajectory, plus the exactness check
    V_pk,gov > V_sys.  φ(t) = ΔP_b,eff/ΔP_b is the delivered fraction of the requested burst amplitude."""
    tr = r.traj
    p = sc.p
    u = sc.schedule[0][1]
    dPb = float(u.dPb) * float(RP63.m_b)  # requested amplitude (m_b = 1)
    z, T, I, V, P, t = tr["z"], tr["T"], tr["I"], tr["V"], tr["P"], tr["t"]
    R0 = R0f(z, T, p)
    Ue = V + R0 * I  # Ũ = U − v1 (from V = Ũ − R0 I)
    Vs = float(p.V_sys)
    safe = Vs * (Ue - Vs) / R0 - P
    cap = np.minimum(dPb, np.maximum(0.0, GOV_FRAC * safe))
    phi = cap / dPb
    Vpk = peak_voltage(Ue, R0, P + cap)
    margin = float(np.min(Vpk - Vs))
    dur = t[-1] - t[0]
    w = np.diff(t)
    mean_phi = float(np.sum(0.5 * (phi[1:] + phi[:-1]) * w) / dur)
    capped = phi < 1.0 - 1e-12
    frac_capped = float(np.sum(w[capped[:-1]]) / dur)
    onset = int(np.argmax(capped)) if capped.any() else None
    s_onset = float(100 * tr["s"][onset]) if onset is not None else None
    t_onset = float(t[onset] / 3600.0) if onset is not None else None
    # mean delivered fraction over the capped part only
    if capped.any():
        tc = t[onset:]
        pc = phi[onset:]
        mean_phi_capped = float(np.sum(0.5 * (pc[1:] + pc[:-1]) * np.diff(tc)) / (tc[-1] - tc[0]))
    else:
        mean_phi_capped = None
    return dict(burst_requested_W=dPb, mean_delivered_frac=mean_phi, time_frac_capped=frac_capped,
                cap_onset_s_pct=s_onset, cap_onset_t_h=t_onset, mean_delivered_frac_while_capped=mean_phi_capped,
                delivered_frac_at_shutdown=float(phi[-1]), min_Vpk_minus_Vsys_V=margin, exact_check=margin > 0)


def main():
    t_start = time.time()
    apply()
    res = dict(map={}, map_governed={}, governor={}, governor_map={}, S9={})
    comp = []

    def add(k, qty, model, proto, tol, kind):
        if kind == "event":
            comp.append(dict(key=k, qty=qty, model=model, proto=proto, within_tol=model == proto))
        elif kind == "rel":
            d = 100 * (model / proto - 1)
            comp.append(dict(key=k, qty=qty, model=model, proto=proto, diff_pct=d, tol_pct=100 * tol,
                             within_tol=abs(d) <= 100 * tol))
        else:
            d = model - proto
            comp.append(dict(key=k, qty=qty, model=model, proto=proto, diff_abs=d, tol_abs=tol,
                             within_tol=abs(d) <= tol))

    # ---- (1) stranded-charge map with the burst envelope (m_b = 1) and with the online O1 governor
    p_gov = replace(RP63, burst=False)  # exact equivalent of the online cap (module docstring)
    for T in TGRID:
        for soh in SOHS:
            k = mk(T, soh)
            _, r = run_s4(T, soh)
            res["map"][k] = rec(r)
            s0, t0 = PROTO_MAP[(T, soh)]
            res["map"][k].update(proto_s_off_pct=s0, proto_tte_h=t0)
            add(f"map|{k}", "s_off_pct", 100 * r.s_off, s0, 0.3, "abs")
            add(f"map|{k}", "tte_h", r.tte_h, t0, 0.01, "rel")
            add(f"map|{k}", "event", r.event, "cutoff" if (T, soh) == (25, 1.0) else "brownout", None, "event")
            sc, rg = run_s4(T, soh, p=p_gov, record=True)
            g = rec(rg)
            g.update(governor_cost(sc, rg))
            g["gain_tte_pct"] = 100 * (rg.tte_h / r.tte_h - 1)
            g["s_off_reduction_pp"] = 100 * (r.s_off - rg.s_off)
            res["map_governed"][k] = g
    # ---- (2) m_b sweep (§12.2 O1 quantification) over the whole map
    for mb in MB:
        res["governor_map"][f"{mb:g}"] = {}
        for T in TGRID:
            for soh in SOHS:
                if mb == 1.0:
                    rr = res["map"][mk(T, soh)]
                    res["governor_map"]["1"][mk(T, soh)] = dict(s_off_pct=rr["s_off_pct"], tte_h=rr["tte_h"],
                                                                 event=rr["event"])
                    continue
                _, r = run_s4(T, soh, p=replace(RP63, m_b=mb))
                res["governor_map"][f"{mb:g}"][mk(T, soh)] = dict(s_off_pct=100 * r.s_off, tte_h=r.tte_h,
                                                                  event=r.event)
    for mb in MB:
        g = dict(res["governor_map"][f"{mb:g}"][mk(-10, 0.8)])
        g["burst_delivered_frac"] = mb
        g["gain_tte_pct"] = 100 * (g["tte_h"] / res["map"][mk(-10, 0.8)]["tte_h"] - 1)
        res["governor"][f"{mb:g}"] = g
        if mb in PROTO_GOV:
            s0, t0 = PROTO_GOV[mb]
            g.update(proto_s_off_pct=s0, proto_tte_h=t0)
            add(f"governor|{mb:g}", "s_off_pct", g["s_off_pct"], s0, 0.3, "abs")
            add(f"governor|{mb:g}", "tte_h", g["tte_h"], t0, 0.01, "rel")
    og = res["map_governed"][mk(-10, 0.8)]
    res["governor"]["online_O1"] = dict(s_off_pct=og["s_off_pct"], tte_h=og["tte_h"], event=og["event"],
                                        burst_delivered_frac=og["mean_delivered_frac"],
                                        gain_tte_pct=og["gain_tte_pct"],
                                        note="cap 0.9*dPb_safe evaluated on the state; identical TTE/s_off to "
                                             "m_b = 0, but delivers the burst amplitude in full until the cap binds")
    # ---- (3) S9 cold hike grid (explicit 30 s BU bursts vs V_cut; TRK envelope 1 W) with / without governor
    for T in TGRID:
        for soh in SOHS:
            k = mk(T, soh)
            sc = scenario("S9", Tamb_C=float(T), soh=float(soh))
            r = sc.run()
            scg = scenario("S9", Tamb_C=float(T), soh=float(soh), p=p_gov)
            rg = scg.run()
            d = rec(r)
            d.update(T0_C=20.0, tte_governed_h=rg.tte_h, s_off_governed_pct=100 * rg.s_off, event_governed=rg.event)
            if (T, soh) in PROTO_S9:
                t0, s0 = PROTO_S9[(T, soh)]
                d.update(proto_tte_h=t0, proto_s_off_pct=s0)
                tol_t, tol_s = (0.02, 1.0) if (T, soh) == (-20, 0.8) else (0.01, 0.3)
                add(f"S9|{k}", "tte_h", r.tte_h, t0, tol_t, "rel")
                add(f"S9|{k}", "s_off_pct", 100 * r.s_off, s0, tol_s, "abs")
                add(f"S9|{k}", "event", r.event, "cutoff", None, "event")
            res["S9"][k] = d

    runtime = time.time() - t_start
    out = dict(meta=dict(spec=SPEC, seed=None, **git_meta(),
                         solver="LSODA (rtol 1e-7, atol [1e-9,1e-7,1e-5], max_step 60 s), model.run",
                         h_s=None, runtime_s=round(runtime, 1),
                         scenario="S10a = S4 social feed (SO, ΔP_b = 3 W) from s0 = 100 %, cold-soaked (T0 = "
                                  "T_amb); S9 from T0 = 20 °C",
                         key_format="'{Tamb degC}|{SOH:.1f}'; governor keys are m_b",
                         governor="online O1: dPb_eff = min(dPb, max(0, 0.9*dPb_safe)); exact equivalent: envelope "
                                  "off (see run_stranded.py docstring); mean_delivered_frac = time-averaged "
                                  "dPb_eff/dPb"),
               **res, proto_comparison=comp)
    RES.mkdir(exist_ok=True)
    (RES / "stranded.json").write_text(json.dumps(clean(out), indent=1))
    fig_map(res)
    table(res)
    bad = [c for c in comp if not c["within_tol"]]
    chk = all(v["exact_check"] for v in res["map_governed"].values())
    print(f"stranded: {len(comp)} PROTO comparisons, {len(bad)} outside tolerance; governor exactness check "
          f"{'passed' if chk else 'FAILED'}; {time.time() - t_start:.1f} s")
    for c in bad:
        print("  OUTSIDE:", c)


# ------------------------------------------------------------------------------------------- figure
def _heat(ax, M, sub, norm, title):
    im = ax.imshow(M, cmap=cmap_seq, norm=norm, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            dark = norm(v) > 0.55
            col = "white" if dark else INK
            ax.text(j, i - 0.12, f"{v:.1f} %", ha="center", va="center", fontsize=7.2, color=col,
                    fontweight="bold")
            ax.text(j, i + 0.24, sub[i][j], ha="center", va="center", fontsize=5.8,
                    color="white" if dark else INK2)
    ax.set_xticks(range(len(SOHS)), [f"{s:.1f}" for s in SOHS])
    ax.set_yticks(range(len(TGRID)), [f"{t} °C" for t in TGRID])
    ax.set_xlabel("state of health")
    ax.tick_params(length=0)
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(title, fontsize=8, loc="left", color=INK, pad=4)
    return im


def fig_map(res):
    fig = plt.figure(figsize=(W1, 2.75))
    gs = fig.add_gridspec(1, 5, width_ratios=[1, 1, 0.055, 0.22, 1.35], wspace=0.22)
    a, b, cax, c = (fig.add_subplot(gs[0, i]) for i in (0, 1, 2, 4))
    norm = Normalize(0, 25)
    M = np.array([[res["map"][mk(T, s)]["s_off_pct"] for s in SOHS] for T in TGRID])
    sub = [[f"{res['map'][mk(T, s)]['tte_h']:.1f} h" for s in SOHS] for T in TGRID]
    _heat(a, M, sub, norm, "no governor")
    Mg = np.array([[res["map_governed"][mk(T, s)]["s_off_pct"] for s in SOHS] for T in TGRID])
    subg = [[f"{res['map_governed'][mk(T, s)]['tte_h']:.1f} h" for s in SOHS] for T in TGRID]
    im = _heat(b, Mg, subg, norm, "online O1 governor")
    b.set_yticklabels([])
    cb = fig.colorbar(im, cax=cax)
    cax.set_title("$s_{off}$ (%)", fontsize=7, color=INK2, pad=4)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=6.8, length=2)
    panel_label(a, "(a)")
    panel_label(b, "(b)")
    # (c) trade-off at SOH 0.8: stranded charge vs delivered burst amplitude; fixed derating m_b (lines) vs the
    # online O1 cap (diamonds).  TTE as direct labels on the -10 °C series ([PROTO] case of §12.2).
    cols = {-20: SEQ_BLUE[11], -10: SEQ_BLUE[6]}
    for T in (-20, -10):
        y = [res["governor_map"][f"{mb:g}"][mk(T, 0.8)]["s_off_pct"] for mb in MB]
        c.plot(MB, y, "-o", color=cols[T], lw=1.5, ms=3.0, label=f"{T:+d} °C, fixed $m_b$".replace("-", "−"))
        og = res["map_governed"][mk(T, 0.8)]
        c.plot([og["mean_delivered_frac"]], [og["s_off_pct"]], marker="D", ms=4.8, color=cols[T], mec=INK,
               mew=0.6, ls="none")
    c.plot([], [], marker="D", ms=4.8, color="white", mec=INK, mew=0.8, ls="none", label="online O1 cap")
    for mb, off, ha, va in ((1.0, (-4, 2), "right", "bottom"), (0.5, (3, -5), "left", "top")):
        g = res["governor"][f"{mb:g}"]
        c.annotate(f"{g['tte_h']:.2f} h", (mb, g["s_off_pct"]), xytext=off, textcoords="offset points",
                   ha=ha, va=va, fontsize=6.2, color=INK2)
    og = res["governor"]["online_O1"]
    c.annotate(f"{og['tte_h']:.2f} h", (og["burst_delivered_frac"], og["s_off_pct"]), xytext=(-3, 11),
               textcoords="offset points", ha="center", va="bottom", fontsize=6.2, color=INK2)
    c.set_xlabel("delivered burst amplitude (fraction)")
    c.set_ylabel("$s_{off}$ (%) at SOH 0.8")
    c.set_xlim(-0.04, 1.06)
    c.set_ylim(0, 25)
    c.legend(loc="upper left", fontsize=6.3, handlelength=1.4)
    c.set_title("derating vs online cap", fontsize=8, loc="left", color=INK, pad=4)
    panel_label(c, "(c)")
    save(fig, "fig_stranded_map")


def table(res):
    L = [r"\begin{tabular}{lrrrrrr}", r"\toprule",
         r" & \multicolumn{3}{c}{no governor} & \multicolumn{3}{c}{O1 governor} \\",
         r"\cmidrule(lr){2-4}\cmidrule(lr){5-7}",
         r"$T_\mathrm{amb}$ & SOH 1.0 & SOH 0.9 & SOH 0.8 & SOH 1.0 & SOH 0.9 & SOH 0.8 \\", r"\midrule"]
    for T in TGRID:
        cells = []
        for key in ("map", "map_governed"):
            for s in SOHS:
                r = res[key][mk(T, s)]
                cells.append(f"{r['s_off_pct']:.1f} ({r['tte_h']:.2f})")
        L.append(f"${T}\\,^\\circ$C & " + " & ".join(cells) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}"]
    TAB.mkdir(parents=True, exist_ok=True)
    (TAB / "tab_stranded.tex").write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
