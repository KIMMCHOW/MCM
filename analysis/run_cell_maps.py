"""AN2: cell maps (FINAL_SPEC §2.5, §2.6, §2.9, F4; §14).

Outputs
  results/cell_maps.json   ocv_table (U(z,T)), R0_grid_mOhm, R1_mOhm/tau1_s vs T, Pcut_W (cut-off-limited power
                           V_cut(U − V_cut)/R_Σ), Pcut_sys_W (burst-limited power V_sys(U − V_sys)/R0), each for the
                           fresh pack and SOH 0.8 (r_growth = γ_R·0.2), plus operating_points (L2 LSODA shutdowns
                           of S4, S6, S7) and the F4 brown-out roots for the same conditions.
  paper/figures/fig_cell_maps.pdf
      (a) OCV U(z,T) and the quasi-static loaded voltage at the S4 load (1.06 W) for six temperatures;
      (b) R0(z,T);  (c) P_cut,sys(z,T) map with the brown-out boundaries P + ΔP_b = P_cut,sys of S4, S6, S7 and
      the L2 shutdown points (z, T) of those scenarios at several ambients.

Maps use v1 = 0 (Ũ = U), as in the F4 P_cut table of §7; the L2 shutdown happens with v1 > 0, so the brown-out
points sit slightly to the right of the v1 = 0 boundary.  T in the maps is the *battery* temperature.
"""
from __future__ import annotations

import time

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D

from basics_common import f, meta, write_json
from phonebatt import RP63, aged, analytic, cell, power, scenario
from phonebatt.params import K0
from style import GRID, INK, INK2, MUTED, SEQ_BLUE, W1, apply, cmap_seq, panel_label, save

T_LINES = (40, 25, 10, 0, -10, -20)  # °C, warm -> cold
T_COLOR = {Tc: SEQ_BLUE[i] for Tc, i in zip(T_LINES, (3, 5, 7, 9, 11, 12))}  # light = warm, dark = cold
SOH_AGED = 0.8
RUNS = [("S4", 25), ("S4", 0), ("S4", -10), ("S4", -20), ("S6", 25), ("S6", 0), ("S6", -10), ("S6", -20),
        ("S7", 35), ("S7", 25), ("S7", 15), ("S7", 0), ("S7", -10)]
PRIM = {"S4": "SO", "S6": "NV", "S7": "GM"}
MARK = {"S4": "o", "S6": "s", "S7": "^"}


def grids(p, z, Tc):
    """U, R0, R1, P_cut, P_cut,sys on the (T, z) grid (rows T, columns z); v1 = 0."""
    T = (np.asarray(Tc, float) + K0)[:, None]
    zz = np.asarray(z, float)[None, :]
    U = cell.ocv(zz, T, p)
    r0 = cell.R0(zz, T, p)
    r1 = cell.R1(T, p) + 0 * zz
    Pc = cell.p_cut(U, r0 + r1, p.V_cut)
    Ps = cell.p_cut(U, r0, p.V_sys)
    return U, r0, r1, Pc, Ps


def main() -> None:
    t0 = time.time()
    pa = aged(RP63, 1 - SOH_AGED)
    # ------------------------------------------------------------------------- JSON grids (coarse)
    zj = np.round(np.linspace(0, 1, 51), 4)
    Tj = np.arange(-20, 46, 5)
    Uf, R0f, R1f, Pcf, Psf = grids(RP63, zj, Tj)
    _, R0a, R1a, Pca, Psa = grids(pa, zj, Tj)
    out = dict(meta=meta("run_cell_maps.py", solver="LSODA (operating points)", h_s=None,
                         note="grids: rows = T_C, columns = z (coulombic SOC); v1 = 0; T is battery temperature; "
                              "negative P_cut values mean U < V* (no deliverable power)"))
    out["ocv_table"] = dict(z=zj, T_C=Tj, U_V=Uf, U_Tref_V=cell.ocv(zj, 298.15, RP63),
                            max_entropic_shift_mV=1e3 * float(np.max(np.abs(Uf - cell.ocv(zj, 298.15, RP63)))))
    out["R0_grid_mOhm"] = dict(z=zj, T_C=Tj, fresh=1e3 * R0f, soh_0p8=1e3 * R0a)
    out["R1_mOhm"] = dict(T_C=Tj, fresh=1e3 * R1f[:, 0], soh_0p8=1e3 * R1a[:, 0])
    out["tau1_s"] = dict(T_C=Tj, value=cell.tau1(Tj + K0, RP63))
    out["Pcut_W"] = dict(z=zj, T_C=Tj, fresh=Pcf, soh_0p8=Pca, V_star="V_cut", R="R0 + R1")
    out["Pcut_sys_W"] = dict(z=zj, T_C=Tj, fresh=Psf, soh_0p8=Psa, V_star="V_sys", R="R0")
    out["z_e"] = f(cell.z_empty(RP63))

    # ------------------------------------------------------------------------- operating points (L2 LSODA)
    ops = []
    for sc, Ta in RUNS:
        r = scenario(sc, Tamb_C=Ta).run()
        u = power.PRIMITIVES[PRIM[sc]]
        dPb = float(u.dPb * RP63.m_b)
        P = r.P_end_W
        Te = r.T_end_C + K0
        Ue = f(cell.ocv(r.z_end, Te, RP63))
        Ps_map = f(cell.p_cut(Ue, cell.R0(r.z_end, Te, RP63), RP63.V_sys))
        br = analytic.brownout_root(P, dPb, Te, RP63)
        ops.append(dict(scenario=sc, Tamb_C=Ta, event=r.event, tte_h=r.tte_h, s_off_pct=100 * r.s_off,
                        z_end=r.z_end, T_end_C=r.T_end_C, P_end_W=P, dPb_W=dPb, P_pk_W=P + dPb,
                        Pcut_sys_map_at_end_W=Ps_map, z_star_F4=f(br["z_star"]), s_off_F4_pct=100 * f(br["s_off"])))
    out["operating_points"] = ops
    # [PROTO] references for the same runs (§5.4 main table, stranded map, navigation note)
    out["proto_refs"] = {"S4|25": [16.64, "cutoff", 0.40], "S4|0": [16.43, "brownout", 1.4],
                         "S4|-10": [16.18, "brownout", 2.8], "S4|-20": [15.54, "brownout", 6.4],
                         "S6|25": [11.83, "cutoff", 0.56], "S6|-10": [11.43, None, None],
                         "S7|35": [10.29, "cutoff", 0.57], "S7|25": [5.68, "cutoff", 1.07],
                         "S7|15": [4.40, "brownout", 1.89]}
    devs = {}
    for o in ops:
        k = f"{o['scenario']}|{o['Tamb_C']}"
        if k in out["proto_refs"]:
            t, ev, so = out["proto_refs"][k]
            devs[k] = dict(tte_dev_pct=100 * (o["tte_h"] / t - 1), event_match=None if ev is None else ev == o["event"],
                           s_off_dev_pp=None if so is None else o["s_off_pct"] - so)
    out["proto_check"] = devs
    write_json("cell_maps.json", out)

    figure(pa, ops)
    print(f"max entropic OCV shift over -20..45 C: {out['ocv_table']['max_entropic_shift_mV']:.1f} mV")
    for o in ops:
        print(f"  {o['scenario']} {o['Tamb_C']:>4} C: {o['event']:9s} TTE {o['tte_h']:6.2f} h  s_off {o['s_off_pct']:5.2f} %"
              f"  z_end {o['z_end']:.4f} (F4 z* {o['z_star_F4']:.4f})  T_end {o['T_end_C']:5.1f}  P_pk {o['P_pk_W']:.2f} W")
    for k, d in devs.items():
        print(f"  vs PROTO {k}: {d}")
    print(f"done in {time.time() - t0:.1f} s")


def figure(pa, ops) -> None:
    apply()
    fig = plt.figure(figsize=(W1, 2.55))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.25], wspace=0.5, left=0.075, right=0.9, top=0.76, bottom=0.16)
    ax1, ax2, ax3 = (fig.add_subplot(gs[i]) for i in range(3))
    cax = ax3.inset_axes([1.04, 0.0, 0.05, 1.0])
    z = np.linspace(0.0, 1.0, 401)

    # (a) OCV and loaded voltage --------------------------------------------------------------------
    P_S4 = 1.06
    for Tc in T_LINES:
        T = Tc + K0
        ax1.plot(z, cell.ocv(z, T, RP63), color=T_COLOR[Tc], lw=1.1)
        ax1.plot(z, analytic.v_qs(z, P_S4, T, RP63), color=T_COLOR[Tc], lw=1.1, ls=(0, (3, 1.5)))
    ax1.axhline(RP63.V_cut, color=MUTED, lw=0.8, ls=(0, (1, 1.5)))
    ax1.axhline(f(RP63.V_sys), color=MUTED, lw=0.8, ls=(0, (1, 1.5)))
    ax1.text(0.98, RP63.V_cut + 0.012, "$V_\\mathrm{cut}$", ha="right", va="bottom", fontsize=6.5, color=INK2)
    ax1.text(0.98, f(RP63.V_sys) - 0.012, "$V_\\mathrm{sys}$", ha="right", va="top", fontsize=6.5, color=INK2)
    ax1.text(0.40, 3.93, "OCV $U(z,T)$", fontsize=6.5, color=INK2, ha="right")
    ax1.annotate("loaded, 1.06 W", xy=(0.55, f(analytic.v_qs(0.55, P_S4, 253.15, RP63))), xytext=(0.62, 3.58),
                 fontsize=6.5, color=INK2, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))
    ax1.set_xlim(0, 1)
    ax1.set_ylim(3.2, 4.5)
    ax1.set_xlabel("coulombic SOC $z$")
    ax1.set_ylabel("voltage (V)")
    panel_label(ax1, "(a)")

    # (b) R0 ---------------------------------------------------------------------------------------
    for Tc in T_LINES:
        ax2.plot(z, 1e3 * cell.R0(z, Tc + K0, RP63), color=T_COLOR[Tc], lw=1.3)
    ax2.plot(z, 1e3 * cell.R0(z, -20 + K0, pa), color=T_COLOR[-20], lw=1.1, ls=(0, (3, 1.5)))
    ax2.set_yscale("log")
    ax2.set_yticks([40, 60, 100, 200, 400])
    ax2.set_yticklabels(["40", "60", "100", "200", "400"])
    ax2.minorticks_off()
    ax2.set_ylim(35, 650)
    ax2.set_xlim(0, 1)
    ax2.text(0.97, 1e3 * f(cell.R0(0.97, -20 + K0, pa)) * 1.07, "−20 °C, SOH 0.8", ha="right", va="bottom",
             fontsize=6.3, color=INK2)
    ax2.set_xlabel("coulombic SOC $z$")
    ax2.set_ylabel("$R_0$ (mΩ)")
    panel_label(ax2, "(b)")

    # (c) P_cut,sys map ------------------------------------------------------------------------------
    zm = np.geomspace(0.02, 1.0, 241)
    Tm = np.linspace(-20, 45, 131)
    _, _, _, _, Ps = grids(RP63, zm, Tm)
    Ps = np.where(Ps > 0.5, Ps, 0.5)
    levels = np.geomspace(0.5, 64, 22)
    cf = ax3.contourf(zm, Tm, Ps, levels=levels, cmap=cmap_seq, norm=LogNorm(0.5, 64), extend="max")
    cb = fig.colorbar(cf, cax=cax, ticks=[1, 2, 4, 8, 16, 32, 64])
    cb.ax.set_yticklabels(["1", "2", "4", "8", "16", "32", "64"])
    cb.ax.tick_params(labelsize=6.5, length=2)
    cb.set_label("$P_\\mathrm{cut,sys}$ (W)", fontsize=7.5, color=INK2, labelpad=2)
    cb.outline.set_visible(False)
    peaks = {}
    for o in ops:  # peak demand used for each boundary: the run's own end power + burst
        peaks.setdefault(o["scenario"], o["P_pk_W"])
    peaks["S7"] = max(o["P_pk_W"] for o in ops if o["scenario"] == "S7")  # unthrottled (cool ambient)
    styles = {"S4": (0, (4, 1.5)), "S6": (0, (1.2, 1.2)), "S7": "solid"}
    for sc, lev in peaks.items():
        ax3.contour(zm, Tm, Ps, levels=[lev], colors=[INK], linewidths=1.0, linestyles=[styles[sc]])
    for o in ops:
        brown = o["event"] == "brownout"
        ax3.plot(o["z_end"], o["T_end_C"], marker=MARK[o["scenario"]], ms=4.0, mec=INK, mew=0.8,
                 mfc="white" if brown else INK, ls="none", zorder=5)
    ax3.set_xscale("log")
    ax3.set_xlim(0.02, 1.0)
    ax3.set_xticks([0.02, 0.05, 0.1, 0.2, 0.5, 1.0])
    ax3.set_xticklabels(["0.02", "0.05", "0.1", "0.2", "0.5", "1"])
    ax3.minorticks_off()
    ax3.set_ylim(-20, 45)
    ax3.set_xlabel("coulombic SOC $z$ (log scale)")
    ax3.set_ylabel("battery temperature (°C)", labelpad=1)
    ax3.grid(False)
    panel_label(ax3, "(c)")
    ax3.set_title("burst-limited power, fresh pack", loc="left", fontsize=7.5, fontweight="normal", color=INK2,
                  pad=3)
    # contour labels placed at the bottom of each boundary
    # shared legends ----------------------------------------------------------------------------------
    hT = [Line2D([], [], color=T_COLOR[Tc], lw=1.5, label=f"{Tc} °C") for Tc in (40, 0, 25, -10, 10, -20)]
    fig.legend(handles=hT, loc="upper left", ncol=3, bbox_to_anchor=(0.06, 1.0), fontsize=6.6, handlelength=1.3,
               columnspacing=0.9, handletextpad=0.3, title=None)
    hS = [Line2D([], [], color=INK, lw=0.9, ls=styles[s], marker=MARK[s], ms=4, mec=INK, mfc="white",
                 label=f"{s} ({peaks[s]:.1f} W)") for s in ("S4", "S6", "S7")]
    hS += [Line2D([], [], marker="o", ms=4, mec=INK, mfc="white", ls="none", label="brown-out"),
           Line2D([], [], marker="o", ms=4, mec=INK, mfc=INK, ls="none", label="cut-off")]
    hS = [hS[0], hS[3], hS[1], hS[4], hS[2]]  # column-major fill -> rows: S4 S6 S7 / brown-out cut-off
    fig.legend(handles=hS, loc="upper right", ncol=3, bbox_to_anchor=(0.935, 1.0), fontsize=6.4, handlelength=2.0,
               columnspacing=0.8, handletextpad=0.25)
    save(fig, "fig_cell_maps")


if __name__ == "__main__":
    main()
