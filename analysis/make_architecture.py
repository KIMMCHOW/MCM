"""AN16: model architecture block diagram (FINAL_SPEC §1, §2; §14).

Output: paper/figures/fig_architecture.pdf (+ PNG preview).  Pure matplotlib patches; no data.

Flow: usage u(t) -> component load model -> constant-power ECM (OCV, R0, RC) <-> lumped thermal node; skin
temperature -> throttle feedback to the SoC power; burst envelope -> peak voltage -> brown-out; slow aging states
-> Q_max, R; outputs SOC / TTE / events.  Feedback arrows (T -> R, T_s -> P_SoC, SOH -> Q, R) are drawn in the
accent colour.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from style import AXIS, INK, INK2, NEUTRAL, ORANGE, VIOLET, W1, apply, save

FB = ORANGE  # electro-thermal feedback
AG = VIOLET  # slow aging coupling


def box(ax, x, y, w, h, title, lines=(), fc=NEUTRAL, ec=AXIS, title_size=7.4, size=6.5, dy=2.55, top_pad=2.3):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.25,rounding_size=1.0", fc=fc, ec=ec, lw=0.8,
                                zorder=1))
    ax.text(x + w / 2, y + h - top_pad, title, ha="center", va="center", fontsize=title_size, fontweight="bold",
            color=INK, zorder=3)
    yy = y + h - top_pad - 3.0
    for ln in lines:
        ax.text(x + w / 2, yy, ln, ha="center", va="center", fontsize=size, color=INK2, zorder=3)
        yy -= dy


def arrow(ax, p0, p1, color=INK2, ls="solid", lw=1.0, rad=0.0, label=None, lpos=None, lha="center", lva="center",
          lsize=6.3, lcolor=INK2):
    a = FancyArrowPatch(p0, p1, arrowstyle="-|>,head_length=3.2,head_width=1.8", color=color, lw=lw, ls=ls,
                        connectionstyle=f"arc3,rad={rad}", shrinkA=0, shrinkB=0, zorder=2, mutation_scale=1.0)
    ax.add_patch(a)
    if label:
        x, y = lpos if lpos else ((p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2)
        ax.text(x, y, label, ha=lha, va=lva, fontsize=lsize, color=lcolor, zorder=4,
                bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none"))


def poly_arrow(ax, pts, color=INK2, ls="solid", lw=1.0):
    """Orthogonal connector through ``pts`` with an arrow head on the last segment."""
    xs, ys = zip(*pts[:-1])
    ax.plot(xs, ys, color=color, lw=lw, ls=ls, solid_capstyle="butt", zorder=2)
    arrow(ax, pts[-2], pts[-1], color=color, ls=ls, lw=lw)


def main() -> None:
    apply()
    fig = plt.figure(figsize=(W1, 3.0))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 44.8)
    ax.axis("off")

    # ---------------------------------------------------------------------------------------------- blocks
    box(ax, 1.0, 13.5, 15.5, 19.0, "Usage input $\\boldsymbol{u}(t)$",
        ["deterministic schedule", "of primitives, or a", "usage CTMC path", "(screen, $b$, APL, $x_c$,",
         "radio, GNSS, camera)"])
    box(ax, 20.5, 13.5, 21.0, 19.0, "Component load model",
        ["OLED: $P_{d0}$ + area·$b^{\\gamma}$·APL", "SoC: DVFS, superlinear in $x_c$",
         "radios: exact Poisson tails", "floor, GNSS, camera, audio", "$P=\\Sigma_k P_k$ at the battery"])
    box(ax, 46.0, 13.5, 22.5, 19.0, "Constant-power ECM",
        ["states $z,\\ v_1$", "OCV $U_{oc}(z,T)$, $R_0(z,T)$, $R_1$, $\\tau_1$",
         "$I=2P/(\\tilde U+\\sqrt{D})$, $V=\\tilde U-R_0I$", "burst peak $V_{pk}$ at $P+\\Delta P_b$",
         "($v_1$ frozen, sub-second)"])
    box(ax, 46.0, 1.0, 22.5, 8.6, "Lumped thermal node $T$",
        ["Bernardi heat + Newton cooling"], top_pad=2.2)
    box(ax, 20.5, 1.0, 21.0, 8.6, "Skin $T_s$ and throttle",
        ["hot spot $T_s$, logistic cap $\\sigma(T_s)$"], top_pad=2.2)
    box(ax, 46.0, 36.0, 22.5, 7.8, "Slow aging states (L3)",
        ["$q_{cal},\\ q_H,\\ q_L$ $\\rightarrow$ SOH, $r_g$"], top_pad=2.0)
    box(ax, 74.0, 1.0, 25.0, 42.8, "Outputs",
        ["SOC: coulombic $s(t)$,", "energy fraction $\\varsigma(t)$", "",
         "TTE = first event:", "charge  $z\\leq 0$", "cut-off  $V\\leq V_{cut}$", "brown-out  $V_{pk}\\leq V_{sys}$",
         "thermal  $T\\geq T_{shut}$", "infeasible  $P>P_{max}$", "",
         "stranded charge $s_{off}$", "$T,\\ T_s,\\ \\rho_{perf},\\ P_k(t)$", "",
         "TTE distribution (CTMC)", "SOH, TTE vs. age (L3)"], dy=2.5, top_pad=2.5)

    # ---------------------------------------------------------------------------------------------- flow
    arrow(ax, (16.9, 23.0), (20.1, 23.0))
    arrow(ax, (41.9, 26.0), (45.6, 26.0), label="$P(t)$", lpos=(43.75, 28.0), lsize=6.5)
    arrow(ax, (41.9, 17.0), (45.6, 17.0), label="$\\Delta P_b$", lpos=(43.75, 19.0), lsize=6.5)
    arrow(ax, (68.9, 26.0), (73.6, 26.0), label="$z, V$", lpos=(71.25, 28.0), lsize=6.5)
    arrow(ax, (68.9, 17.0), (73.6, 17.0), label="$V_{pk}$", lpos=(71.25, 19.0), lsize=6.5)
    arrow(ax, (68.9, 5.3), (73.6, 5.3), label="$T$", lpos=(71.25, 7.2), lsize=6.5)

    # electro-thermal coupling
    arrow(ax, (52.5, 13.1), (52.5, 10.0), label="heat $q$", lpos=(51.6, 11.55), lha="right", lsize=6.3)
    arrow(ax, (62.0, 10.0), (62.0, 13.1), color=FB, lw=1.2, label="$T\\rightarrow R_0, R_1, \\tau_1$",
          lpos=(62.9, 11.55), lha="left", lsize=6.3)
    # thermal -> skin -> throttle -> SoC power
    arrow(ax, (45.6, 5.3), (41.9, 5.3), color=FB, lw=1.2, label="$T$", lpos=(43.75, 7.2), lsize=6.5)
    arrow(ax, (31.0, 10.0), (31.0, 13.1), color=FB, lw=1.2, label="$P_{SoC}$ capped", lpos=(31.9, 11.55),
          lha="left", lsize=6.3)
    # aging coupling
    arrow(ax, (52.5, 35.6), (52.5, 32.9), color=AG, lw=1.2, label="$Q_{max}\\downarrow,\\ R\\uparrow$",
          lpos=(51.6, 34.25), lha="right", lsize=6.3)
    arrow(ax, (62.0, 32.9), (62.0, 35.6), color=AG, lw=1.0, ls=(0, (2.5, 1.5)), label="$I, T, z$ over days",
          lpos=(62.9, 34.25), lha="left", lsize=6.3)

    # legend for the feedback colours (text in ink)
    ax.plot([2.0, 5.0], [41.5, 41.5], color=FB, lw=1.2)
    ax.text(5.8, 41.5, "electro-thermal feedback", fontsize=6.3, color=INK2, va="center")
    ax.plot([2.0, 5.0], [38.3, 38.3], color=AG, lw=1.2)
    ax.text(5.8, 38.3, "aging coupling (slow)", fontsize=6.3, color=INK2, va="center")
    ax.plot([2.0, 5.0], [35.1, 35.1], color=INK2, lw=1.0)
    ax.text(5.8, 35.1, "signal flow", fontsize=6.3, color=INK2, va="center")
    save(fig, "fig_architecture")
    print("fig_architecture written")


if __name__ == "__main__":
    main()
