"""Shared figure style for all paper figures (static PDF, light surface).

Rules (dataviz method):
* Categorical hues are assigned in the fixed order of ``CAT`` and follow the
  entity, never its rank.  Never cycle past 8; fold extra series into "Other"
  or use small multiples.  Scatter/bubble/small-multiple forms use at most the
  first 3 slots.
* Sequential magnitude = one hue light->dark (``SEQ_BLUE`` / ``cmap_seq``).
* Diverging polarity = blue <-> red with a neutral gray midpoint (``cmap_div``).
* One y-axis per panel - never a twin/dual axis.  Different measures -> separate
  panels.
* Thin marks (lines 1.6-2 pt), recessive grid, legend for >= 2 series,
  selective direct labels only.
* Text is in ink colors, never in series colors.

Usage::

    from style import apply, CAT, save
    apply()
    fig, ax = plt.subplots(figsize=(W1, 2.6))
    ...
    save(fig, "fig_name")      # -> paper/figures/fig_name.pdf (+ .png preview)
"""
from __future__ import annotations

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "paper" / "figures"

# Categorical slots (fixed order; validated light-mode palette)
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
       "#e87ba4", "#008300", "#6250d6", "#e34948"]
BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED = CAT

# Ink and chrome
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SURFACE = "#fcfcfb"
NEUTRAL = "#f0efec"

# Status colors (reserved; always with a label)
GOOD, WARNING, SERIOUS, CRITICAL = "#0ca30c", "#fab219", "#ec835a", "#d03b3b"

SEQ_BLUE = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
            "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
SEQ_ORANGE = ["#fde3d6", "#f9c4ab", "#f4a47f", "#ef8558", "#eb6834", "#c9531f", "#a34015", "#7d300f"]
RED_ARM = ["#f7d4d3", "#f0aaa9", "#e98180", "#e34948", "#bd3534", "#962726"]

cmap_seq = LinearSegmentedColormap.from_list("seq_blue", SEQ_BLUE)
cmap_seq_orange = LinearSegmentedColormap.from_list("seq_orange", SEQ_ORANGE)
cmap_div = LinearSegmentedColormap.from_list(
    "div_blue_red", list(reversed(SEQ_BLUE[2:11:2])) + [NEUTRAL] + RED_ARM[:5])

# Figure widths (inches) for the LaTeX layout (text width ~ 6.7 in)
W1 = 6.7     # full width
W_HALF = 3.25


def apply() -> None:
    """Apply the global rcParams."""
    mpl.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": 8.5,
        "axes.titlesize": 9,
        "axes.titleweight": "bold",
        "axes.labelsize": 8.5,
        "axes.labelcolor": INK2,
        "axes.edgecolor": AXIS,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "axes.axisbelow": True,
        "grid.color": GRID,
        "grid.linewidth": 0.5,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelcolor": INK2,
        "ytick.labelcolor": INK2,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 7.5,
        "legend.frameon": False,
        "lines.linewidth": 1.8,
        "lines.markersize": 4.5,
        "axes.prop_cycle": mpl.cycler(color=CAT),
        "text.color": INK,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.03,
        "figure.dpi": 150,
    })


def save(fig, name: str, png: bool = True) -> Path:
    """Save ``fig`` to paper/figures/<name>.pdf (and a PNG preview)."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    out = FIG_DIR / f"{name}.pdf"
    fig.savefig(out)
    if png:
        prev = FIG_DIR / "png"
        prev.mkdir(exist_ok=True)
        fig.savefig(prev / f"{name}.png", dpi=110)
    plt.close(fig)
    return out


def panel_label(ax, text: str) -> None:
    """Bold panel tag (a), (b) ... at the top-left, in ink color."""
    ax.text(-0.02, 1.04, text, transform=ax.transAxes, fontsize=9,
            fontweight="bold", color=INK, ha="right", va="bottom")
