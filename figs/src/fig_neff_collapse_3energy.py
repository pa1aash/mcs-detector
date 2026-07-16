#!/usr/bin/env python3
"""The three-energy N_eff^-1 collapse with an explicit closure panel.

The derived C=N_eff^-1 line is a reference behind the data; the pooled OLS fit is
secondary.  Energy points receive only the disclosed +/-4% display offset.  Two
QC categories remain visible and semantically distinct: window-clipped Voronoi
points carry an x overlay, while near-zero diamond controls remain grey diamonds.
For diamond intervals that cross zero, a downward arrow replaces the impossible
negative part of the interval on the log/positive axes.
"""
from __future__ import annotations

import numpy as np
import figstyle as fs

fs.set_style()
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter

L = 10.0
J2 = fs.load("e2_results_combined.json")
J1 = fs.load("e2_results_m3_1000.json")
rows = J2["rows"] + J1["rows"]
aE = {int(k): v for k, v in {**J2["a_eff_of_E"], **J1["a_eff_of_E"]}.items()}
fit = J2["exponent_fit"]

EGREY = "#3a3a3a"
TOP_YMIN, TOP_YMAX = 0.0045, 1.35
RES_YMIN, RES_YMAX = 0.05, 1.72
XMIN, XMAX = 1.7, 23.0
EDODGE = {200: 0.96, 500: 1.0, 1000: 1.04}


def c_value(row):
    return row["dk4"] / (
        3.0 * aE[int(row["E"])] ** 2
        * row["f_built"] * (1.0 - row["f_built"]) * L ** 2
    )


def draw_point(ax, x, y, lo, hi, row, floor):
    """Draw one point, preserving the QC state and zero-crossing information."""
    topo = row["topology"]
    color = fs.COLORS[topo]
    is_diamond = topo == "diamond"
    is_window = bool(row.get("thin")) and not is_diamond
    marker_size = 5.0 if not is_diamond else 5.2
    zorder = 5 if not (is_diamond or is_window) else 4

    if lo > 0:
        lower = y - lo
        upper = hi - y
        ax.errorbar(
            x, y, yerr=[[max(lower, 0)], [max(upper, 0)]],
            fmt=fs.MARKERS[topo], ms=marker_size, mew=1.0, zorder=zorder,
            **fs.estyle(row["E"], color), **fs.ebar(color),
        )
    else:
        # Only the positive upper half can be represented on these axes.  The
        # downward arrow explicitly records that the 95% interval crosses zero.
        ax.errorbar(
            x, y, yerr=[[0], [max(hi - y, 0)]],
            fmt=fs.MARKERS[topo], ms=marker_size, mew=1.0, zorder=zorder,
            **fs.estyle(row["E"], color), **fs.ebar(color),
        )
        start = max(floor * 1.22, y * 0.78)
        ax.annotate(
            "", xy=(x, floor * 1.04), xytext=(x, start),
            arrowprops=dict(
                arrowstyle="-|>", color=color, lw=0.8,
                shrinkA=0, shrinkB=0, mutation_scale=6
            ), zorder=3,
        )

    if is_window:
        ax.plot(x, y, marker="x", ls="none", ms=5.2, mew=1.0,
                color=fs.HERO, zorder=6)


def main():
    positive = [
        r for r in rows
        if r.get("dk4", 0) > 0
        and np.isfinite(r.get("N_eff", np.nan))
        and r["N_eff"] > 0
    ]
    shown = [r for r in positive if XMIN <= r["N_eff"] <= XMAX]

    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(fs.FULL, fs.FULL * 0.78), sharex=True,
        height_ratios=[2.45, 1.05]
    )

    # Reference and fit are laid down first, behind every observation.
    xx = np.logspace(np.log10(XMIN), np.log10(XMAX), 160)
    ax.plot(xx, 1.0 / xx, "-", color=fs.HERO, lw=2.0, zorder=1)
    ax.plot(
        xx, np.exp(fit["intercept"]) * xx ** fit["slope"],
        ls=fs.DASH_FIT, color=fs.MUTE, lw=fs.LW_FIT, zorder=2
    )

    for row in shown:
        C = c_value(row)
        C_lo = row["dk4_lo"] / row["dk4"] * C
        C_hi = row["dk4_hi"] / row["dk4"] * C
        x = row["N_eff"] * EDODGE[int(row["E"])]
        draw_point(ax, x, C, C_lo, C_hi, row, TOP_YMIN)
        draw_point(
            axr, x, C * row["N_eff"], C_lo * row["N_eff"],
            C_hi * row["N_eff"], row, RES_YMIN
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(XMIN, XMAX)
    ax.set_ylim(TOP_YMIN, TOP_YMAX)
    ax.set_yticks([0.01, 0.03, 0.1, 0.3, 1.0])
    ax.set_yticklabels(["0.01", "0.03", "0.1", "0.3", "1"])
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_ylabel(r"normalised excess fourth cumulant $C$")
    fs.despine(ax)
    fs.panel_title(ax, "a", "Three-energy scale collapse")

    # One structured legend: the second column also discloses both QC states.
    header = lambda label: Line2D([], [], marker="none", ls="none",
                                  label=rf"\textbf{{{label}}}")
    blank = Line2D([], [], marker="none", ls="none", label="")
    energy_h = [
        Line2D([], [], marker="o", ls="none", ms=5.8,
               label=rf"\SI{{{E}}}{{MeV}}", **fs.estyle(E, EGREY))
        for E in (200, 500, 1000)
    ]
    topo_h = [
        Line2D([], [], marker=fs.MARKERS[t], ls="none", mfc="white",
               mec=fs.COLORS[t], mew=1.1, ms=5.8, label=fs.TOPO_LABEL[t])
        for t in ("rectilinear", "schwarzp", "gyroid", "voronoi")
    ]
    diamond_h = Line2D(
        [], [], marker=fs.MARKERS["diamond"], ls="none", mfc="white",
        mec=fs.COLORS["diamond"], mew=1.1, ms=5.8,
        label="diamond control (excl.)"
    )
    window_h = Line2D(
        [], [], marker="x", ls="none", color=fs.HERO, mew=1.0, ms=5.8,
        label="window-clipped (excl.)"
    )
    # Seven entries per column (Matplotlib fills columns top-to-bottom).
    handles = [header("Energy"), *energy_h, blank, blank, blank,
               header("Topology / QC"), *topo_h, diamond_h, window_h]
    leg_main = ax.legend(
        handles=handles, ncol=2, loc="upper right", fontsize=8.5,
        columnspacing=1.0, handletextpad=0.42, labelspacing=0.25,
        borderaxespad=0.25,
    )
    ax.add_artist(leg_main)
    # The two reference lines are labelled here (horizontally) rather than with
    # rotated in-plot text.
    line_handles = [
        Line2D([], [], color=fs.HERO, lw=2.0, ls="-",
               label=r"derived $C=N_\mathrm{eff}^{-1}$"),
        Line2D([], [], color=fs.MUTE, lw=fs.LW_FIT, ls=fs.DASH_FIT,
               label=r"pooled fit"),
    ]
    ax.legend(handles=line_handles, loc="lower left", fontsize=8.5,
              frameon=False, handlelength=2.6, labelspacing=0.3,
              borderaxespad=0.6)

    axr.axhspan(0.9, 1.1, color=fs.FAINT, alpha=0.28, lw=0, zorder=0)
    axr.axhline(1.0, color=fs.HERO, lw=1.3, zorder=1)
    axr.text(
        0.012, 1.105, r"$\pm10\%$ reference",
        transform=axr.get_yaxis_transform(), fontsize=8.5,
        color=fs.MUTE, ha="left", va="bottom"
    )
    axr.set_xscale("log")
    axr.set_xticks([2, 3, 5, 10, 20])
    axr.set_xticklabels(["2", "3", "5", "10", "20"])
    axr.minorticks_off()
    axr.set_ylim(RES_YMIN, RES_YMAX)
    axr.set_yticks([0.1, 0.5, 1.0, 1.5])
    axr.set_xlabel(r"effective cell count $N_\mathrm{eff}=f(1-f)L^2/\mathrm{Var}(t)$")
    axr.set_ylabel(r"$C\,N_\mathrm{eff}$")
    fs.despine(axr)
    fs.panel_title(axr, "b", "Closure about the derived law")

    fs.save(fig, "fig_neff_collapse_3energy")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
