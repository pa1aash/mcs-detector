#!/usr/bin/env python3
"""The homogeneous approximation is half-valid: same width, different shape.

Panel (a) shows the width-equivalent material fraction after normalising by each
run's realised mean PLA path.  This removes the nominal/as-built fill mismatch and
tests the actual width law; every campaign configuration closes within 1.9%.

Panel (b) tests the complementary shape law directly, plotting the observed
deconvolved Delta-kappa4 against 3 a_eff^2 Var(t_pla) for every positive,
QC-resolved campaign row at 200, 500 and 1000 MeV.
"""
from __future__ import annotations

import collections
import numpy as np
import figstyle as fs

fs.set_style()
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

J2 = fs.load("e2_results_combined.json")
J1 = fs.load("e2_results_m3_1000.json")
rows = J2["rows"] + J1["rows"]
a_eff = {int(k): v for k, v in {**J2["a_eff_of_E"], **J1["a_eff_of_E"]}.items()}

TOPOS = ["rectilinear", "schwarzp", "gyroid", "diamond", "voronoi"]
FILLS = [0.20, 0.30, 0.40, 0.50]
DIAMOND_C = "#565656"


def topo_color(topo):
    return DIAMOND_C if topo == "diamond" else fs.COLORS[topo]


def corrected_width(row):
    """Width-equivalent fill divided by the realised mean path fraction."""
    return row["f_width"] / (row["tpla_mean"] / 10.0)


def main():
    width_rows = [
        r for r in rows
        if r.get("f_designed", 0) > 0
        and np.isfinite(r.get("f_width", np.nan))
        and np.isfinite(r.get("tpla_mean", np.nan))
        and r["tpla_mean"] > 0
    ]
    shape_rows = [
        r for r in rows
        if r.get("dk4", 0) > 0
        and not r.get("thin")
        and np.isfinite(r.get("tpla_var", np.nan))
        and r["tpla_var"] > 0
    ]

    fig, (axA, axB) = plt.subplots(
        2, 1, figsize=(fs.FULL, fs.FULL * 0.90),
        height_ratios=[0.90, 1.45]
    )

    # ---------------- (a) width after realised-path correction ----------------
    S = 0.52
    FS = 1.82
    GGAP = 2.05
    GBASE = (len(FILLS) - 1) * FS + GGAP
    marker_size = 5.2

    axA.axhspan(0.98, 1.02, color=fs.FAINT, alpha=0.30, lw=0, zorder=0)
    axA.axhline(1.0, color=fs.HERO, lw=1.1, zorder=1)

    clusters = collections.defaultdict(list)
    for row in width_rows:
        clusters[(row["topology"], round(row["f_designed"], 2))].append(row)

    point_x = {}
    for (topo, fill), group in clusters.items():
        ti, fi = TOPOS.index(topo), FILLS.index(fill)
        centre = ti * GBASE + fi * FS
        group = sorted(group, key=lambda r: int(r["E"]))
        for k, row in enumerate(group):
            x = centre + (k - (len(group) - 1) / 2) * S
            point_x[(topo, fill, int(row["E"]))] = x
            axA.plot(
                x, corrected_width(row), fs.MARKERS[topo],
                ms=marker_size, mew=1.0, zorder=4,
                **fs.estyle(row["E"], topo_color(topo)),
            )

    for ti, topo in enumerate(TOPOS):
        if ti:
            axA.axvline(ti * GBASE - GGAP / 2, color=fs.FAINT, lw=0.6, zorder=0)
        axA.text(
            ti * GBASE + 1.5 * FS, 1.0217, fs.TOPO_LABEL[topo],
            ha="center", va="bottom", fontsize=9.5, color=topo_color(topo)
        )

    # The +/-2% band plus the caption's exact "within 1.9%" carry the closure
    # bound; no per-point deviation label is needed on the canvas.
    axA.text(
        0.012, 1.0158, r"$\pm2\%$ reference", transform=axA.get_yaxis_transform(),
        fontsize=8.5, color=fs.MUTE, ha="left", va="center"
    )

    tick_x = [ti * GBASE + fi * FS for ti in range(len(TOPOS)) for fi in range(len(FILLS))]
    axA.set_xticks(tick_x)
    axA.set_xticklabels([f"{f:.1f}" for _ in TOPOS for f in FILLS], fontsize=8.5)
    axA.set_xlim(-1.05, (len(TOPOS) - 1) * GBASE + 3 * FS + 1.05)
    axA.set_ylim(0.976, 1.026)
    axA.set_yticks([0.98, 0.99, 1.00, 1.01, 1.02])
    axA.set_xlabel(r"designed fill fraction $f$")
    axA.set_ylabel(r"$f_\mathrm{width}/(\langle t_\mathrm{pla}\rangle/L)$")
    fs.despine(axA)
    fs.panel_title(axA, "a", "Width after mean-path correction")

    # ------------------ (b) fourth-cumulant shape closure -------------------
    predictions = np.array([
        3.0 * a_eff[int(r["E"])] ** 2 * r["tpla_var"] for r in shape_rows
    ])
    observed = np.array([r["dk4"] for r in shape_rows])
    lo_lim = min(predictions.min(), observed.min()) / 1.45
    hi_lim = max(predictions.max(), observed.max()) * 1.45
    axB.plot([lo_lim, hi_lim], [lo_lim, hi_lim], color=fs.HERO,
             lw=1.5, zorder=1)
    axB.text(0.79, 0.82, r"$y=x$", transform=axB.transAxes,
             color=fs.MUTE, fontsize=9, rotation=36, ha="center", va="center")

    for row, pred in zip(shape_rows, predictions):
        topo = row["topology"]
        color = topo_color(topo)
        y = row["dk4"]
        ylo, yhi = row["dk4_lo"], row["dk4_hi"]
        axB.errorbar(
            pred, y,
            yerr=[[y - ylo], [yhi - y]],
            fmt=fs.MARKERS[topo], ms=5.2, mew=1.0, zorder=4,
            **fs.estyle(row["E"], color), **fs.ebar(color),
        )

    axB.set_xscale("log")
    axB.set_yscale("log")
    axB.set_xlim(lo_lim, hi_lim)
    axB.set_ylim(lo_lim, hi_lim)
    axB.set_xlabel(
        r"scale-mixture prediction $3a_\mathrm{eff}^2\,\mathrm{Var}(t_\mathrm{pla})$ [rad$^4$]"
    )
    axB.set_ylabel(r"observed $\Delta\kappa_4$ [rad$^4$]")

    header = lambda label: Line2D([], [], marker="none", ls="none",
                                  label=rf"\textbf{{{label}}}")
    blank = Line2D([], [], marker="none", ls="none", label="")
    energy_h = [
        Line2D([], [], marker="o", ls="none", ms=5.8,
               label=rf"\SI{{{E}}}{{MeV}}", **fs.estyle(E, fs.HERO))
        for E in (200, 500, 1000)
    ]
    topo_h = [
        Line2D([], [], marker=fs.MARKERS[t], ls="none", ms=5.8,
               mfc="white", mec=topo_color(t), mew=1.1,
               label=fs.TOPO_LABEL[t])
        for t in TOPOS
    ]
    # Matplotlib fills multi-column legends by column; six entries per column.
    handles = [header("Energy"), *energy_h, blank, blank,
               header("Topology"), *topo_h]
    axB.legend(
        handles=handles, ncol=2, loc="upper left", fontsize=8.5,
        columnspacing=1.2, handletextpad=0.45, labelspacing=0.28,
        borderaxespad=0.5,
    )
    fs.despine(axB)
    fs.panel_title(axB, "b", "Scale-mixture shape law")

    fs.save(fig, "fig_width_invariance")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
