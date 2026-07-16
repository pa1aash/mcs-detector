#!/usr/bin/env python3
"""Engine and transport-tool validation.

Panel (a) shows the Geant4 Highland core-width residual across the control-slab
thickness scan.  The offset is step-converged but has a visible thickness trend,
so the figure reports the full -5.2% to -3.1% range rather than calling it flat.

Panel (b) expresses the transport comparison as the directly readable closure
ratio Geant4/transport.  Open markers are the six geometry-induced Delta-kappa4
checks at 0.5, 1.0 and 2.5 mm; the filled point is the independent, floor-free
path-variance closure at the 0.2 mm foam-scale rectilinear cell.  All error bars
are Geant4 bootstrap 95% confidence intervals.
"""
from __future__ import annotations

import numpy as np
import figstyle as fs

fs.set_style()
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullLocator

TG = fs.load("transport_vs_g4.json")

# Foam-scale path-variance check, committed in transport_vs_g4.md.
RATIO_TA = 0.868
RATIO_G4 = 0.854
RATIO_G4_ERR = 0.068 / 2.301

# Highland core-width residual (sigma_core/theta_0 - 1) [%].
# SOURCE: VALIDATION.md section 1 (committed, locked).
HRES = {
    200:  ([3, 4, 8, 16], [-5.21, -5.14, -4.55, -3.07]),
    500:  ([3, 4, 8, 16], [-5.06, -4.87, -4.55, -3.65]),
    1000: ([3, 4, 8, 16], [-4.49, -4.28, -3.94, -3.31]),
}


def main():
    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=(fs.FULL, fs.FULL * 0.47), width_ratios=[1.08, 0.92]
    )

    # ------------------------- (a) Highland trend -------------------------
    axA.axhline(0, color=fs.REF, lw=fs.LW_REF, zorder=1)
    energy_ls = {200: "-", 500: (0, (4.5, 2.0)), 1000: (0, (1.4, 1.5))}
    handles = []
    for E in (200, 500, 1000):
        t, residual = HRES[E]
        t = np.asarray(t, float)
        residual = np.asarray(residual, float)
        axA.plot(t, residual, color=fs.HERO, lw=1.3, ls=energy_ls[E], zorder=2)
        axA.plot(
            t, residual, fs.MARKERS["rectilinear"], ms=fs.MS, mew=1.1,
            zorder=3, **fs.estyle(E, fs.HERO)
        )
        handles.append(Line2D(
            [0], [0], color=fs.HERO, lw=1.3, ls=energy_ls[E],
            marker=fs.MARKERS["rectilinear"], ms=fs.MS, mew=1.1,
            **fs.estyle(E, fs.HERO)))
    axA.legend(handles, [r"\SI{200}{MeV}", r"\SI{500}{MeV}", r"\SI{1000}{MeV}"],
               loc="upper right", bbox_to_anchor=(0.99, 0.86), fontsize=8.5,
               frameon=False, handlelength=2.6, labelspacing=0.35,
               borderaxespad=0.0)
    axA.set_xscale("log")
    axA.set_xlim(2.5, 19)
    axA.set_ylim(-5.65, 0.35)
    axA.set_xticks([3, 4, 8, 16])
    axA.set_xticklabels(["3", "4", "8", "16"])
    axA.xaxis.set_minor_locator(NullLocator())
    axA.set_xlabel(r"control-slab thickness $t$ [mm]")
    axA.set_ylabel(r"core-width residual $\sigma_\mathrm{core}/\theta_0-1$ [\%]")
    fs.despine(axA)
    fs.panel_title(axA, "a", "Highland width trend")

    # ----------------------- (b) transport closure ------------------------
    axB.axhline(1.0, color=fs.HERO, lw=1.2, zorder=1)
    for row in TG:
        topo = row["name"]
        color = fs.COLORS[topo]
        pred = row["dk4_ta"]
        ratio = row["dk4_g4"] / pred
        lo = row["dk4_g4_lo"] / pred
        hi = row["dk4_g4_hi"] / pred
        axB.errorbar(
            row["cell"], ratio,
            yerr=[[ratio - lo], [hi - ratio]],
            fmt=fs.MARKERS[topo], ms=6.1, mew=1.1, zorder=4,
            **fs.estyle(200, color), **fs.ebar(color),
        )

    # The former inset becomes a normal-size seventh validation point.
    foam_ratio = RATIO_G4 / RATIO_TA
    foam_err = RATIO_G4_ERR / RATIO_TA
    axB.errorbar(
        0.2, foam_ratio, yerr=foam_err,
        fmt=fs.MARKERS["rectilinear"], ms=6.5, mfc=fs.COLORS["rectilinear"],
        mec=fs.COLORS["rectilinear"], mew=1.0, zorder=5,
        **fs.ebar(fs.COLORS["rectilinear"]),
    )

    axB.set_xscale("log")
    axB.set_xlim(0.16, 3.2)
    axB.set_ylim(0.45, 1.28)
    axB.set_xticks([0.2, 0.5, 1.0, 2.5])
    axB.set_xticklabels(["0.2", "0.5", "1", "2.5"])
    axB.xaxis.set_minor_locator(NullLocator())
    axB.set_yticks([0.5, 0.75, 1.0, 1.25])
    axB.set_xlabel(r"validation cell size [mm]")
    axB.set_ylabel(r"Geant4 / transport")

    rect_h = Line2D([], [], marker=fs.MARKERS["rectilinear"], ls="none", ms=6,
                    mfc="white", mec=fs.COLORS["rectilinear"], mew=1.1,
                    label=r"rectilinear")
    gyr_h = Line2D([], [], marker=fs.MARKERS["gyroid"], ls="none", ms=6,
                   mfc="white", mec=fs.COLORS["gyroid"], mew=1.1,
                   label=r"gyroid")
    foam_h = Line2D([], [], marker=fs.MARKERS["rectilinear"], ls="none", ms=6,
                    mfc=fs.COLORS["rectilinear"], mec=fs.COLORS["rectilinear"],
                    label=r"foam")
    # placed in the empty bottom-left corner (all data sits at cell >= 0.5 mm).
    axB.legend(handles=[rect_h, gyr_h, foam_h], loc="lower left", fontsize=8.5,
               frameon=False, handletextpad=0.4, labelspacing=0.3,
               borderaxespad=0.4)
    fs.despine(axB)
    fs.panel_title(axB, "b", "Transport closure")

    fs.save(fig, "fig_validation")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
