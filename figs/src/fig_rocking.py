#!/usr/bin/env python3
"""fig_rocking.py -- orientation dependence.

(a) N_eff(alpha) for rectilinear and gyroid: curves = tilted-chord transport
    prediction; a bottom rug flags the angles at which Geant4 was run. The
    gyroid rocking curve peaks at the chord-lattice commensurability tan(alpha)=1/4.
    The window is the scan range the caption describes (rise + peak); the analytic
    prediction is rendered on a fine angular grid so it reads smooth.
(b) the collapse prefactor C*N_eff(alpha) from the Geant4 tilt scan with bootstrap
    95% error bars (results/M2/collapse_theta.json); dashed line at 1.

Reads results/M2/collapse_theta.json; recomputes the smooth N_eff(alpha) curve via
analysis/m2_tilt.neff_tilted. Presentation only -- data and physics are unchanged.
"""
from __future__ import annotations
import json, os, sys
import numpy as np
import figstyle as fs

fs.set_style()
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(fs.ROOT, "analysis"))
import m2_tilt  # noqa

AMAX = 18.0          # display window: rise + commensurability peak (14.0 deg)
APEAK = 14.04        # tan(alpha) = 1/4
ASHADE = 10.0        # Geant4 was run up to 10 deg; beyond is prediction only


def main():
    cj = fs.loadr("M2", "collapse_theta.json")
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(fs.FULL, fs.FULL * 0.49))

    # ---- (a) N_eff(alpha): smooth prediction curves + Geant4-run ticks ----------
    grid = np.arange(0.0, AMAX + 0.001, 0.25)   # fine grid -> smooth render
    for topo in ("gyroid", "rectilinear"):
        ne = np.array([m2_tilt.neff_tilted(topo, float(t))["N_eff"] for t in grid])
        axL.plot(grid, ne, ls=fs.LINESTYLES[topo], color=fs.COLORS[topo],
                 lw=1.8, zorder=4,
                 solid_capstyle="round")
        runs = sorted(r["theta"] for r in cj["rows"] if r["topo"] == topo)
        rug_y = 0.035 if topo == "rectilinear" else 0.075
        axL.plot(runs, [rug_y] * len(runs), marker=fs.MARKERS[topo],
                 transform=axL.get_xaxis_transform(), color=fs.COLORS[topo],
                 ms=4.8, mfc="white", mew=1.0, ls="", zorder=5, clip_on=False)

    axL.set_yscale("log")
    axL.set_xlim(-0.5, AMAX + 0.3)
    axL.set_ylim(1.7, 900)
    axL.set_xticks([0, 5, 10, 15])
    axL.set_xlabel(r"incidence tilt $\alpha$ [deg]")
    axL.set_ylabel(r"$N_\mathrm{eff}(\alpha)$")

    # one subtle band: beyond 10 deg is transport prediction only (caption says so)
    axL.axvspan(ASHADE, AMAX + 0.3, color=fs.FAINT, alpha=0.16, lw=0, zorder=0)
    # a very short in-band cue so the band reads without the caption
    axL.text(17.7, 660, "prediction only", fontsize=8.5, color=fs.MUTE,
             ha="right", va="top", zorder=1)

    # label sits directly above the (self-evident) commensurability peak; no
    # arrowhead/caret -- the peak is the gyroid maximum right below it.
    axL.text(APEAK, 250, r"$\tan\alpha=1/4$", fontsize=9, color=fs.MUTE,
             ha="center", va="bottom")
    # direct curve labels (no legend in this figure)
    axL.text(3.4, 8.4, "gyroid", fontsize=9.5, color=fs.COLORS["gyroid"],
             ha="left", va="bottom")
    axL.text(12.5, 3.55, "rectilinear", fontsize=9.5,
             color=fs.COLORS["rectilinear"], ha="center", va="top")
    fs.despine(axL)
    fs.panel_title(axL, "a", "Predicted orientation response")

    # ---- (b) prefactor C*N_eff(alpha) with bootstrap 95% errors -----------------
    YLO, YHI = 0.1, 1.72
    lab_pos = {"gyroid": (2.4, 1.32, "bottom"),
               "rectilinear": (7.2, 0.82, "top")}
    for topo in ("rectilinear", "gyroid"):
        runs = sorted([r for r in cj["rows"] if r["topo"] == topo and r["dk4"] > 0],
                      key=lambda r: r["theta"])
        col = fs.COLORS[topo]
        th = np.array([r["theta"] for r in runs], float) + \
            (-0.12 if topo == "rectilinear" else 0.12)
        v = np.array([r["ratio"] for r in runs])
        lo = np.array([r.get("ratio_lo", r["ratio"]) for r in runs])
        hi = np.array([r.get("ratio_hi", r["ratio"]) for r in runs])
        axR.errorbar(th, v, yerr=[v - lo, hi - v], fmt=fs.MARKERS[topo],
                     color=col, ms=5.8, mfc="white", mew=1.1, ls="",
                     elinewidth=0.9, capsize=2.4, capthick=0.8, zorder=4)
        x, y, va = lab_pos[topo]
        axR.text(x, y, topo, fontsize=9.5, color=col, ha="center", va=va)

    axR.axhline(1.0, color=fs.MUTE, lw=0.9, ls=(0, (4, 3)), zorder=2)
    axR.set_xlim(-0.6, 11.2)
    axR.set_ylim(YLO, YHI)
    axR.set_xticks([0, 2, 4, 6, 8, 10])
    axR.set_yticks([0.5, 1.0, 1.5])
    axR.set_xlabel(r"incidence tilt $\alpha$ [deg]")
    axR.set_ylabel(r"$C\,N_\mathrm{eff}(\alpha)$")
    fs.despine(axR)
    fs.panel_title(axR, "b", "Geant4 closure under tilt")

    fs.save(fig, "fig_rocking")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
