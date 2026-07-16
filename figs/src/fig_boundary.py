#!/usr/bin/env python3
"""fig_boundary.py -- practical and extrapolated validity boundaries.

(a) gamma2 profile (evidence): x = cell size, y = excess kurtosis gamma2.
    Transport curves (rectilinear, gyroid; 200 MeV) over the computed sweep
    5 um - 3 mm, with full-Geant4 anchors (open markers, WentzelVI) at the
    2.5 mm printable cell (both topologies) and the 0.2 mm foam-scale cell
    (rectilinear).
(b) boundary map (message): x = proton kinetic energy, y = cell.
    c_break lines (both topologies), 100 MeV Geant4-confirmed onsets (open
    markers), and an explicitly hatched order-of-magnitude region for the
    threshold-defined, extrapolated rectilinear shape-recovery scale.  The 5 um
    sweep floor is drawn so the extrapolation cannot be mistaken for a CI.

Reads homog_boundary.json, e5_msc_systematic.json (foam anchor), e0_break.json.
Geant4 gamma2 at the 2.5 mm cell (200 MeV): rectilinear 2.22, gyroid 1.28
  % SOURCE: data/analysis/homog_boundary.md, C1 table.
"""
from __future__ import annotations
import numpy as np
import figstyle as fs

fs.set_style()
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple

H = fs.load("homog_boundary.json")
E5 = fs.load("e5_msc_systematic.json")
B = fs.load("e0_break.json")

fig, (axA, axB) = plt.subplots(1, 2, figsize=(fs.FULL, fs.FULL * 0.54),
                               gridspec_kw={"width_ratios": [1.05, 1.25]})
FOAM_UM = (fs.FOAM_BAND[0] * 1e3, fs.FOAM_BAND[1] * 1e3)   # 100-500 um

# ---------- (a) gamma2 profile: x = cell, y = gamma2 ----------
axA.axvspan(*FOAM_UM, color=fs.FOAM, alpha=0.32, lw=0, zorder=0)
for t in ("rectilinear", "gyroid"):
    cur = H["gamma2_curves"][t]["200"]
    cells = np.array(sorted(float(k) for k in cur))
    axA.plot(cells * 1e3, [cur[f"{x:g}"] for x in cells],
             ls=fs.LINESTYLES[t], color=fs.COLORS[t], lw=1.8, zorder=3)

# Geant4 anchors loaded from the seeded campaign, all 200 MeV.
J = fs.load("e2_results_combined.json")
g4_25 = {r["topology"]: r["dk4"] / r["k2"] ** 2 for r in J["rows"]
         if r["E"] == 200 and abs(r["f_designed"] - 0.40) < 1e-9
         and r["topology"] in ("rectilinear", "gyroid")}
# white halo so open markers punch cleanly out of the foam/print band.
def _halo(ax, x, y, marker):
    ax.scatter(x, y, s=95, marker=marker, facecolor="white",
               edgecolor="none", zorder=4)
for t, g in g4_25.items():
    _halo(axA, [2500], [g], fs.MARKERS[t])
    fs.phys_scatter(axA, [2500], [g], fs.COLORS[t], "WentzelVI",
                    marker=fs.MARKERS[t], size=38, zorder=5, filled=False)
g_w = E5["foam"]["locked_gamma2_g4"]     # 1.965, 0.2 mm rectilinear
_halo(axA, [200], [g_w], fs.MARKERS["rectilinear"])
fs.phys_scatter(axA, [200], [g_w], fs.COLORS["rectilinear"], "WentzelVI",
                marker=fs.MARKERS["rectilinear"], size=38, zorder=5, filled=False)

# direct topology labels and practical band
axA.text(1800, 2.27, "rectilinear", fontsize=9.5, color=fs.COLORS["rectilinear"],
         ha="right", va="bottom")
axA.text(700, 1.12, "gyroid", fontsize=9.5, color=fs.COLORS["gyroid"],
         ha="left", va="top")
axA.text(np.sqrt(FOAM_UM[0] * FOAM_UM[1]), 0.13, "foam / print",
         fontsize=8.5, color=fs.FOAM_EDGE, ha="center", va="bottom")

# one method legend (the figure's only legend): line = transport, open = Geant4.
# both Geant4 marker shapes are shown (circle = rectilinear, triangle = gyroid,
# in their topology colours) so the open triangle is defined, not just the circle.
g4_rect = Line2D([0], [0], color="none", marker=fs.MARKERS["rectilinear"],
                 mfc="white", mec=fs.COLORS["rectilinear"], mew=1.1, ms=6.5)
g4_gyr = Line2D([0], [0], color="none", marker=fs.MARKERS["gyroid"],
                mfc="white", mec=fs.COLORS["gyroid"], mew=1.1, ms=6.5)
proxies = [Line2D([0], [0], color="#4d4d4d", lw=1.7), (g4_rect, g4_gyr)]
leg = axA.legend(proxies, ["transport", "Geant4"], title=r"\SI{200}{MeV}",
                 loc="upper left", bbox_to_anchor=(0.015, 0.985),
                 fontsize=9, title_fontsize=9, handlelength=1.5,
                 labelspacing=0.5, handletextpad=0.6, borderaxespad=0.0,
                 handler_map={tuple: HandlerTuple(ndivide=None, pad=0.7)})
leg.get_title().set_ha("left")

axA.axhline(0.0, color=fs.MUTE, ls=(0, (4, 2)), lw=0.9, zorder=1)
axA.set_xscale("log")
axA.set_xlim(3.5, 3600)
axA.set_ylim(-0.05, 2.55)
axA.set_xticks([5, 10, 100, 1000])
axA.set_xticklabels(["5", "10", "100", "1000"])
axA.set_xlabel(r"cell size $c$ [\si{\micro\metre}]")
axA.set_ylabel(r"excess kurtosis $\gamma_2$")
fs.despine(axA); fs.panel_title(axA, "a", "Shape persists across practical cells")

# ---------- (b) boundary map: x = energy, y = cell ----------
axB.set_yscale("log")
# Crop the empty extrapolation tail: the shape-recovery region only matters up to
# its ~1 um upper edge, so showing it as a labelled strip (rather than 3 empty
# decades) gives the c_break-vs-foam-band message most of the panel.
axB.set_ylim(7e-2, 1.3e3)
axB.set_ylabel(r"cell size $c$ [\si{\micro\metre}]")
axB.axhspan(*FOAM_UM, color=fs.FOAM, alpha=0.32, lw=0, zorder=0)
Pm = np.array(B["momenta_MeV"], float)
cb = {t: np.array([B["cell_break_mm"][t][str(int(p))] for p in Pm]) * 1e3
      for t in ("rectilinear", "gyroid")}
ch_lo = np.array([H["cell_homog"]["rectilinear"][str(int(p))]["0.1"] for p in Pm]) * 1e3
ch_hi = np.array([H["cell_homog"]["rectilinear"][str(int(p))]["1.0"] for p in Pm]) * 1e3

# foam / print pore band label (placed where c_break has dropped below it)
axB.text(640, 224, "foam / print pores", fontsize=9, color=fs.FOAM_EDGE,
         va="center", ha="center", zorder=7)

# c_break lines + 100 MeV Geant4 onsets (open markers)
for t in ("rectilinear", "gyroid"):
    axB.plot(Pm, cb[t], ls=fs.LINESTYLES[t], color=fs.COLORS[t], lw=1.8, zorder=4)
    fs.phys_scatter(axB, [100], [cb[t][0]], fs.COLORS[t], "WentzelVI",
                    marker=fs.MARKERS[t], size=34, zorder=6, filled=False)
axB.text(1040, 57, r"$c_\mathrm{break}$", fontsize=10, color="#3a3a3a",
         ha="left", va="center")
# One organised legend for panel (b): the two c_break curves (by topology) and
# the shared Geant4 100 MeV onset markers, in the clear zone above the floor.
g4_rb = Line2D([0], [0], color="none", marker=fs.MARKERS["rectilinear"],
               mfc="white", mec=fs.COLORS["rectilinear"], mew=1.1, ms=6)
g4_gb = Line2D([0], [0], color="none", marker=fs.MARKERS["gyroid"],
               mfc="white", mec=fs.COLORS["gyroid"], mew=1.1, ms=6)
cb_handles = [
    Line2D([0], [0], color=fs.COLORS["rectilinear"],
           ls=fs.LINESTYLES["rectilinear"], lw=1.8),
    Line2D([0], [0], color=fs.COLORS["gyroid"],
           ls=fs.LINESTYLES["gyroid"], lw=1.8),
    (g4_rb, g4_gb),
]
cb_labels = ["rectilinear", "gyroid", r"Geant4, \SI{100}{MeV}"]
axB.legend(cb_handles, cb_labels, loc="upper left", bbox_to_anchor=(0.02, 0.74),
           fontsize=8.5, frameon=True, framealpha=0.92, edgecolor="#cfcfcf",
           handlelength=2.1, labelspacing=0.4, handletextpad=0.6,
           borderaxespad=0.0,
           handler_map={tuple: HandlerTuple(ndivide=None, pad=0.5)})

# Shape-recovery estimate: the full envelope of the rectilinear gamma2=1 to 0.1
# threshold crossings.  This is a definition range and extrapolation, not a CI.
rec_lo = max(float(np.min(ch_lo)), 1.3e-3)
rec_hi = float(np.max(ch_hi))
axB.axhspan(rec_lo, rec_hi, facecolor="#eeeeee", edgecolor=fs.MUTE,
            hatch="////", linewidth=0.6, zorder=1)
axB.axhline(5.0, color=fs.MUTE, lw=0.8, ls=(0, (3, 2)), zorder=2)
axB.text(104, 6.4, r"\SI{5}{\micro\metre} sweep floor", fontsize=8.5,
         color=fs.MUTE, ha="left", va="bottom")
axB.text(105, 0.12, "shape recovery (extrapolated)", fontsize=8.5,
         color=fs.MUTE, ha="left", va="bottom",
         bbox=dict(facecolor="#ededed", edgecolor="none", pad=2.0, alpha=0.92))

# the payoff: the >= 2-decade gap between the pore band and c_homog (dark = hero)
axB.annotate("", xy=(340, rec_hi), xytext=(340, FOAM_UM[0]),
             arrowprops=dict(arrowstyle="<->", color=fs.HERO, lw=1.1,
                             shrinkA=0, shrinkB=0))
axB.text(372, 9.5, r"$\geq 10^{2}\times$", fontsize=10.5, color=fs.HERO,
         ha="left", va="center")

axB.set_xscale("log")
axB.set_xlim(92, 1300)
axB.set_xticks([100, 200, 500, 1000])
axB.set_xticklabels(["100", "200", "500", "1000"])
from matplotlib.ticker import NullLocator  # noqa: E402
axB.xaxis.set_minor_locator(NullLocator())   # drop x minors; keep y log minors
axB.set_xlabel(r"proton kinetic energy $E$ [\si{MeV}]")
fs.despine(axB); fs.panel_title(axB, "b", "Practical validity map")

fs.save(fig, "fig_boundary")
