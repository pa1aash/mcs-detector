#!/usr/bin/env python3
"""fig_affordability.py -- the cost and the reach of the shape channel.

(a) Statistical cost: fractional 95% CI on the deconvolved Delta_kappa4 vs proton
    count, scaled from the simulated campaign anchors (one uniform filled circle per
    topology) by the 1/sqrt(N) moment-estimator law; beam-time top axis at 1 MHz;
    20/10% targets. Curves are direct-labelled; a single frameless key defines
    the anchor marker so all three filled circles read as one class.
(b) Metrology reach: relative 1-sigma bootstrap width of the recovered (f, N_eff)
    vs proton count for the anchored rectilinear f=0.30 200 MeV run. Filled markers
    are GENUINE bootstrap resamples at that size (D9: no law-scaled points); the
    open marker at 1e7 is an explicit sqrt-N extrapolation (dashed segment).
Reads e3_affordability.json + results/N3/posterior.json.
"""
from __future__ import annotations
import numpy as np
from matplotlib.lines import Line2D
import figstyle as fs

fs.set_style()
import matplotlib.pyplot as plt

E3 = fs.load("e3_affordability.json")
PN = fs.loadr("N3", "posterior.json")
flux = 1.0e6                                     # protons/s (johnson2017)

fig, (ax, axb) = plt.subplots(1, 2, figsize=(fs.FULL, fs.FULL * 0.49))

# ================= (a) statistical cost =================
NG_LO, NG_HI = 1.0e5, 5.0e8
Ng = np.logspace(np.log10(NG_LO), np.log10(NG_HI), 200)

# fixed draw order: rectilinear -> gyroid -> Voronoi (bottom to top of the fan)
order = ["rectilinear", "gyroid", "voronoi"]
cfg = {c["topology"]: c for c in E3["configs"] if c["E"] == 200}
# where to park each direct label relative to its line's right end (log-space nudge)
lab_pos = {"rectilinear": (1.22, "bottom"), "gyroid": (1.28, "bottom"),
           "voronoi": (1.20, "bottom")}

for topo in order:
    c = cfg[topo]
    col = fs.COLORS[topo]
    scale = np.sqrt(c["n_protons_measured"] / Ng)
    line = c["cif_dk4"] * 100 * scale
    scale_ls = {"rectilinear": (0, (6, 2.2)),
                "gyroid": (0, (4, 1.8, 1.2, 1.8)),
                "voronoi": (0, (1.4, 1.5))}[topo]
    ax.plot(Ng, line, ls=scale_ls, lw=1.7, color=col, zorder=3,
            solid_capstyle="round")
    # simulated campaign anchor: one uniform filled circle for every topology, so
    # the shape carries no meaning (topology is already the line colour + label).
    # White halo + high zorder let the circle punch cleanly through any gridline.
    yend = c["cif_dk4"] * 100
    ax.plot(c["n_protons_measured"], yend, "o", ms=7,
            color=fs.HERO, mec="white", mew=1.1, zorder=6)
    # direct curve label at the right terminus
    fac, va = lab_pos[topo]
    ax.text(NG_HI * 1.08, (c["cif_dk4"] * 100 * np.sqrt(c["n_protons_measured"] / NG_HI)) * fac,
            fs.TOPO_LABEL[topo], color=col, fontsize=9.5, ha="left", va=va,
            fontweight="bold")

# frameless key: filled dot = the simulated campaign anchor; the colour lines are the
# 1/sqrt(N) moment-estimator scaling from that anchor (anchored vs extrapolated made
# explicit, as in panel b).
anchor_key = Line2D([], [], marker="o", ms=7, color=fs.HERO, mec="white",
                    mew=1.1, ls="none")
scale_key = Line2D([], [], color=fs.MUTE, lw=fs.LW_DATA, ls=(0, (4, 2)))
ax.legend([anchor_key, scale_key],
          ["Geant4 campaign", r"$N^{-1/2}$ extrapolation"], loc="lower left",
          frameon=False, fontsize=9, handletextpad=0.5, labelspacing=0.35,
          borderaxespad=0.9)

# Two decision-relevant target precisions; additional values are given in the text.
for p in (20, 10):
    ax.axhline(p, color=fs.FAINT, lw=0.8, ls=(0, (1, 2)), zorder=1)
    ax.text(1.18e5, p * 0.94, rf"{p}\%", fontsize=8.5, color="#7a7a7a",
            ha="left", va="top")

ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(8.5e4, 8.5e8)
ax.set_ylim(0.55, 550)
ax.set_xlabel(r"protons $N$")
ax.set_ylabel(r"95\% CI full width$/|\Delta\kappa_4|$ [\%]")

sx = ax.secondary_xaxis("top", functions=(lambda n: n / flux, lambda t: t * flux))
sx.set_xlabel(r"beam time at \SI{1}{MHz} [s]")
sx.tick_params(labelsize=9)
fs.despine(ax, keep=("left", "bottom"))
ax.spines["top"].set_visible(True)
fs.panel_title(ax, "a", "Statistical cost", y=1.31)

# ================= (b) metrology reach =================
SER = (("N_eff", fs.HERO, "o", r"$N_\mathrm{eff}$"),
       ("f", fs.MUTE, "s", r"$f$"))
post = PN["posteriors"]
Ns = np.array([r["N"] for r in post], float)
meas = np.array([bool(r["measured"]) for r in post])

for key, col, mk, lab in SER:
    rel = np.array([100.0 * r[key][1] / r[key][0] for r in post])
    # measured span
    axb.plot(Ns[meas], rel[meas], "-", color=col, lw=1.6, zorder=3)
    axb.plot(Ns[meas], rel[meas], mk, ms=7, mfc=col, mec="white", mew=0.8,
             ls="", zorder=5, label=lab)
    # sqrt(N) extrapolation to 1e7, drawn open + dashed
    if (~meas).any():
        xe = np.r_[Ns[meas][-1], Ns[~meas]]
        ye = np.r_[rel[meas][-1], rel[~meas]]
        axb.plot(xe, ye, "--", color=col, lw=1.2, zorder=2)
        axb.plot(Ns[~meas], rel[~meas], mk, ms=7, mfc="white", mec=col,
                 mew=1.3, ls="", zorder=5)

# Reference N^{-1/2} slope, normalised to the first genuine N_eff bootstrap.
neff_rel = np.array([100.0 * r["N_eff"][1] / r["N_eff"][0] for r in post])
guide_x = np.array([Ns[meas][0], Ns[meas][-1]])
guide_y = neff_rel[meas][0] * np.sqrt(guide_x[0] / guide_x)
axb.plot(guide_x, guide_y, color=fs.FAINT, lw=0.9, ls=(0, (2, 2)), zorder=1)

axb.text(0.97, 0.96, r"open: $\sqrt{N}$ extrapolation", transform=axb.transAxes,
         fontsize=8.5, color="#6f6f6f", ha="right", va="top")
axb.set_xscale("log"); axb.set_yscale("log")
axb.set_xlim(7.0e4, 1.5e7)
axb.set_ylim(0.06, 25)
axb.set_xlabel(r"protons $N$")
axb.set_ylabel(r"relative $1\sigma$ width [\%]")
axb.legend(loc="lower left", fontsize=9.5, handletextpad=0.4)
fs.despine(axb)
fs.panel_title(axb, "b", r"Recovery of $f$ and $N_\mathrm{eff}$", y=1.31)

fs.save(fig, "fig_affordability")
