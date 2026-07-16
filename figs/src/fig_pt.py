#!/usr/bin/env python3
"""fig_pt.py -- line-integral probability at a common mean.

Three horizontal small multiples (rectilinear, gyroid, Voronoi at f = 0.40)
use one shared probability scale.  The continuum is probability per
fixed 0.2 mm bin; endpoint atoms P(t=0) and P(t=L) use the SAME vertical units,
so an atom and a continuum bin can be compared without an arbitrary mapping.
The shared dashed line marks the designed mean fL: nearly the same mean, very
different spread and hence N_eff.
Reads data/analysis/pt_hist.json (tools/make_pt_hist.py; raw campaign voxels).
"""
from __future__ import annotations
import numpy as np
import figstyle as fs

fs.set_style()
import matplotlib.pyplot as plt

D = fs.load("pt_hist.json")
bins = np.array(D["bins_mm"])
dw = bins[1] - bins[0]                     # native bin width (0.1 mm)
fL = D["f_designed"] * D["L_mm"]
L = D["L_mm"]

K = 2                                      # display rebin: 0.1 mm -> 0.2 mm
bins_d = bins[::K]
ctr = 0.5 * (bins_d[1:] + bins_d[:-1])
BINW = K * dw

ORDER = ("rectilinear", "gyroid", "voronoi")
fig, axes = plt.subplots(1, 3, sharex=True, sharey=True,
                         figsize=(fs.FULL, fs.FULL * 0.34))

for ax, t in zip(axes, ORDER):
    d = D["topologies"][t]
    col = fs.COLORS[t]
    # Native density integrates to one.  Convert it to probability per native bin,
    # subtract the exact endpoint masses, then add adjacent bins to 0.2 mm.
    p_native = np.array(d["density"], float) * dw
    p_void = float(d.get("p_void", 0.0))
    p_solid = float(d.get("p_solid", 0.0))
    p_native[0] = max(0.0, p_native[0] - p_void)
    p_native[-1] = max(0.0, p_native[-1] - p_solid)
    prob = p_native.reshape(-1, K).sum(axis=1)

    ax.fill_between(ctr, 0, prob, step="mid", color=col, alpha=0.18,
                    lw=0, zorder=2)
    ax.step(ctr, prob, where="mid", color=col, lw=1.6, zorder=3)
    ax.axvline(fL, color=fs.MUTE, lw=0.9, ls=(0, (5, 4)), zorder=1)

    # Point masses share the bin-probability scale; no density-to-height mapping.
    def lollipop(x, p, pos, ha, xlab):
        ax.plot([x, x], [0, p], color=col, lw=2.2,
                solid_capstyle="butt", zorder=4)
        ax.plot([x], [p], marker="o", ms=5.5,
                mfc=col, mec=col, zorder=5)
        ax.text(xlab, p + 0.012, rf"$P(t{{=}}{pos})={p:.2f}$",
                fontsize=8.5, color=col, ha=ha, va="bottom")

    if p_void > 0.005:
        lollipop(0.0, p_void, "0", "left", 0.35)
    if p_solid > 0.005:
        lollipop(L, p_solid, "L", "right", L - 0.35)

    neff = D["f_designed"] * (1.0 - D["f_designed"]) * L**2 / d["var_t_mm2"]
    fs.panel_title(ax, "abc"[ORDER.index(t)],
                   rf"{fs.TOPO_LABEL[t]}\quad $N_\mathrm{{eff}}={neff:.1f}$")
    ax.set_xlim(-0.45, L + 0.45)
    ax.set_xticks([0, 2, 4, 6, 8, 10])
    fs.despine(ax)

axes[0].set_ylim(0, 0.295)
axes[0].set_yticks([0, 0.1, 0.2])
axes[0].set_ylabel(rf"probability per \SI{{{BINW:g}}}{{mm}} bin")
axes[2].text(fL, 0.278, r"designed $fL$", fontsize=8.5,
             color=fs.MUTE, ha="center", va="top")
fig.supxlabel(r"material line-integral $t$ [mm]", fontsize=10)

fs.save(fig, "fig_pt")
