#!/usr/bin/env python3
"""fig_mechanism.py -- ARCHIVAL PROTOTYPE, NOT the manuscript source.

The figure included in the paper is built from fig_mechanism.tex (TikZ); see
paper/Makefile 'figures' and figs/README.md. This Python prototype uses DIFFERENT
illustrative panel-(a) values (geometric chord lengths, not the (2.2,5.0,9.5) of the
TikZ figure) and writes to fig_mechanism_prototype.pdf so it can never overwrite the
real fig_mechanism.pdf. Kept for design history only; do not cite its numbers.

the scale-mixture mechanism (concept figure, no data).

(a) A proton pencil enters a porous slab; rays at different transverse positions
    accumulate different solid line-integrals t, drawn from p(t) with the SAME mean
    <t>=fL (accumulator bars below the slab). Each ray's projected kink is Gaussian
    of width sqrt(a t): short t -> narrow, long t -> wide.
(b) The marginal kink-angle distribution is the p(t)-weighted sum of those
    Gaussians (thin, colour-matched to the rays): a peaked, heavy-shouldered
    (leptokurtic, gamma_2>0) curve, above a same-variance Gaussian (dashed). Inset:
    the same comparison on a log axis, making the heavy tails explicit.
Schematic; illustrative numbers only. Style matches figstyle (newtx, palette).
"""
from __future__ import annotations
import numpy as np
import figstyle as fs

fs.set_style()
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Ellipse, FancyArrowPatch, Polygon

RAMP = ["#7ea4c9", "#4a76a6", "#243f5c"]      # short -> long t (light -> dark)


def main():
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(fs.FULL, fs.FULL * 0.44),
                                   width_ratios=[1.18, 0.82])

    # ===================== (a) rays through the porous slab =====================
    # Reading: one pencil enters a porous slab; three rays cross different amounts
    # of SOLID (t, thick coloured segments), and each exits into a scattering CONE
    # whose width grows with t (theta ~ N(0, a t)). The pure-void channel (grey,
    # t=0) exits straight. Below: the same t values as bars, sharing the mean fL.
    axA.set_xlim(-3.2, 17.4); axA.set_ylim(-3.5, 8.1); axA.axis("off")
    axA.set_aspect("equal")
    Lz, Hx = 10.0, 6.0
    axA.add_patch(Rectangle((0, 0), Lz, Hx, fc="#eceef1", ec=fs.HERO, lw=1.0,
                            zorder=1))
    # a few large, faint pores: a near-continuous top lane (the pure-void channel)
    # plus scattered pores lower down, so lower rays cut through more solid.
    voids = [(2.6, 5.42, 2.75, 0.50), (7.5, 5.42, 2.75, 0.50),
             (2.5, 4.05, 1.75, 0.60), (6.1, 4.20, 1.95, 0.64), (9.2, 4.05, 1.15, 0.58),
             (4.4, 2.62, 2.05, 0.58), (8.4, 2.66, 1.45, 0.54),
             (2.7, 1.05, 1.35, 0.42), (7.7, 1.02, 1.25, 0.40)]
    for cx, cy, rx, ry in voids:
        axA.add_patch(Ellipse((cx, cy), 2 * rx, 2 * ry, fc="white",
                              ec="#d4d8dd", lw=0.5, zorder=2))

    def solid_segments(y0):
        """solid x-intervals along a horizontal ray at height y0 (void complement)."""
        vs = []
        for cx, cy, rx, ry in voids:
            d = 1.0 - ((y0 - cy) / ry) ** 2
            if d > 0:
                hw = rx * np.sqrt(d)
                vs.append((max(0, cx - hw), min(Lz, cx + hw)))
        vs.sort()
        merged = []
        for a, b in vs:
            if merged and a <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], b))
            else:
                merged.append((a, b))
        solids, x = [], 0.0
        for a, b in merged:
            if a > x:
                solids.append((x, a))
            x = max(x, b)
        if x < Lz:
            solids.append((x, Lz))
        return solids, sum(b - a for a, b in solids)

    # three rays: top crosses least solid (short t), bottom the most (long t).
    ys = [4.20, 2.65, 1.05]
    ts = [solid_segments(y)[1] for y in ys]
    tmax = max(ts)
    Lf = 3.5                                    # exit-cone length

    for y0, col in zip(ys, RAMP):
        segs, t = solid_segments(y0)
        axA.plot([-2.9, 0], [y0, y0], "-", color=col, lw=1.5, zorder=4)   # incoming
        axA.plot([0, Lz], [y0, y0], "-", color=col, lw=0.7, alpha=0.30,
                 zorder=3)                                                 # faint chord
        for a, b in segs:                                                  # solid crossed
            axA.plot([a, b], [y0, y0], "-", color=col, lw=2.8, zorder=5,
                     solid_capstyle="butt")
        # exit scattering cone: half-width grows with sqrt(t) (theta ~ N(0, a t))
        half = 0.16 + 0.66 * np.sqrt(t / tmax)
        axA.add_patch(Polygon([(Lz, y0), (Lz + Lf, y0 + half), (Lz + Lf, y0 - half)],
                              closed=True, fc=col, ec="none", alpha=0.20, zorder=3))
        axA.plot([Lz, Lz + Lf], [y0, y0 + half], "-", color=col, lw=0.8, zorder=4)
        axA.plot([Lz, Lz + Lf], [y0, y0 - half], "-", color=col, lw=0.8, zorder=4)
        axA.add_patch(FancyArrowPatch((Lz, y0), (Lz + Lf, y0), arrowstyle="-|>",
                                      mutation_scale=7, color=col, lw=1.0, zorder=5))

    # the pure-void channel (t=0): threads the top lane, exits straight (no scatter)
    # -- the population that makes the theta = 0 spike (cf. panel b).
    yv, vcol = 5.5, "#8f969e"
    axA.plot([-2.9, Lz], [yv, yv], "-", color=vcol, lw=1.4, zorder=4)
    axA.add_patch(FancyArrowPatch((Lz, yv), (Lz + Lf, yv), arrowstyle="-|>",
                                  mutation_scale=7, color=vcol, lw=1.4, zorder=5))

    # labels: entry (left); the kink law as a header; void / narrow / wide at right
    axA.text(-2.9, 6.85, r"proton pencil", fontsize=8.4, color=fs.HERO,
             ha="left", va="center")
    axA.text(Lz + 0.5 * Lf, 7.6, r"exit kink $\theta\sim\mathcal{N}(0,\,a\,t)$",
             fontsize=8.6, color=fs.HERO, ha="center", va="center")
    xlab = Lz + Lf + 0.25
    axA.text(xlab, yv, r"void channel, $t=0$", fontsize=8.0, color="#5f666d",
             ha="left", va="center")
    axA.text(xlab, ys[0], r"narrow", fontsize=8.0, color=RAMP[0],
             ha="left", va="center")
    axA.text(xlab, ys[2], r"wide", fontsize=8.0, color=RAMP[2],
             ha="left", va="center")

    # accumulator bars: the solid path t of each ray, same mean fL, different spread
    bt = [0.0] + ts                             # void channel (t=0) then the rays
    bc = [vcol] + list(RAMP)
    ybar = [-0.75, -1.35, -1.95, -2.55]
    axA.text(-0.5, -1.65, "solid\npath $t$", fontsize=8.2, color=fs.HERO,
             ha="right", va="center", linespacing=1.2)
    tmean = float(np.mean(bt))
    for t, col, yb in zip(bt, bc, ybar):
        axA.plot([0, Lz], [yb, yb], "-", color="#e3e6ea", lw=3.4,
                 solid_capstyle="round", zorder=2)
        if t > 0:
            axA.plot([0, t], [yb, yb], "-", color=col, lw=3.4,
                     solid_capstyle="round", zorder=3)
        else:
            axA.plot([0.06], [yb], "o", ms=4.4, mfc=col, mec=col, zorder=3)
    axA.plot([tmean, tmean], [ybar[-1] - 0.5, ybar[0] + 0.6], "--",
             color=fs.HERO, lw=0.9, zorder=4)
    axA.text(tmean, ybar[0] + 0.78, r"$\langle t\rangle=fL$", fontsize=8.2,
             color=fs.HERO, ha="center", va="bottom")
    fs.panel(axA, "a", dy=1.0)

    # ===================== (b) the leptokurtic marginal =====================
    a = 1.0
    tvals = np.array([2.2, 5.0, 9.5])       # short / mean-ish / long (illustrative)
    w = np.array([0.42, 0.40, 0.18])        # p(t) weights
    th = np.linspace(-9, 9, 1201)

    def gauss(x, var):
        return np.exp(-x ** 2 / (2 * var)) / np.sqrt(2 * np.pi * var)

    comps = [w[i] * gauss(th, a * tvals[i]) for i in range(3)]
    marg = sum(comps)
    var_tot = float(np.sum(w * a * tvals))   # kappa_2 = a<t>
    ref = gauss(th, var_tot)

    for i in range(3):
        axB.plot(th, comps[i], "-", color=RAMP[i], lw=1.0, zorder=3)
    axB.plot(th, ref, "--", color=fs.MUTE, lw=1.3, zorder=4)
    axB.plot(th, marg, "-", color=fs.HERO, lw=2.0, zorder=5)
    axB.set_xlim(-9, 9); axB.set_ylim(0, marg.max() * 1.18)
    axB.set_yticks([])
    axB.set_xlabel(r"projected kink angle $\theta$")
    axB.set_ylabel(r"density")
    axB.annotate(r"scale mixture", xy=(0.35, marg.max() * 0.97),
                 xytext=(-8.4, marg.max() * 1.05), fontsize=8.4, color=fs.HERO,
                 ha="left", va="center",
                 arrowprops=dict(arrowstyle="-", color=fs.HERO, lw=0.7,
                                 shrinkA=2, shrinkB=2))
    # parked in the empty mid-left band: both curves are below y=0.13 for
    # theta < -4.2, so neither line can be sliced by a curve here.
    axB.text(-8.5, marg.max() * 0.79, r"same-variance", fontsize=7.8,
             color=fs.MUTE, ha="left", va="bottom")
    axB.text(-8.5, marg.max() * 0.705, r"Gaussian (dashed)", fontsize=7.8,
             color=fs.MUTE, ha="left", va="bottom")
    fs.despine(axB)
    fs.panel(axB, "b", dy=1.0)

    # inset: heavy tails on a log axis
    axi = axB.inset_axes([0.63, 0.50, 0.35, 0.45])
    axi.semilogy(th, marg, "-", color=fs.HERO, lw=1.4, zorder=4)
    axi.semilogy(th, ref, "--", color=fs.MUTE, lw=1.2, zorder=3)
    axi.set_xlim(0, 8.5); axi.set_ylim(1e-4, marg.max() * 1.4)
    axi.set_xticks([0, 4, 8]); axi.tick_params(labelsize=7.0, length=2.2)
    axi.set_yticks([1e-3, 1e-1])
    # the tick labels land in the MAIN panel's data region: give them an opaque
    # backing so the mixture curve cannot slice through "10^-3" / the "0".
    for lab in axi.get_yticklabels() + axi.get_xticklabels():
        lab.set_bbox(dict(facecolor="white", edgecolor="none", pad=0.7))
    axi.set_title(r"heavy tails", fontsize=7.8, pad=2)
    fs.despine(axi)

    fs.save(fig, "fig_mechanism_prototype")   # never overwrite the real TikZ figure
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
