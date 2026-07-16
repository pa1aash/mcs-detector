#!/usr/bin/env python3
"""fig_impact.py -- M1 reconstruction impact.

(a) rectilinear-struct against matched-slab reconstructed kink distributions,
    normalised to their fitted sigma, 0.1-sigma step bins, log-y, unit Gaussian for
    reference. The narrow central component of the structured curve is the physical
    void-channel population (t = 0 for 34% of protons through the rectilinear f=0.30
    lattice; measured from the run's per-proton path lengths), annotated as such.
(b) tail fraction P(|theta| > k sigma) against integer k, struct against slab, with
    the structured/slab ratio annotated above each structured point.

Reads the tracked results/M1/f6_hist.json (derived from the raw event arrays by
tools/make_impact_hist.py; the raw f6_data.npz is regenerable via
analysis/tracking/run_m1_analysis.py).

Presentation notes: ONE legend for the whole figure (panel a); panel (b) reuses the
same fixed red/grey encoding and carries only the required per-k ratio labels. All
on-canvas text >= 9 pt. Hero = rectilinear (red, heavier); slab recedes (grey).
"""
from __future__ import annotations
import json, os
import numpy as np
import figstyle as fs

R = os.path.join(fs.ROOT, "results", "M1", "f6_hist.json")
# p_void = P(t==0) is read from the JSON (computed by tools/make_impact_hist.py from
# the pencil run's tpla branch), not hard-coded here.


def main():
    if not os.path.exists(R):
        print("F6 histogram missing (run tools/make_impact_hist.py first):", R)
        return 1
    d = json.load(open(R))
    bins = np.array(d["bins"])
    ctr = 0.5 * (bins[1:] + bins[:-1])
    bw = bins[1] - bins[0]
    fs.set_style()
    fig, (axL, axR) = fs.plt.subplots(1, 2, figsize=(fs.FULL, fs.FULL * 0.49))

    # ---- (a) normalised kink distributions, 0.1-sigma bins, log-y ----------------
    # visual hierarchy: rectilinear = hero (red, heavier); slab recedes (grey dashed);
    # unit Gaussian is a faint reference guide.
    hsl = np.array(d["counts_slab"], float) / (d["n_slab"] * bw)
    hre = np.array(d["counts_rect"], float) / (d["n_rect"] * bw)
    h_slab = axL.step(ctr, np.where(hsl > 0, hsl, np.nan), where="mid",
                      color=fs.MUTE, lw=1.3, ls=(0, (4.0, 1.6)), zorder=2)[0]
    h_rect = axL.step(ctr, np.where(hre > 0, hre, np.nan), where="mid",
                      color=fs.COLORS["rectilinear"], lw=1.7, ls="-", zorder=3)[0]
    g = np.exp(-0.5 * bins ** 2) / np.sqrt(2 * np.pi)
    h_g = axL.plot(bins, g, color=fs.FAINT, lw=1.0, ls=":", zorder=1)[0]
    for lo, hi in ((-6, -3), (3, 6)):
        axL.axvspan(lo, hi, color=fs.FAINT, alpha=0.10, lw=0, zorder=0)

    # the single figure legend (hero first), placed clear of curves and annotation
    axL.legend([h_rect, h_slab, h_g],
               ["rectilinear lattice", "matched slab", "unit Gaussian"],
               loc="lower center", fontsize=9, labelspacing=0.35,
               handlelength=1.9, borderaxespad=0.6)

    # void-channel spike pointer: mark the t=0 mass on-canvas (p_void from the JSON)
    pv = d.get("p_void", 0.342) * 100
    axL.annotate(rf"void channel, {pv:.0f}\%", xy=(0.07, 1.9), xytext=(1.35, 4.3),
                 fontsize=9, color=fs.COLORS["rectilinear"], ha="left", va="center",
                 arrowprops=dict(arrowstyle="-", color=fs.COLORS["rectilinear"],
                                 lw=0.8, shrinkA=2, shrinkB=2))

    axL.set_yscale("log")
    axL.set_ylim(1e-5, 1.5e1)
    axL.set_xlim(-6, 6)
    axL.set_xticks([-6, -4, -2, 0, 2, 4, 6])
    axL.set_xlabel(r"$\theta_\mathrm{kink}/\sigma$")
    axL.set_ylabel(r"probability density")
    fs.despine(axL)
    fs.panel_title(axL, "a", "Same width, different shape")

    # ---- (b) tail fraction P(|theta| > k sigma) at integer k ---------------------
    ks = np.array([2, 3, 4, 5])
    ts = np.array([d["tail_slab"][str(k)] for k in ks])
    tr = np.array([d["tail_rect"][str(k)] for k in ks])
    # same fixed encoding as (a): red = rectilinear (hero), grey = matched slab,
    # slab DASHED as in (a) so the pair separates in greyscale too (G4)
    axR.plot(ks, ts, "o", ls=(0, (4.0, 1.6)), color=fs.MUTE, ms=5.5, mfc="white",
             mew=1.2, lw=1.3, zorder=2)
    axR.plot(ks, tr, "-" + fs.MARKERS["rectilinear"], color=fs.COLORS["rectilinear"],
             ms=5.5, mfc="white", mew=1.4, lw=1.7, zorder=3)
    # required per-k ratio labels -- consistently centred just above each red marker
    for k, a, b in zip(ks, tr, ts):
        axR.annotate(rf"$\times{a/b:.1f}$", (k, a), textcoords="offset points",
                     xytext=(0, 9), ha="center", va="bottom", fontsize=9,
                     color=fs.COLORS["rectilinear"])

    axR.set_yscale("log")
    from matplotlib.ticker import FixedLocator, NullFormatter, FuncFormatter
    axR.yaxis.set_major_locator(FixedLocator([0.005, 0.01, 0.02, 0.05]))
    axR.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    axR.yaxis.set_minor_formatter(NullFormatter())
    axR.set_ylim(3.2e-3, 9.5e-2)
    axR.set_xlim(1.6, 5.4)
    axR.set_xticks(ks)
    axR.set_xlabel(r"$k$ (threshold in $\sigma$)")
    axR.set_ylabel(r"$P(|\theta_\mathrm{kink}| > k\sigma)$")
    fs.despine(axR)
    fs.panel_title(axR, "b", "Excess large-angle kinks")

    fs.save(fig, "fig_impact")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
