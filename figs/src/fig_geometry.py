#!/usr/bin/env python3
"""Solid-path-length maps of the five campaign geometries.

For each topology the panel shows t(x,y): the PLA line-integral a beam-direction
(+z) proton crosses as a function of its transverse entry point, through one
2.5 mm cell of the *unfiltered* as-built f=0.40 campaign voxel field.  This is
exactly the geometric input to the scale-mixture law -- p(t) (Figure fig:pt) is
its histogram, and the spatial variance of each map fixes Var(t) and N_eff.

White = open channel (t~0); dark = a full 2.5 mm of solid.  For the four
periodic lattices the map mean is f times the cell depth (1.0 mm at f=0.40); the
central crop is cropped directly from the raw field with no smoothing.
"""
from __future__ import annotations

import os
import shutil
import subprocess

import numpy as np
import figstyle as fs

fs.set_style()
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap, Normalize

VOX = os.path.join(fs.ROOT, "data", "geom_stats", "voxel")
TOPOS = ["rectilinear", "schwarzp", "gyroid", "diamond", "voronoi"]
CELL = 2.5
BEAM = "#b23a48"

# void = white, more solid = darker ink (same solid/void reading as before).
CMAP = LinearSegmentedColormap.from_list(
    "solidpath", ["#ffffff", "#aeb7c2", "#4a5766", "#232c37"])


def load_cell_cube(name: str) -> tuple[np.ndarray, float]:
    """Aligned central nominal-cell crop from the campaign raw field."""
    stem = os.path.join(VOX, f"{name}_f40_c2.5_camp_vox")
    meta, raw = stem + ".raw.meta", stem + ".raw"
    if not os.path.exists(raw):
        raise FileNotFoundError(
            f"missing {raw}; regenerate with "
            f"`python geom/make_campaign_voxels.py {name} 0.40`"
        )
    Nx, Ny, Nz, voxel, _ = open(meta).read().split()
    Nx, Ny, Nz = int(Nx), int(Ny), int(Nz)
    voxel = float(voxel)
    n = int(round(CELL / voxel))
    if not np.isclose(n * voxel, CELL, atol=2e-5):
        raise ValueError(f"{name}: cell/voxel is not integral ({CELL}/{voxel})")
    chi = np.fromfile(raw, dtype=np.uint8).reshape(Nx, Ny, Nz)
    nc = (Nx // n, Ny // n, Nz // n)
    x0, y0, z0 = tuple(((k - 1) // 2) * n for k in nc)
    cube = chi[x0:x0 + n, y0:y0 + n, z0:z0 + n]
    if cube.shape != (n, n, n):
        raise ValueError(f"{name}: incomplete central crop {cube.shape}")
    return np.ascontiguousarray(cube), voxel


def main() -> int:
    with plt.rc_context({"figure.constrained_layout.use": False}):
        fig = plt.figure(figsize=(fs.FULL, fs.FULL * 0.72))
    gs = fig.add_gridspec(2, 3, hspace=0.22, wspace=0.12)
    axes = [fig.add_subplot(gs[i // 3, i % 3]) for i in range(5)]

    for i, (ax, name) in enumerate(zip(axes, TOPOS)):
        cube, voxel = load_cell_cube(name)
        tmap = cube.sum(axis=2) * voxel               # solid path along +z [mm]
        ax.imshow(tmap.T, cmap=CMAP, vmin=0, vmax=CELL, origin="lower",
                  interpolation="nearest", extent=[0, CELL, 0, CELL],
                  rasterized=True)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_color(fs.HERO); s.set_linewidth(0.6)
        ax.text(0.5, 1.06, fs.TOPO_LABEL[name], transform=ax.transAxes,
                ha="center", va="bottom", fontsize=10.0, color=fs.HERO)
        ax.text(0.03, 0.965, rf"\textbf{{({chr(97 + i)})}}",
                transform=ax.transAxes, ha="left", va="top",
                fontsize=10.5, color=fs.HERO)

    # legend / colourbar / beam key
    lx = fig.add_subplot(gs[1, 2]); lx.axis("off")
    lx.set_xlim(0, 1); lx.set_ylim(0, 1)
    lx.text(0.02, 0.93, r"map: solid path $t(x,y)$", fontsize=9,
            color=fs.HERO, va="top")
    lx.text(0.02, 0.80, r"a $+z$ proton crosses", fontsize=9,
            color=fs.HERO, va="top")
    lx.annotate("", xy=(0.93, 0.94), xytext=(0.93, 0.68),
                arrowprops=dict(arrowstyle="-|>", mutation_scale=12,
                                color=BEAM, lw=1.35))
    lx.text(0.885, 0.81, r"beam $+z$", fontsize=9, color=BEAM,
            ha="right", va="center")
    cax = lx.inset_axes([0.04, 0.58, 0.80, 0.075])
    sm = cm.ScalarMappable(cmap=CMAP, norm=Normalize(0, CELL))
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_ticks([0, 1.0, 2.5])
    cb.ax.tick_params(labelsize=8.5, length=2.5)
    cb.set_label(r"solid path length $t$ [\si{mm}]", fontsize=9, labelpad=2)
    lx.text(0.02, 0.30, r"white: open channel ($t\!\approx\!0$)", fontsize=8.4,
            color=fs.HERO, va="center")
    lx.text(0.02, 0.18, r"dark: full \SI{2.5}{mm} solid", fontsize=8.4,
            color=fs.HERO, va="center")
    lx.text(0.02, 0.06, r"\SI{2.5}{mm} cell, $f=0.40$ nominal", fontsize=8.4,
            color=fs.HERO, va="center")

    with plt.rc_context({"figure.constrained_layout.use": False}):
        fig.subplots_adjust(left=0.02, right=0.985, bottom=0.055, top=0.92)

    out = os.path.join(fs.FIGS, "fig_geometry.pdf")
    fig.savefig(out, facecolor="white", transparent=False, dpi=600)
    os.makedirs(fs.PREVIEW, exist_ok=True)
    pdftoppm = shutil.which("pdftoppm")
    if pdftoppm:
        subprocess.run([
            pdftoppm, "-png", "-r", "300", "-singlefile", out,
            os.path.join(fs.PREVIEW, "fig_geometry"),
        ], check=True)
    print("wrote", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
