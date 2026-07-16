#!/usr/bin/env python3
"""make_foam_voxels.py -- S6 Stage-2 Phase-1: foam-scale (0.2 mm cell) voxel blocks for
the Geant4 gamma2 spot-check, at COARSE voxel resolution so the Geant4 step stays ABOVE
the S3-validated 0.02 mm maxStep floor (a 0.2 mm cell at the campaign res=32 would give
~6 um steps, in the corrupted MSC-vs-single-scattering regime that motivated the
transport tool). Two resolutions (res=4 -> 0.05 mm voxel/step; res=8 -> 0.025 mm) bracket
the step-regime guard: gamma2 must be voxel-resolution-stable to serve as a clean anchor.

Usage:  python geom/make_foam_voxels.py            # rect+gyroid, f0.4, res {4,8}
"""
from __future__ import annotations
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import topologies as topo  # noqa

OUT = os.path.join(os.path.dirname(HERE), "data", "geom_stats", "voxel")
CELL = 0.2          # foam-scale cell [mm]
L = 10.0
NCELL_XY = 5        # 1.0 mm wide block; spot 0.2 mm + tens-of-um wandering stays interior
INFILL = 0.40
TOPOS = ("rectilinear", "gyroid")
RESES = (4, 8)      # voxel = CELL/res = 0.05, 0.025 mm  (both > 0.02 mm step floor)


def stem_for(name, res):
    return os.path.join(OUT, f"{name}_f{int(INFILL*100):02d}_c{CELL:g}_res{res}_foam_vox")


def build(name, res):
    os.makedirs(OUT, exist_ok=True)
    voxel = CELL / res
    Nx = Ny = int(round(NCELL_XY * CELL / voxel))
    Nz = int(round(L / voxel))
    p, _ = topo.tune_analytic(name, CELL, INFILL, n=80, tol=8e-4)
    gx = (np.arange(Nx) + 0.5) * voxel
    gz = (np.arange(Nz) + 0.5) * voxel
    X, Y, Z = np.meshgrid(gx, gx, gz, indexing="ij")
    chi = topo.chi_analytic(name, X, Y, Z, CELL, p).astype(np.uint8)
    fach = float(chi.mean())
    stem = stem_for(name, res)
    np.ascontiguousarray(chi).tofile(stem + ".raw")
    with open(stem + ".raw.meta", "w") as fh:
        fh.write(f"{Nx} {Ny} {Nz} {voxel:.6f} {CELL:.6f}\n")
    return stem, fach, (Nx, Ny, Nz), voxel


def main():
    print(f"{'topo':12} {'res':>3} {'voxel[mm]':>9} {'f_ach':>6} {'grid':>14} {'Mvox':>6}")
    for name in TOPOS:
        for res in RESES:
            stem, fach, (Nx, Ny, Nz), voxel = build(name, res)
            print(f"{name:12} {res:>3} {voxel:9.4f} {fach:6.3f} {Nx}x{Ny}x{Nz:>4} "
                  f"{Nx*Ny*Nz/1e6:6.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
