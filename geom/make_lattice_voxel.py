"""make_lattice_voxel.py -- S5.5 / S5-rebuilt: voxelise an analytic lattice
(rect/gyroid) into the .raw/.meta format VoxelParam reads. Non-cubic L-deep
blocks for the Phase-B4/D Geant4 validation: ncell_xy transverse cells (so a
wandering proton sees correct periodic neighbours) x ncell_z = depth/cell cells
deep. Voxel grid (no coincident inter-cell surfaces) => 0 stuck tracks (the
S5.5 decision: the proton campaign builds VOXELS).

meta format: "Nx Ny Nz voxel_mm cell_mm".
Usage: python geom/make_lattice_voxel.py <topology> <cell_mm> [ncell_xy] [depth_mm] [res]
"""
from __future__ import annotations
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import topologies as topo  # noqa: E402

OUT = os.path.join(os.path.dirname(HERE), "data", "geom_stats", "voxel")
F = 0.40


def make(name, cell, ncell_xy=3, depth=10.0, res=32):
    os.makedirs(OUT, exist_ok=True)
    voxel = cell / res
    Nx = Ny = int(round(ncell_xy * cell / voxel))
    Nz = int(round(depth / voxel))
    p, _ = topo.tune_analytic(name, cell, F, n=80, tol=8e-4)
    gx = (np.arange(Nx) + 0.5) * voxel
    gz = (np.arange(Nz) + 0.5) * voxel
    X, Y, Z = np.meshgrid(gx, gx, gz, indexing="ij")
    chi = topo.chi_analytic(name, X, Y, Z, cell, p).astype(np.uint8)  # [ix,iy,iz]
    raw = np.ascontiguousarray(chi)                                   # ix*Ny*Nz+iy*Nz+iz
    stem = os.path.join(OUT, f"{name}_f40_c{cell:g}_L{depth:g}_vox")
    raw.tofile(stem + ".raw")
    with open(stem + ".raw.meta", "w") as fh:
        fh.write(f"{Nx} {Ny} {Nz} {voxel:.6f} {cell:.6f}\n")
    print(f"{name} c{cell} L{depth}: {Nx}x{Ny}x{Nz} voxel={voxel:.4f}mm "
          f"f_vox={chi.mean():.4f} -> {os.path.relpath(stem,os.path.dirname(HERE))}.raw")
    return stem, float(chi.mean())


if __name__ == "__main__":
    name = sys.argv[1] if len(sys.argv) > 1 else "rectilinear"
    cell = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0
    ncell_xy = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    depth = float(sys.argv[4]) if len(sys.argv) > 4 else 10.0
    res = int(sys.argv[5]) if len(sys.argv) > 5 else 32
    make(name, cell, ncell_xy, depth, res)
