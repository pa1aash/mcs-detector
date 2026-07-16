"""make_campaign_voxels.py -- S6 Phase 0: build the campaign voxel geometries.

5 topologies x 4 infills at the campaign cell (2.5 mm, L=10 mm), as VoxelParam
.raw/.meta blocks (res=32 -> voxel 0.078 mm, the S5.5-validated resolution; 3x3
transverse cells so a wandering proton sampled across the central cell always has
periodic neighbours). N_eff is varied by topology+infill at FIXED cell.

Analytic topologies (rectilinear + TPMS gyroid/schwarzp/diamond) use tune_analytic;
voronoi uses voronoi_field (stochastic, seeded for reproducibility). Reports the
achieved f for the provenance/closure check.

Usage:
  python geom/make_campaign_voxels.py            # build all 20
  python geom/make_campaign_voxels.py <topo> <infill>          # build one @ 2.5 mm cell
  python geom/make_campaign_voxels.py <topo> <infill> <cell>   # build one @ a given cell (M3 wide-N_eff)
"""
from __future__ import annotations
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import topologies as topo  # noqa: E402

OUT = os.path.join(os.path.dirname(HERE), "data", "geom_stats", "voxel")
CELL = 2.5
L = 10.0
NCELL_XY = 3
# voxels per cell: rectilinear's thin cylindrical struts need res=48 to converge f
# to <1% at low infill (res=32 over-resolves them by ~4%); the smooth TPMS and the
# stochastic voronoi are <1% at res=32 (validated S5.5).
RES_BY_TOPO = {"rectilinear": 48}
RES_DEFAULT = 32
INFILLS = (0.20, 0.30, 0.40, 0.50)
TOPOS = ("rectilinear", "gyroid", "schwarzp", "diamond", "voronoi")
SEED = 6006                 # voronoi reproducibility


def stem_for(name, infill):
    return os.path.join(OUT, f"{name}_f{int(round(infill*100)):02d}_c{CELL:g}_camp_vox")


def build(name, infill):
    os.makedirs(OUT, exist_ok=True)
    voxel = CELL / RES_BY_TOPO.get(name, RES_DEFAULT)
    Nx = Ny = int(round(NCELL_XY * CELL / voxel))
    Nz = int(round(L / voxel))
    box = (Nx * voxel, Ny * voxel, Nz * voxel)
    if name == "voronoi":
        rng = np.random.default_rng(SEED + int(round(infill * 100)))
        chi, _, fach = topo.voronoi_field(CELL, box, infill, (Nx, Ny, Nz), rng)
        chi = chi.astype(np.uint8)
    else:
        p, _ = topo.tune_analytic(name, CELL, infill, n=80, tol=8e-4)
        gx = (np.arange(Nx) + 0.5) * voxel
        gz = (np.arange(Nz) + 0.5) * voxel
        X, Y, Z = np.meshgrid(gx, gx, gz, indexing="ij")
        chi = topo.chi_analytic(name, X, Y, Z, CELL, p).astype(np.uint8)
        fach = float(chi.mean())
    stem = stem_for(name, infill)
    np.ascontiguousarray(chi).tofile(stem + ".raw")
    with open(stem + ".raw.meta", "w") as fh:
        fh.write(f"{Nx} {Ny} {Nz} {voxel:.6f} {CELL:.6f}\n")
    return stem, fach, (Nx, Ny, Nz)


def main():
    global CELL
    rows = []
    if len(sys.argv) >= 4:                       # M3 wide-N_eff: override the cell size
        CELL = float(sys.argv[3])
    if len(sys.argv) >= 3:
        jobs = [(sys.argv[1], float(sys.argv[2]))]
    else:
        jobs = [(t, f) for t in TOPOS for f in INFILLS]
    print(f"{'topo':>12} {'f_des':>6} {'f_ach':>6} {'f_err%':>7} {'grid':>14} {'Mvox':>6}")
    for name, infill in jobs:
        stem, fach, (Nx, Ny, Nz) = build(name, infill)
        ferr = fach / infill - 1.0
        mvox = Nx * Ny * Nz / 1e6
        rows.append((name, infill, fach, ferr, mvox))
        print(f"{name:>12} {infill:>6.2f} {fach:>6.3f} {ferr*100:>+6.2f}% "
              f"{Nx}x{Ny}x{Nz:>4} {mvox:>6.2f}")
    nbad = sum(abs(r[3]) > 0.02 for r in rows)
    print(f"\n{len(rows)} configs; |f_err|>2%: {nbad}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
