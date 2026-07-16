"""make_voronoi_voxel.py -- S5.5 Task 5 prep: one Voronoi foam realisation as
BOTH a voxel field (for /mcs/det/geom voxel) and a tessellated STL (for
/mcs/det/geom stl), so the voxel-vs-tessellated kink agreement can be checked on
the SAME geometry. Voronoi is the theory-critical stochastic geometry (it anchors
the L/2 l_int scaling), and its ~750k-facet tessellation is why the S4 decision
was to voxelise it -- this validates that the voxelisation is artifact-free.

Outputs (data/geom_stats/voxel/):
  voronoi_block.raw   : N^3 uint8 (ix*N*N+iy*N+iz), 1=PLA 0=void
  voronoi_block.raw.meta : "N voxel_mm cell_mm"
  voronoi_block.stl   : ASCII STL (marching cubes of the SAME field)

Run inside g4highland:  python geom/make_voronoi_voxel.py
"""
from __future__ import annotations
import os
import sys

import numpy as np
import trimesh

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import topologies as topo  # noqa: E402
from generate import write_ascii_stl, assert_ascii, _mesh  # noqa: E402

OUT = os.path.join(os.path.dirname(HERE), "data", "geom_stats", "voxel")
S = 5.0            # block side (mm)
CELL = 2.5         # Voronoi cell (mm) -- printable
N = 128            # voxel grid edge -> voxel = S/N ~ 0.039 mm
TARGET_F = 0.40
SEED = 4040


def main():
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(SEED)
    box = (S, S, S)
    gs = (N, N, N)
    chi, wall, f = topo.voronoi_field(CELL, box, TARGET_F, gs, rng)  # chi[ix,iy,iz]
    voxel = S / N

    # voxel field (C-order ravel = ix*N*N + iy*N + iz, matches VoxelParam)
    raw = np.ascontiguousarray(chi.astype(np.uint8))
    raw.tofile(os.path.join(OUT, "voronoi_block.raw"))
    with open(os.path.join(OUT, "voronoi_block.raw.meta"), "w") as fh:
        fh.write(f"{N} {voxel:.6f} {CELL:.6f}\n")

    # tessellated STL from the SAME field (centred at origin)
    field = chi.astype(np.float32)
    verts, faces = _mesh(field, voxel)
    verts = verts - np.array([S / 2, S / 2, S / 2])
    stl = os.path.join(OUT, "voronoi_block.stl")
    write_ascii_stl(stl, verts, faces, name="voronoi_block")
    assert_ascii(stl)
    m = trimesh.load(stl)

    print(f"Voronoi block: S={S} mm  cell={CELL} mm  N={N}  voxel={voxel:.4f} mm")
    print(f"  voxel solid fraction f = {f:.4f}  (target {TARGET_F})")
    print(f"  STL: facets={len(m.faces)}  watertight={m.is_watertight}  "
          f"f_meshed={abs(m.volume)/S**3:.4f}")
    print(f"  wrote {OUT}/voronoi_block.{{raw,raw.meta,stl}}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
