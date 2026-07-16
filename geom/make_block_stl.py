"""make_block_stl.py -- S5-verify Phase 1: single WATERTIGHT block STL of an
analytic lattice (ncell_xy^2 transverse x depth/cell deep), as ONE marching-cubes
surface (NOT tiled unit cells -> 0 stuck, per the S5.5 finding). Correct f at fine
marching-cubes resolution, and -- unlike a voxel grid -- it does NOT force sub-floor
geometric steps: with maxStep=0.1 mm the proton runs in the S3-validated step
regime. This is the geometry rep for the 100 MeV break CONFIRMATION + step-regime
check (voxel rep is for the campaign's navigation robustness).

Usage: python geom/make_block_stl.py <topology> <cell_mm> [ncell_xy] [depth] [spc]
"""
from __future__ import annotations
import os, sys
import numpy as np
import trimesh

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import topologies as topo  # noqa: E402
from generate import write_ascii_stl, assert_ascii, _mesh  # noqa: E402

OUT = os.path.join(os.path.dirname(HERE), "geom", "stl", "blocks")
F = 0.40


def make(name, cell, ncell_xy=3, depth=10.0, spc=28):
    os.makedirs(OUT, exist_ok=True)
    dx = cell / spc
    bx = ncell_xy * cell
    box = (bx, bx, depth)
    gs = (int(round(bx / dx)), int(round(bx / dx)), int(round(depth / dx)))
    p, _ = topo.tune_analytic(name, cell, F, n=80, tol=8e-4)
    field = topo.sample_field(name, box, cell, p, gs)
    verts, faces = _mesh(field, dx)
    f_mesh = abs(float(np.sum(np.einsum("ij,ij->i", verts[faces][:, 0],
              np.cross(verts[faces][:, 1], verts[faces][:, 2])) / 6.0))) / (bx*bx*depth)
    verts = verts - np.array([bx/2, bx/2, depth/2])
    path = os.path.join(OUT, f"{name}_f40_c{cell:g}_block.stl")
    write_ascii_stl(path, verts, faces, name=os.path.basename(path)[:-4])
    assert_ascii(path)
    m = trimesh.load(path)
    print(f"{name} c{cell}: block {bx:g}x{bx:g}x{depth:g}mm spc={spc} "
          f"facets={len(m.faces)} WT={m.is_watertight} f_mesh={f_mesh:.4f} "
          f"-> {os.path.relpath(path, os.path.dirname(HERE))}")
    return path, len(m.faces), f_mesh, bool(m.is_watertight)


if __name__ == "__main__":
    name = sys.argv[1]; cell = float(sys.argv[2])
    ncell = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    depth = float(sys.argv[4]) if len(sys.argv) > 4 else 10.0
    spc = int(sys.argv[5]) if len(sys.argv) > 5 else 28
    make(name, cell, ncell, depth, spc)
