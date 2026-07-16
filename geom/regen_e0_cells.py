"""regen_e0_cells.py -- S5.5 Task 8: regenerate the E0 single-cubic-cell STLs.

Adds the sub-printable break cells {0.2, 0.03, 0.01} mm that S5 lacked, and
re-cuts ALL e0_pair cells at higher marching-cubes resolution (spc=80 vs the S4
spc=40) so the rectilinear strut volume is resolved to <0.2% (the S4 spc=40 cut
under-resolved the cylindrical struts by ~1%, failing the <1% geometry-
consistency gate; gyroid was already exact). One cubic unit cell per (topology,
cell), centred at the origin, ASCII STL (GUARD1 re-assert), to be tiled by
mcs_sim (/mcs/det/geom stl + ntile) to fill L=10 mm.

Run inside g4highland:  python geom/regen_e0_cells.py
"""
from __future__ import annotations
import csv
import os
import sys

import numpy as np
import trimesh

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import topologies as topo  # noqa: E402
from generate import write_ascii_stl, assert_ascii, _mesh  # noqa: E402

F = 0.40
SPC = 80                                    # samples per cell (S4 used 40)
CELLS = (0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0)  # full S5 sweep
TOPOS = ("rectilinear", "gyroid")
STL_DIR = os.path.join(HERE, "stl", "e0_pair")


def regen(name, cell, spc=SPC):
    box = (cell, cell, cell)
    dx = cell / spc
    gs = tuple(int(round(b / dx)) for b in box)
    p, _ = topo.tune_analytic(name, cell, F, n=80, tol=8e-4)
    field = topo.sample_field(name, box, cell, p, gs)
    verts, faces = _mesh(field, dx)
    f_mesh = float(abs(np.sum(np.einsum(
        "ij,ij->i", verts[faces][:, 0],
        np.cross(verts[faces][:, 1], verts[faces][:, 2])) / 6.0)) / cell ** 3)
    verts = verts - np.array([cell / 2, cell / 2, cell / 2])
    os.makedirs(STL_DIR, exist_ok=True)
    path = os.path.join(STL_DIR, f"{name}_f{int(F*100):02d}_c{cell:g}.stl")
    write_ascii_stl(path, verts, faces, name=os.path.basename(path)[:-4])
    assert_ascii(path)                                       # GUARD1
    m = trimesh.load(path)
    return dict(name=name, cell=cell, facets=int(len(m.faces)),
                watertight=bool(m.is_watertight), f_meshed=f_mesh,
                path=os.path.relpath(path, os.path.dirname(HERE)))


def main():
    rows = []
    print(f"{'part':>34} {'facets':>8} {'WT':>3} {'f_mesh':>7} {'dF_analytic':>11}")
    f_ana = {}
    for name in TOPOS:
        # tune AND evaluate at the same cell (=1) so the strut radius / TPMS level
        # is consistent with the sampling grid (absolute-length params).
        p, fa = topo.tune_analytic(name, 1.0, F, n=80, tol=8e-4)
        g = (np.arange(200) + 0.5) / 200
        X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
        f_ana[name] = float(topo.chi_analytic(name, X, Y, Z, 1.0, p).mean())
    for cell in CELLS:
        for name in TOPOS:
            r = regen(name, cell)
            r["f_analytic"] = f_ana[name]
            r["dF"] = r["f_meshed"] / f_ana[name] - 1.0
            rows.append(r)
            print(f"{os.path.basename(r['path']):>34} {r['facets']:>8d} "
                  f"{'Y' if r['watertight'] else 'N':>3} {r['f_meshed']:>7.4f} "
                  f"{r['dF']*100:>+10.2f}%")
    nbadf = sum(abs(r["dF"]) > 0.01 for r in rows)
    nbadwt = sum(not r["watertight"] for r in rows)
    with open(os.path.join(STL_DIR, "manifest.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["part", "topology", "cell", "facets", "watertight",
                    "f_meshed", "f_analytic", "dF"])
        for r in rows:
            w.writerow([os.path.basename(r["path"]), r["name"], r["cell"],
                        r["facets"], r["watertight"], f"{r['f_meshed']:.4f}",
                        f"{r['f_analytic']:.4f}", f"{r['dF']:.4f}"])
    print(f"\nparts={len(rows)}  |dF|>1%:{nbadf}  non-watertight:{nbadwt}  "
          f"spc={SPC}  -> manifest {os.path.join(STL_DIR,'manifest.csv')}")
    return 0 if (nbadf == 0 and nbadwt == 0) else 1


if __name__ == "__main__":
    raise SystemExit(main())
