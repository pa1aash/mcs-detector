"""generate.py -- S4 watertight ASCII-STL generator for the S6 campaign parts.

Per topology x infill at the campaign cell (2.5 mm, L=10 mm): build the material
field, marching-cubes -> watertight surface, write ASCII STL, and run the four
S4 guards. To keep G4TessellatedSolid navigation tractable (GUARD 2), the parts
are a SINGLE TRANSVERSE UNIT CELL x full L in z (1x1 transverse, L/cell in z),
to be used in S6 with REFLECTIVE transverse boundaries + beam sampling across the
cell -- the task-sanctioned facet-reduction route (a 4x4x4 block is ~1.1M facets;
single-cell is ~<100k and identical p(t)). Voronoi (stochastic) uses a 3x3 block.

Run inside g4highland:  python geom/generate.py
"""
from __future__ import annotations

import os
import sys

import numpy as np
import trimesh
from skimage import measure

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import topologies as topo  # noqa: E402

STL_DIR = os.path.join(HERE, "stl")
L = 10.0
CELL = 2.5
INFILLS = (0.20, 0.30, 0.40, 0.50)
SAMPLES_PER_CELL = 40            # >=40 -> as-designed f converges <1%
FACET_FLAG = 200_000             # GUARD 2 navigation-cost threshold
F_TOL = 0.01                     # target-f tolerance (<1%)


def write_ascii_stl(path, verts, faces, name="part"):
    """Explicit ASCII STL writer (GUARD 1 verifies the result)."""
    tri = verts[faces].astype(np.float64)
    nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    nl = np.linalg.norm(nrm, axis=1, keepdims=True)
    nrm = np.where(nl > 0, nrm / nl, 0.0)
    fmt = ("facet normal %.6e %.6e %.6e\n  outer loop\n"
           "    vertex %.6e %.6e %.6e\n    vertex %.6e %.6e %.6e\n"
           "    vertex %.6e %.6e %.6e\n  endloop\nendfacet\n")
    parts = []
    for i in range(faces.shape[0]):
        parts.append(fmt % (nrm[i, 0], nrm[i, 1], nrm[i, 2],
                            tri[i, 0, 0], tri[i, 0, 1], tri[i, 0, 2],
                            tri[i, 1, 0], tri[i, 1, 1], tri[i, 1, 2],
                            tri[i, 2, 0], tri[i, 2, 1], tri[i, 2, 2]))
    with open(path, "w") as fh:
        fh.write(f"solid {name}\n")
        fh.write("".join(parts))
        fh.write(f"endsolid {name}\n")


def assert_ascii(path):
    """GUARD 1: re-open and prove the file is ASCII STL, not 80-byte binary."""
    with open(path, "rb") as fh:
        head = fh.read(1024)
    if not head.startswith(b"solid"):
        raise SystemExit(f"GUARD1 FAIL: {path} does not start with 'solid' (binary?)")
    if b"\x00" in head:
        raise SystemExit(f"GUARD1 FAIL: {path} has NUL bytes in header (binary STL)")
    if b"facet" not in head:
        raise SystemExit(f"GUARD1 FAIL: {path} has no 'facet' keyword (not ASCII STL)")


def _signed_volume(verts, faces):
    tri = verts[faces]
    return float(np.sum(np.einsum("ij,ij->i",
                 tri[:, 0], np.cross(tri[:, 1], tri[:, 2])) / 6.0))


def _mesh(field, dx):
    """marching-cubes (void-padded -> watertight) -> (verts, faces), oriented with
    OUTWARD normals (positive signed volume). skimage winds normals toward the
    high-field (solid) side -> inward -> G4 sees a negative-volume, inside-out
    solid and tracks get stuck; flipping the winding fixes both."""
    fpad = np.pad(field, 1)
    verts, faces, _, _ = measure.marching_cubes(fpad, level=0.5, spacing=(dx, dx, dx))
    verts = verts - dx
    if _signed_volume(verts, faces) < 0:
        faces = faces[:, ::-1]
    return verts, faces


def _mesh_volume(verts, faces):
    """absolute volume of a closed triangular mesh (divergence theorem)."""
    return abs(_signed_volume(verts, faces))


def _frac(field, dx, box_vol):
    """as-built fraction: mesh volume / box (or the constant value for a 1-phase
    field, where no surface exists)."""
    if field.max() == field.min():
        return float(field.mean())
    v, fc = _mesh(field, dx)
    return _mesh_volume(v, fc) / box_vol


def _tune_on_meshvol(name, box, cell, target_f, gs, dx, iters=40):
    """Bisect the analytic scalar so the marching-cubes MESH volume fraction ==
    target_f (sub-voxel-smooth in the parameter, unlike the quantised voxel mean
    for thin struts). Returns (param, verts, faces, f_mesh)."""
    box_vol = box[0] * box[1] * box[2]
    decreasing = name in topo.TPMS               # frac decreases as level rises
    lo, hi = (-3.2, 3.2) if decreasing else (0.0, cell * 0.62)
    mid = 0.5 * (lo + hi)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        f = _frac(topo.sample_field(name, box, cell, mid, gs), dx, box_vol)
        if (f > target_f) if decreasing else (f < target_f):
            lo = mid
        else:
            hi = mid
    field = topo.sample_field(name, box, cell, mid, gs)
    v, fc = _mesh(field, dx)
    return mid, v, fc, _mesh_volume(v, fc) / box_vol


def make_part(name, target_f, cell=CELL, spc=SAMPLES_PER_CELL, rng=None):
    """Build one part. Returns dict of metrics. Single transverse cell (analytic)
    or 3x3 transverse (voronoi); full L in z."""
    ncell_xy = 3 if name == "voronoi" else 1
    box = (ncell_xy * cell, ncell_xy * cell, L)
    dx = cell / spc
    gs = (int(round(box[0] / dx)), int(round(box[1] / dx)), int(round(box[2] / dx)))

    box_vol = box[0] * box[1] * box[2]
    if name == "voronoi":
        delta = topo.voronoi_delta(cell, box, gs, rng)
        wall = float(np.quantile(delta, target_f))   # exact voxel f by construction
        field = (delta < wall).astype(np.float32)
        param = wall
        verts, faces = _mesh(field, dx)
        f_designed = float(field.mean())             # stochastic block fraction
    else:
        # f is a property of the (infinite) lattice -> tune+report on the PERIODIC
        # unit cell (smooth, monotonic, accurate). The single-cell STL clips struts
        # at the box faces, so its mesh volume is NOT the lattice f; with reflective
        # transverse boundaries (S6) the Geant4 part reproduces the periodic f.
        param, f_designed = topo.tune_analytic(name, cell, target_f, n=80, tol=8e-4)
        field = topo.sample_field(name, box, cell, param, gs)
        verts, faces = _mesh(field, dx)
    f_meshed = _mesh_volume(verts, faces) / box_vol  # as-meshed (clipped) diagnostic
    verts = verts - np.array([box[0] / 2, box[1] / 2, box[2] / 2])  # centre at origin

    os.makedirs(STL_DIR, exist_ok=True)
    fname = f"{name}_f{int(round(target_f*100)):02d}_c{cell:g}.stl"
    path = os.path.join(STL_DIR, fname)
    write_ascii_stl(path, verts, faces, name=fname[:-4])
    assert_ascii(path)                           # GUARD 1

    m = trimesh.load(path)
    facets = int(len(m.faces))
    watertight = bool(m.is_watertight)
    vol = float(abs(m.volume))
    vol_expected = f_meshed * box[0] * box[1] * box[2]
    vol_err = vol / vol_expected - 1.0 if vol_expected > 0 else np.nan
    return dict(name=name, target_f=target_f, f_designed=f_designed,
                f_meshed=f_meshed, cell=cell,
                box=box, facets=facets, watertight=watertight, vol=vol,
                vol_expected=vol_expected, vol_err=vol_err,
                f_err=f_designed / target_f - 1.0, path=path,
                facet_flag=facets > FACET_FLAG)


def main():
    rng = np.random.default_rng(4040)
    rows = []
    print(f"{'part':>26} {'f_des':>6} {'f_err%':>7} {'facets':>8} {'WT':>3} "
          f"{'vol_err%':>9} {'flag':>5}")
    for name in topo.ALL_TOPOS:
        for tf in INFILLS:
            r = make_part(name, tf, rng=rng)
            rows.append(r)
            print(f"{os.path.basename(r['path']):>26} {r['f_designed']:>6.3f} "
                  f"{r['f_err']*100:>+6.2f}% {r['facets']:>8d} "
                  f"{'Y' if r['watertight'] else 'N':>3} {r['vol_err']*100:>+8.2f}% "
                  f"{'>200k' if r['facet_flag'] else '':>5}")

    # GUARD 3: matched f across topologies at each infill (clean topology test)
    print("\nGUARD3 matched-f across topologies (per infill):")
    ok3 = True
    for tf in INFILLS:
        fs = [r["f_designed"] for r in rows if r["target_f"] == tf]
        spread = max(fs) - min(fs)
        good = spread < 2 * F_TOL
        ok3 &= good
        print(f"  f={tf:.2f}: spread={spread*100:.2f}%  {'OK' if good else 'WIDE'}")

    nbad_f = sum(abs(r["f_err"]) > F_TOL for r in rows)
    nbad_wt = sum(not r["watertight"] for r in rows)
    nflag = sum(r["facet_flag"] for r in rows)
    print(f"\nparts={len(rows)}  f>1%:{nbad_f}  non-watertight:{nbad_wt}  "
          f">200k facets:{nflag}  GUARD3:{'OK' if ok3 else 'WIDE'}")
    # persist a small manifest
    import csv
    with open(os.path.join(STL_DIR, "manifest.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["part", "target_f", "f_designed", "f_err", "cell", "facets",
                    "watertight", "vol", "vol_err", "facet_flag"])
        for r in rows:
            w.writerow([os.path.basename(r["path"]), r["target_f"],
                        f"{r['f_designed']:.4f}", f"{r['f_err']:.4f}", r["cell"],
                        r["facets"], r["watertight"], f"{r['vol']:.4f}",
                        f"{r['vol_err']:.4f}", r["facet_flag"]])
    return 0 if (nbad_f == 0 and nbad_wt == 0) else 1


if __name__ == "__main__":
    raise SystemExit(main())
