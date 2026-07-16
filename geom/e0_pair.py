"""e0_pair.py -- S4 / E0 matched pair (rectilinear <-> gyroid) across the
cell-size sweep, cached for the S5 break hunt.

GUARD 4 finding (important): for DETERMINISTIC PERIODIC lattices the integral
correlation length l_int and N_eff are CELL-INDEPENDENT (identical cells stay
correlated over the whole depth L) and TOPOLOGY-FIXED. Measured here:
  rectilinear l_int ~ 4.2 mm (z-columns -> long correlation),
  gyroid      l_int ~ 2.2 mm (interpenetrating),
non-overlapping for ALL (f, cell). So the literally-specified "equal f AND equal
l_int" rect<->gyroid pair is geometrically INFEASIBLE -- l_int cannot be matched
by tuning cell (it does not scale) or f (ranges disjoint). FLAGGED.

Feasible E0 pair shipped instead: matched **f and cell size** (equal areal
density, equal feature scale / printability), topology varied -- the cleanest
break test available for two periodic topologies. The line-integral prediction
(each topology's own N_eff) is what S5 tests against full Geant4 transport as the
cell shrinks toward the lateral-scattering scale. e0_pair.csv reports BOTH
topologies' f, l_int, N_eff per cell so the (un)matched quantities are explicit.

If an EQUAL-N_eff contrast is wanted, S5 should pair a STOCHASTIC foam (N_eff
tunable via cell, ~1/cell) against a periodic lattice at matched f and N_eff --
periodic N_eff is not tunable. Noted for S5.

Run inside g4highland:  python geom/e0_pair.py
"""
from __future__ import annotations

import csv
import os
import sys

import numpy as np
from skimage import measure

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "analysis", "lib"))
import raytrace as rt          # noqa: E402
import topologies as topo      # noqa: E402
from generate import write_ascii_stl, assert_ascii, _mesh  # noqa: E402

L = 10.0
F = 0.40
MOMENTA_MeV = (100, 200, 500, 1000)          # S5 axis (geometry is momentum-independent)
SWEEP = (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0)
STL_SWEEP = (0.02, 0.05, 0.1, 0.5, 1.0, 2.0, 3.0)   # cache single-cell STLs at these
LINT_TOL = 0.05
OUT = os.path.join(os.path.dirname(HERE), "data", "geom_stats")
STL_DIR = os.path.join(HERE, "stl", "e0_pair")


def ref_stats(name, f, cells=(0.5, 1.0, 2.0, 2.5)):
    """l_int and N_eff are cell-independent for periodic lattices: measure at a
    few cells, return the mean and the spread (to prove cell-independence)."""
    p, fa = topo.tune_analytic(name, 2.5, f, n=80, tol=8e-4)
    lints, neffs, fs = [], [], []
    for c in cells:
        p_c, _ = topo.tune_analytic(name, c, f, n=80, tol=8e-4)
        chi, dz = topo.ray_chi_analytic(name, c, p_c, L, nxy=72, dz=c / 32.0)
        s = rt.stats_from_chi(chi, dz, L, corr_frac=0.9)
        lints.append(s.l_int); neffs.append(s.N_eff_exact); fs.append(s.f)
    return dict(l_int=float(np.mean(lints)), l_int_sd=float(np.std(lints)),
                N_eff=float(np.mean(neffs)), f=float(np.mean(fs)))


def cache_stl(name, cell, f, spc=40):
    """ONE cubic unit cell (cell^3), ASCII STL. For the sub-printable E0 cells a
    full-L part would be L/cell cells deep (~millions of facets / GB files), so
    S5 tiles this single cubic cell by REFLECTION in all three axes to fill the
    10 mm depth + transverse beam region."""
    box = (cell, cell, cell)
    dx = cell / spc
    gs = tuple(int(round(b / dx)) for b in box)
    p, _ = topo.tune_analytic(name, cell, f, n=64, tol=1e-3)
    field = topo.sample_field(name, box, cell, p, gs)
    verts, faces = _mesh(field, dx)
    verts = verts - np.array([box[0] / 2, box[1] / 2, box[2] / 2])
    os.makedirs(STL_DIR, exist_ok=True)
    path = os.path.join(STL_DIR, f"{name}_f{int(F*100):02d}_c{cell:g}.stl")
    write_ascii_stl(path, verts, faces, name=os.path.basename(path)[:-4])
    assert_ascii(path)
    return path, len(faces)


def main():
    os.makedirs(OUT, exist_ok=True)
    print("measuring cell-independence of l_int / N_eff (periodic lattices)...")
    rec = ref_stats("rectilinear", F)
    gyr = ref_stats("gyroid", F)
    print(f"  rectilinear: l_int={rec['l_int']:.3f}+-{rec['l_int_sd']:.3f} "
          f"N_eff={rec['N_eff']:.3f} f={rec['f']:.3f}")
    print(f"  gyroid:      l_int={gyr['l_int']:.3f}+-{gyr['l_int_sd']:.3f} "
          f"N_eff={gyr['N_eff']:.3f} f={gyr['f']:.3f}")
    l_int_match = abs(rec["l_int"] - gyr["l_int"]) / gyr["l_int"] < LINT_TOL
    print(f"  l_int matchable within {LINT_TOL:.0%}? {l_int_match}  "
          f"(ranges disjoint -> FLAG, ship matched-(f,cell) pair instead)")

    rows = []
    for cell in SWEEP:
        rows.append(dict(cell=cell, f_rect=rec["f"], f_gyroid=gyr["f"],
                         l_int_rect=rec["l_int"], l_int_gyroid=gyr["l_int"],
                         N_eff_rect=rec["N_eff"], N_eff_gyroid=gyr["N_eff"],
                         f_match_ok=abs(rec["f"] - gyr["f"]) < 0.01,
                         l_int_match_ok=l_int_match))

    with open(os.path.join(OUT, "e0_pair.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["cell_size", "f_rect", "f_gyroid", "l_int_rect", "l_int_gyroid",
                    "N_eff_rect", "N_eff_gyroid", "f_match_ok", "l_int_match_ok",
                    "note"])
        for r in rows:
            note = ("matched f+cell; l_int NOT matchable (periodic, disjoint)"
                    if not r["l_int_match_ok"] else "matched")
            w.writerow([r["cell"], f"{r['f_rect']:.4f}", f"{r['f_gyroid']:.4f}",
                        f"{r['l_int_rect']:.4f}", f"{r['l_int_gyroid']:.4f}",
                        f"{r['N_eff_rect']:.4f}", f"{r['N_eff_gyroid']:.4f}",
                        r["f_match_ok"], r["l_int_match_ok"], note])

    print("\ncaching single-cell E0 STLs (rect + gyroid)...")
    for cell in STL_SWEEP:
        for name in ("rectilinear", "gyroid"):
            path, nf = cache_stl(name, cell, F)
            print(f"  {os.path.basename(path):>28}  facets={nf}")

    # expected break range for S5
    print("\nExpected break (S5): proton lateral-scattering scale ~20-70 um at "
          "100-1000 MeV -> straight-chord prediction expected to break at cell "
          "~0.02-0.1 mm (SUB-PRINTABLE). Sweep momenta", MOMENTA_MeV, "MeV; lower "
          "momentum pushes the break to larger (more printable) cells. If the "
          "break never reaches printable scale (>~0.3 mm) even at 100 MeV, that is "
          "the design-rule result (structure is safely averageable when printable).")
    print(f"\ne0_pair.csv: {len(rows)} cells; f matched (0.40), l_int_match_ok="
          f"{l_int_match} (FLAGGED), STLs cached for {len(STL_SWEEP)} cells x 2 topo.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
