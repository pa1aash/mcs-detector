#!/usr/bin/env python3
"""make_pt_hist.py -- derive the compact input for the p(t) figure.

The plotted quantity is the *straight-chord geometric* line integral through the
as-built f=0.40 voxel fields, not Geant4's wandering ``tpla`` path.  This matches
the definition in Eq. (2.2), keeps 0 <= t <= L exactly, and gives literal endpoint
masses P(t=0) and P(t=L).  Output: data/analysis/pt_hist.json.
"""
from __future__ import annotations
import csv, json, os
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VOX = os.path.join(ROOT, "data", "geom_stats", "voxel")
OUT = os.path.join(ROOT, "data", "analysis", "pt_hist.json")

BINS = np.linspace(0.0, 10.0, 101)      # 0.1 mm bins over the full depth L = 10 mm
out = {"bins_mm": BINS.tolist(), "f_designed": 0.40, "L_mm": 10.0,
       "quantity": "straight-chord line integral through the as-built voxel field",
       "topologies": {}}
for t in ("rectilinear", "gyroid", "voronoi"):
    stem = os.path.join(VOX, f"{t}_f40_c2.5_camp_vox")
    nx, ny, nz, voxel, cell = open(stem + ".raw.meta").read().split()
    nx, ny, nz = int(nx), int(ny), int(nz)
    voxel, cell = float(voxel), float(cell)
    chi = np.fromfile(stem + ".raw", dtype=np.uint8).reshape(nx, ny, nz)
    chords = chi.sum(axis=2, dtype=float).ravel() * voxel
    counts, _ = np.histogram(chords, bins=BINS)
    probability = counts / counts.sum()
    h = probability / np.diff(BINS)
    out["topologies"][t] = {
        "density": h.tolist(),
        "p_void": float(np.mean(chords == 0.0)),
        "p_solid": float(np.mean(np.isclose(chords, out["L_mm"], atol=voxel / 2))),
        "mean_t_mm": float(chords.mean()),
        "var_t_mm2": float(chords.var()),
        "n": int(chords.size),
        "source": os.path.relpath(stem + ".raw", ROOT),
        "voxel_mm": voxel,
        "cell_mm": cell,
    }
    print(f"{t:12s} n={chords.size:.3g}  <t>={chords.mean():.3f} mm  "
          f"Var(t)={chords.var():.3f} mm^2  P0={np.mean(chords == 0.0):.4f}  "
          f"PL={out['topologies'][t]['p_solid']:.4f}")

# realised-mean check for the width-invariance residuals (referee point): for every
# available campaign run, the as-built mean line-integral relative to nominal fL.
# The kappa2 width ratio tracks <t>_realised/(fL); the Fig-3 spread is as-built fill,
# not a width-law violation.
import glob
try:
    import uproot
except ImportError:
    uproot = None
out["mean_t_over_fL"] = {}
for path in ([] if uproot is None else
             sorted(glob.glob(os.path.join(ROOT, "data", "runs", "camp_*_E*.root")))):
    tag = os.path.basename(path)[:-5]
    try:
        f_nom = int(tag.split("_f")[1].split("_")[0]) / 100.0
        tpla = uproot.open(path)["kinks"]["tpla"].array(library="np")
        out["mean_t_over_fL"][tag] = round(float(tpla.mean() / (f_nom * 10.0)), 5)
    except Exception as e:            # unreadable/partial file: skip, never fail
        print("skip", tag, e)

json.dump(out, open(OUT, "w"))
print("wrote", OUT, f"({len(out['mean_t_over_fL'])} mean-t entries)")
