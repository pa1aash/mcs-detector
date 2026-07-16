#!/usr/bin/env python3
"""geom_consistency.py -- S5.5 Task 4 (the STOP gate) + helpers.

Fires straight, zero-scattering GEANTINO chords through the production mcs_sim
STL-lattice build (CADMesh unit cell, periodically tiled) and confirms the
Geant4-SEEN solid fraction f and the chord variance Var(t) match the ray-tracer
values (data/geom_stats/<topology>_f40.npz, summary.csv) the line-integral
prediction is built on. If they disagree beyond tolerance the detector is NOT
building the geometry the theory assumes -> STOP.

Geantinos transport through the geometry with no physics, so the per-primary PLA
path length (ntuple column `tpla`) is exactly t = integral chi dz for a straight
ray -> <t>/L = f and Var(t) directly comparable to the ray-tracer.
"""
from __future__ import annotations
import json
import os
import subprocess
import sys

import numpy as np
import uproot

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BIN = os.path.join(ROOT, "sim", "build", "mcs_sim")
GEOM = os.path.join(ROOT, "data", "geom_stats")
STL = os.path.join(ROOT, "geom", "stl", "e0_pair")
TMP = "/tmp/gc"
L = 10.0
F_TOL = 0.01        # <1% on f  (the gate)
VAR_TOL = 0.03      # <3% on Var(t) (ray-tracer estimator itself carries ~2%)
NEVT = 200000
NTILE = 5
THREADS = 6
os.makedirs(TMP, exist_ok=True)


def run_geantino_probe(name, cell):
    stl = os.path.join(STL, f"{name}_f40_c{cell:g}.stl")
    out = os.path.join(TMP, f"probe_{name}_c{cell:g}")
    mac = out + ".mac"
    with open(mac, "w") as f:
        f.write(f"""/mcs/det/geom stl
/mcs/det/stl {stl}
/mcs/det/cell {cell} mm
/mcs/det/depth {L} mm
/mcs/det/ntile {NTILE}
/mcs/det/topology {name}
/mcs/det/infill 0.40
/mcs/det/maxStep 0.1 mm
/run/setCut 0.05 mm
/run/initialize
/gun/particle geantino
/gun/energy 1 GeV
/mcs/gun/spotXY {cell} mm
/analysis/setFileName {out}
/run/beamOn {NEVT}
""")
    subprocess.run([BIN, mac, str(THREADS)], check=True,
                   stdout=open(out + ".log", "w"), stderr=subprocess.STDOUT)
    tpla = uproot.open(out + ".root")["kinks"]["tpla"].array(library="np")
    return tpla


def ray_targets(name):
    d = np.load(os.path.join(GEOM, f"{name}_f40.npz"))
    return float(d["t_mean"]), float(d["var_t"]), float(d["f"])


def main():
    if not os.path.exists(BIN):
        sys.exit("build mcs_sim first")
    cell = 2.0          # printable cell for the consistency check
    rows, gate = [], True
    for name in ("rectilinear", "gyroid"):
        tpla = run_geantino_probe(name, cell)
        t_mean, var_t = float(tpla.mean()), float(tpla.var(ddof=1))
        f_g4 = t_mean / L
        tgt_tmean, tgt_var, tgt_f = ray_targets(name)
        df = f_g4 / tgt_f - 1.0
        dvar = var_t / tgt_var - 1.0
        ok = abs(df) < F_TOL and abs(dvar) < VAR_TOL
        gate &= ok
        rows.append(dict(name=name, cell=cell, n=int(tpla.size),
                         f_g4=f_g4, f_ray=tgt_f, df=df,
                         var_g4=var_t, var_ray=tgt_var, dvar=dvar,
                         tmean_g4=t_mean, tmean_ray=tgt_tmean, ok=ok))
    out = dict(check="geometry_consistency", cell_mm=cell, ntile=NTILE,
               n_geantino=NEVT, f_tol=F_TOL, var_tol=VAR_TOL,
               rows=rows, gate="PASS" if gate else "FAIL")
    with open(os.path.join(ROOT, "data", "analysis",
                           "geom_consistency.json"), "w") as f:
        json.dump(out, f, indent=2)
    for r in rows:
        print(f"{r['name']:12s} c{r['cell']:g}: f_G4={r['f_g4']:.4f} "
              f"f_ray={r['f_ray']:.4f} ({r['df']*100:+.2f}%)  "
              f"Var_G4={r['var_g4']:.3f} Var_ray={r['var_ray']:.3f} "
              f"({r['dvar']*100:+.2f}%)  {'OK' if r['ok'] else 'FAIL'}")
    print(f"\nGEOMETRY CONSISTENCY: {'PASS' if gate else 'FAIL (STOP)'}")
    return 0 if gate else 1


if __name__ == "__main__":
    raise SystemExit(main())
