#!/usr/bin/env python3
"""m2_tilt.py -- M2 orientation dependence: N_eff(theta) rocking curve + (when tilted
Geant4 runs exist) the Delta_kappa4(theta) collapse check.

N_eff(theta) is computed from TILTED chords through the analytic material field: a beam
tilted by theta about x drifts in y as z*tan(theta), with path element ds = dz/cos(theta)
and total path L/cos(theta). N_eff(theta) = f(1-f)(L/cos)^2 / Var(t_tilted). Reuses
topologies.chi_analytic + raytrace.stats_from_chi. Writes results/M2/neff_theta.json.
"""
from __future__ import annotations
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "geom"))
sys.path.insert(0, os.path.join(ROOT, "analysis", "lib"))
import topologies as topo   # noqa
import raytrace as rt       # noqa

OUT = os.path.join(ROOT, "results", "M2")
CELL, L, F = 2.5, 10.0, 0.30
TILTS = [0, 5, 15, 30]
FINE = [0, 1, 2, 3]


def neff_tilted(name, theta_deg, nxy=80, chunk=200):
    th = np.radians(theta_deg); ct = np.cos(th); tt = np.tan(th)
    p, _ = topo.tune_analytic(name, CELL, F, n=80, tol=8e-4)
    dz = CELL / 80.0
    nz = int(round(L / dz))
    z = (np.arange(nz) + 0.5) * dz
    g = (np.arange(nxy) + 0.5) / nxy * CELL
    XX, YY = np.meshgrid(g, g, indexing="ij")
    xr, yr = XX.ravel(), YY.ravel()
    n = xr.size
    chi = np.empty((n, nz), np.uint8)
    for i in range(0, n, chunk):
        xc = xr[i:i + chunk][:, None]
        yt = yr[i:i + chunk][:, None] + z[None, :] * tt      # tilted chord: y drifts with z
        chi[i:i + chunk] = topo.chi_analytic(name, xc, yt, z[None, :], CELL, p)
    s = rt.stats_from_chi(chi, dz / ct, L / ct, corr_frac=0.9)  # ds=dz/cos, path=L/cos
    return dict(theta=theta_deg, f=float(s.f), var_t=float(s.var_t),
                N_eff=float(s.N_eff_exact), tpath=float(L / ct))


def main():
    os.makedirs(OUT, exist_ok=True)
    out = {"cell": CELL, "L": L, "f": F, "rock": {}, "fine": {}}
    for name in ("rectilinear", "gyroid"):
        out["rock"][name] = [neff_tilted(name, t) for t in TILTS]
        print(f"  {name} rocking N_eff(theta):")
        for r in out["rock"][name]:
            print(f"    theta={r['theta']:2d} deg  f={r['f']:.3f}  N_eff={r['N_eff']:.2f}  "
                  f"path={r['tpath']:.2f} mm")
    out["fine"]["rectilinear"] = [neff_tilted("rectilinear", t) for t in FINE]
    json.dump(out, open(os.path.join(OUT, "neff_theta.json"), "w"), indent=1)
    print(f"wrote {OUT}/neff_theta.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
