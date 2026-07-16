#!/usr/bin/env python3
"""m2_collapse.py -- M2 orientation collapse: does Delta_kappa4(theta) follow N_eff(theta)^-1?

For each tilted Geant4 run, Delta_kappa4 = kappa4(struct,W) - <kappa_M(tpla)> (all-order floor
from the M2 PLA solids), and C = Delta_kappa4 / (3 a_eff^2 f(1-f) L_path^2) with L_path=L/cos(theta).
Theory: C = 1/N_eff(theta). N_eff(theta) from the tilted-chord integrator (m2_tilt.neff_tilted) at
the exact run angle. Writes results/M2/collapse_theta.json + F4 data.
"""
from __future__ import annotations
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "analysis", "lib"))
sys.path.insert(0, HERE)
import kink_stats as ks   # noqa
import m2_tilt            # noqa  (neff_tilted)

M2 = os.path.join(ROOT, "data", "runs_m2")
CAMP = os.path.join(ROOT, "data", "runs")     # campaign solids (full t range, D8)
OUT = os.path.join(ROOT, "results", "M2")
W, X0, L, F, E = 37.84e-3, 315.423, 10.0, 0.30, 200
RUNS = {"rectilinear": [0, 1, 2, 3, 5, 10], "gyroid": [0, 5, 10]}
FLOOR_THICKS = (2, 3, 4, 5, 8, 16)            # covers the full tpla range (max = L = 10)


def k2k4(tag, base=None):
    a = ks.load_run(os.path.join(base or M2, tag + ".root")).angles
    return (*ks.cumulants_in_window(a, W),)


def load_tpla(tag):
    import uproot
    with uproot.open(os.path.join(M2, tag + ".root")) as f:
        return f["kinks"]["tpla"].array(library="np")


def build_floor():
    # D8: fit over the full traversed-path range using the campaign solid controls
    # (t up to 16 mm); the old M2-local fit (t<=5, clipped) under-covered the 13.5%
    # of rectilinear paths with tpla > 5 mm and biased dk4 high by ~30%.
    kt, kv = [0.0], [0.0]
    for t in FLOOR_THICKS:
        _, k4 = k2k4(f"solid_E{E}_t{t}", base=CAMP)
        kt.append(float(t)); kv.append(k4)
    kt, kv = np.array(kt), np.array(kv)
    b, c = np.linalg.lstsq(np.vstack([kt, kt ** 2]).T, kv, rcond=None)[0]
    tmax = kt.max()
    return lambda tt: b * np.clip(tt, 0, tmax) + c * np.clip(tt, 0, tmax) ** 2


def a_eff():
    return float(np.mean([k2k4(f"solid_E{E}_t{t}", base=CAMP)[0] / t for t in (2, 3, 4, 5)]))


def main():
    os.makedirs(OUT, exist_ok=True)
    kM = build_floor(); aeff = a_eff()
    print(f"a_eff(200) = {aeff:.4e}")
    rows = []
    for topo, angles in RUNS.items():
        for th in angles:
            tag = f"m2_{topo}_tilt{th}"
            if not os.path.exists(os.path.join(M2, tag + ".root")):
                continue
            k2, k4 = k2k4(tag)
            tpla = load_tpla(tag)
            floor_mean = float(np.mean(kM(tpla)))
            dk4 = float(k4 - floor_mean)
            # bootstrap 95% CI on kappa4 (>=2000 resamples, Ground Rule 13); floor held
            # fixed (its sampling error is subdominant), matching the e2 convention
            a = ks.load_run(os.path.join(M2, tag + ".root")).angles
            _, k4_lo, k4_hi = ks.bootstrap_kappa4(a, W, n_boot=2000, seed=1234)
            dk4_lo, dk4_hi = float(k4_lo - floor_mean), float(k4_hi - floor_mean)
            ct = np.cos(np.radians(th))
            Lp = L / ct
            g = m2_tilt.neff_tilted(topo, th)          # tilted-chord N_eff(theta)
            neff = g["N_eff"]
            denom = 3 * aeff ** 2 * F * (1 - F) * Lp ** 2
            C = dk4 / denom
            rows.append(dict(topo=topo, theta=th, N_eff=neff, dk4=dk4,
                             dk4_lo=dk4_lo, dk4_hi=dk4_hi, C=C,
                             C_theory=1.0 / neff, ratio=C * neff,
                             ratio_lo=dk4_lo / denom * neff,
                             ratio_hi=dk4_hi / denom * neff,
                             tpla=float(tpla.mean()), Lpath=Lp))
            print(f"  {topo:11s} th={th:2d}  N_eff={neff:7.2f}  dk4={dk4:.3e}  "
                  f"C={C:.3e}  1/N_eff={1/neff:.3e}  C*N_eff={C*neff:.2f}")
    json.dump({"a_eff": aeff, "rows": rows}, open(os.path.join(OUT, "collapse_theta.json"), "w"),
              indent=1)
    # summary: does C track 1/N_eff? (C*N_eff ~ prefactor k, ~const across theta if collapse holds)
    kvals = [r["ratio"] for r in rows if r["N_eff"] < 50 and r["dk4"] > 0]  # exclude near-degenerate
    if kvals:
        print(f"\nC*N_eff (prefactor k) over resolved tilts: mean={np.mean(kvals):.2f} "
              f"sd={np.std(kvals):.2f}  -> collapse holds if ~const (theta=0 anchor {rows[0]['ratio']:.2f})")
    print(f"wrote {OUT}/collapse_theta.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
