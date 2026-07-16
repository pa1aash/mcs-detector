#!/usr/bin/env python3
"""calib_kappa4.py -- S6 Phase 0c: bootstrap-kappa4-CI-vs-N calibration.

Given ONE full-statistics run (Voronoi 0.30 @ 500 MeV), subsample at increasing N
and report kappa4, its bootstrap 95% CI width, and the fractional CI width
(CI/|kappa4|) vs N. This shows the N at which kappa4 is "stable enough" to resolve
the N_eff trend, sizing the per-config statistics budget. Uses the LOCKED W(500).
"""
from __future__ import annotations
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HERE, "lib"))
import kink_stats as ks  # noqa
RUNS = os.path.join(ROOT, "data", "runs")
W500 = 16.22e-3


def main():
    tag = sys.argv[1] if len(sys.argv) > 1 else "camp_voronoi_f30_E500"
    rp = os.path.join(RUNS, tag + ".root")
    if not os.path.exists(rp):
        print("missing", rp); return 1
    a = ks.load_run(rp).angles                      # 2 entries per primary (x,y)
    Ntot = a.size
    print(f"{tag}: {Ntot} angle samples (= 2 x primaries)")
    print(f"{'N_prim':>9} {'kappa4':>11} {'CI95 width':>11} {'CI/|k4|':>8}")
    rng = np.random.default_rng(7)
    out = []
    for Nprim in [50000, 100000, 200000, 500000, Ntot // 2]:
        n = min(2 * Nprim, Ntot)
        sub = a[rng.permutation(Ntot)[:n]]
        _, k4 = ks.cumulants_in_window(sub, W500)
        k4m, lo, hi = ks.bootstrap_kappa4(sub, W500, n_boot=400)
        ci = hi - lo
        frac = ci / abs(k4m) if k4m else float("nan")
        out.append((Nprim, k4m, ci, frac))
        print(f"{Nprim:>9} {k4m:>11.3e} {ci:>11.2e} {frac:>7.1%}")
    # crude target: N where fractional CI < 30%
    good = [o[0] for o in out if o[3] < 0.30]
    print(f"\nN_prim for fractional kappa4 CI < 30%: "
          f"{min(good) if good else '> '+str(out[-1][0])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
