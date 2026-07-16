"""validate_estimator.py -- GATE 4 stop condition (S4 task 4(i)).

Generate a synthetic telegraph (exponential-chord Markov) field of known lambda_c,
run it through the raytrace.py l_int + Var(t) estimator, and confirm it reproduces
theory.py's closed-form B10 Var(t) and l_int = lambda_c within tolerance.

This validates the estimator + theory.py independent of any printed geometry.
If it disagrees badly, STOP -- the correlation-length estimator or theory.py is
wrong.  Run inside g4highland.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "analysis", "lib"))
import raytrace as rt          # noqa: E402
import theory as th            # noqa: E402

TOL = 0.05  # 5% relative tolerance on Var(t) and l_int (finite N + discretisation)


def main() -> int:
    rng = np.random.default_rng(20260620)
    L = 10.0
    cases = [(0.30, 0.20), (0.40, 0.50), (0.40, 1.00), (0.30, 2.00)]
    print(f"{'f':>5} {'lam_c':>6} {'Var_t meas':>12} {'Var_t thy':>12} "
          f"{'dVar%':>7} {'l_int meas':>11} {'l_int thy':>10} {'dl%':>7} "
          f"{'Neff_ex':>8} {'Neff_thy':>8}")
    worst = 0.0
    rows = []
    for f, lam_c in cases:
        dz = lam_c / 30.0
        chi = rt.make_telegraph(f, lam_c, L, dz, n_rays=40000, rng=rng)
        s = rt.stats_from_chi(chi, dz, L, corr_frac=0.9)
        var_thy = th.var_t_markov(f, L, lam_c)
        lint_thy = th.l_int_markov(lam_c)
        neff_thy = th.N_eff_exact(f, L, var_thy)
        dvar = s.var_t / var_thy - 1.0
        dl = s.l_int / lint_thy - 1.0
        worst = max(worst, abs(dvar), abs(dl))
        rows.append((f, lam_c, dvar, dl))
        print(f"{f:>5.2f} {lam_c:>6.2f} {s.var_t:>12.4e} {var_thy:>12.4e} "
              f"{dvar*100:>+6.2f}% {s.l_int:>11.4f} {lint_thy:>10.4f} "
              f"{dl*100:>+6.2f}% {s.N_eff_exact:>8.2f} {neff_thy:>8.2f}")
    ok = worst <= TOL
    print(f"\nworst |residual| = {worst*100:.2f}%   tolerance = {TOL*100:.0f}%   "
          f"-> {'PASS' if ok else 'FAIL (STOP)'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
