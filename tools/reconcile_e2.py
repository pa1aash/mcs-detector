#!/usr/bin/env python3
"""reconcile_e2.py -- P7/D7 reconciliation: compare the freshly regenerated
e2_results_combined.json (+ per-energy fits) against the committed headline values.
Any quantity outside its committed CI fails the check (exit 1) and must be logged in
DEVIATIONS.md before any manuscript number changes.
"""
from __future__ import annotations
import json, math, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
A = os.path.join(ROOT, "data", "analysis")

COMMITTED = {
    "pooled_slope": (-0.983, -1.052, -0.914),
    "slope_200":    (-0.951, -1.040, -0.866),
    "slope_500":    (-1.001, -1.095, -0.908),
    "width_maxdev_pct_lt": 8.0,          # max |ratio-1| must stay below 8%
    "a_eff_200": 3.8743e-6,              # +-1.5%
    "a_eff_500": 7.2405e-7,
}


def main():
    J = json.load(open(os.path.join(A, "e2_results_combined.json")))
    fit = J["exponent_fit"]
    ok = True

    def check(name, val, lo, hi):
        nonlocal ok
        good = lo <= val <= hi
        ok &= good
        print(f"{'PASS' if good else 'FAIL'}  {name}: fresh {val:+.3f}  committed CI [{lo:+.3f}, {hi:+.3f}]")

    check("pooled slope", fit["slope"], COMMITTED["pooled_slope"][1], COMMITTED["pooled_slope"][2])
    for E, key in ((200, "slope_200"), (500, "slope_500")):
        p = os.path.join(A, f"e2_results_{E}.json")
        if os.path.exists(p):
            s = json.load(open(p))["exponent_fit"]["slope"]
            check(f"slope {E}", s, COMMITTED[key][1], COMMITTED[key][2])

    rows = [r for r in J["rows"] if r["topology"] != "diamond" and not r["thin"]
            and r.get("f_width") and math.isfinite(r["f_width"])]
    maxdev = max(abs(r["f_width"] / r["f_designed"] - 1) for r in rows) * 100
    good = maxdev < COMMITTED["width_maxdev_pct_lt"]
    ok &= good
    print(f"{'PASS' if good else 'FAIL'}  width max deviation: {maxdev:.1f}% (< 8% required)")

    for E in (200, 500):
        a = J["a_eff_of_E"][str(E)]
        ref = COMMITTED[f"a_eff_{E}"]
        good = abs(a / ref - 1) < 0.015
        ok &= good
        print(f"{'PASS' if good else 'FAIL'}  a_eff({E}): fresh {a:.4e}  committed {ref:.4e}  ({(a/ref-1)*100:+.2f}%)")

    print("RECONCILIATION:", "GREEN" if ok else "RED — log in DEVIATIONS.md before touching manuscript numbers")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
