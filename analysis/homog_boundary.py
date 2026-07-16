#!/usr/bin/env python3
"""homog_boundary.py -- S6 Stage-2C: the HOMOGENEOUS-approximation boundary cell_homog(p),
defined from the excess kurtosis gamma2 (the error the homogeneous/Highland model cannot
reproduce), distinct from the line-integral boundary cell_break(p) (Var ratio = 0.90).

For the projected-angle scale mixture, gamma2 = kappa4/kappa2^2 = 3 Var(t)/<t>^2, using the
ACTUAL (wandering) path length t the proton accumulates. The validated transport ray-tracer
(geom/transport_raytrace.py, Geant4-confirmed incl. 100 MeV in S5-verify) gives t_actual.
  - Large cell: gamma2 -> geometry maximum (proton resolves chords; homogeneous FAILS).
  - Fine cell : gamma2 -> 0 (MCS wandering averages the structure; homogeneous valid).
cell_homog(p) = the cell where gamma2 crosses a STATED threshold (default 0.10; also report
the angular-resolution-motivated value). cell_break(p) = Var(t_act)/Var(t_str) = 0.90.

C3 numerical-floor check: at fine cells, vary the transport step (steps_per_cell) and N to
test whether the residual gamma2 is physical (stable) or a transport-tool floor (shrinks).

Outputs data/analysis/homog_boundary.json. Run inside g4highland.
"""
from __future__ import annotations
import argparse, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "geom"))
import transport_raytrace as tr  # noqa

OUT = os.path.join(ROOT, "data", "analysis")
L = 10.0
F = 0.40
TOPOS = ["rectilinear", "gyroid"]
MOMENTA = [100, 200, 500, 1000]
# gamma2 stays near its geometry maximum until the cell approaches the lateral
# wandering scale y_rms(p) (~tens of um), so the sweep must reach few-um cells to
# capture the gamma2 threshold crossing (cell_homog << cell_break).
# Feasible cell range: gamma2 stays ~= its geometry max until the cell approaches the
# few-um wandering scale, so the gamma2-threshold crossing (cell_homog) is at ~um. The
# fine-cell transport cost ~ 1/cell, so the sweep stops at 5 um and cell_homog at the
# strict threshold is reported by interpolation/extrapolation of the measured curve.
CELLS = [3.0, 2.5, 2.0, 1.0, 0.5, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01, 0.007, 0.005]
GAMMA2_THRESHS = [1.0, 0.5, 0.1]   # report cell_homog at several stated thresholds
BREAK_RATIO = 0.90             # cell_break: Var(t_act)/Var(t_str) departs 10%
N_PROTON = 8000                # Var(t)/gamma2 well-determined at 8k rays
STEPS_PER_CELL = 12
N_SEED = 2


def metrics_at(name, cell, p, n_proton=N_PROTON, steps_per_cell=STEPS_PER_CELL, n_seed=N_SEED):
    """Return mean gamma2, its sd, and the Var ratio, over n_seed transport runs."""
    a = tr.A_OF_P[p]
    g2s, ratios = [], []
    for s in range(n_seed):
        r = tr.transport_trace(name, cell, F, a, L=L, n_proton=n_proton,
                               steps_per_cell=steps_per_cell,
                               rng=np.random.default_rng(1000 * p + s))
        ta = r["t_actual"]
        g2 = 3.0 * np.var(ta, ddof=1) / ta.mean() ** 2
        g2s.append(float(g2))
        ratios.append(float(np.var(ta, ddof=1) / np.var(r["t_straight"], ddof=1)))
    return float(np.mean(g2s)), float(np.std(g2s)), float(np.mean(ratios))


def _cross(curve, target):
    """curve: list of (cell, metric) sorted DESC in cell. Cell where metric drops below
    `target` going to smaller cells (log-cell interpolation). If it never crosses within
    the swept range, LOG-LINEAR extrapolate from the two finest cells and flag it."""
    c = sorted(curve, key=lambda t: -t[0])
    for i in range(1, len(c)):
        (c0, m0), (c1, m1) = c[i - 1], c[i]
        if m0 >= target and m1 < target:
            lr = np.log(c0) + (target - m0) * (np.log(c1) - np.log(c0)) / (m1 - m0)
            return float(np.exp(lr)), "interp"
    # not reached within range -> extrapolate from the finest two points (gamma2 still
    # above target at the finest cell)
    (c0, m0), (c1, m1) = c[-2], c[-1]
    if m1 >= target and m1 < m0:
        lr = np.log(c0) + (target - m0) * (np.log(c1) - np.log(c0)) / (m1 - m0)
        return float(np.exp(lr)), "extrap(<%.3gmm)" % c1
    return None, "none"


def build_curves():
    grid, out = {}, {"thresholds_gamma2": GAMMA2_THRESHS, "break_ratio": BREAK_RATIO,
                     "L": L, "f": F, "cells": CELLS, "momenta": MOMENTA,
                     "cell_homog": {}, "cell_homog_method": {}, "cell_break": {},
                     "gamma2_curves": {}, "gamma2_max": {}}
    for name in TOPOS:
        for k in ("cell_homog", "cell_homog_method", "cell_break", "gamma2_curves",
                  "gamma2_max"):
            out[k][name] = {}
        for p in MOMENTA:
            curve_g2, curve_ratio = [], []
            for c in CELLS:
                g2, g2sd, ratio = metrics_at(name, c, p)
                grid[(name, c, p)] = (g2, ratio)
                curve_g2.append((c, g2)); curve_ratio.append((c, ratio))
            out["gamma2_max"][name][str(p)] = curve_g2[0][1]      # coarsest cell
            cb, _ = _cross(curve_ratio, BREAK_RATIO)
            out["cell_break"][name][str(p)] = cb
            out["cell_homog"][name][str(p)] = {}
            out["cell_homog_method"][name][str(p)] = {}
            for thr in GAMMA2_THRESHS:
                ch, method = _cross(curve_g2, thr)
                out["cell_homog"][name][str(p)][str(thr)] = ch
                out["cell_homog_method"][name][str(p)][str(thr)] = method
            out["gamma2_curves"][name][str(p)] = {f"{c:g}": grid[(name, c, p)][0] for c in CELLS}
    return out


def numerical_floor_check():
    """C3: at fine cells, vary steps_per_cell (dz resolution) and N; the residual gamma2
    is PHYSICAL if stable vs both, a transport FLOOR if it drifts. rect @ 100 MeV (the
    flagged case). Light: spc up to 24, N up to 16000 (the 5-h overrun came from spc=80)."""
    res = {}
    for cell in (0.05, 0.02, 0.01):
        row = {}
        for spc in (8, 12, 24):
            g2, sd, _ = metrics_at("rectilinear", cell, 100, n_proton=10000,
                                   steps_per_cell=spc, n_seed=2)
            row[f"spc{spc}"] = round(g2, 4)
        for npr in (8000, 16000):
            g2, sd, _ = metrics_at("rectilinear", cell, 100, n_proton=npr,
                                   steps_per_cell=12, n_seed=2)
            row[f"N{npr}"] = round(g2, 4)
        res[f"cell{cell}mm"] = row
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true", help="tiny run to validate plumbing")
    args = ap.parse_args()
    if args.smoke:
        for c in (2.5, 0.1):
            g2, sd, ratio = metrics_at("rectilinear", c, 200, n_proton=8000, n_seed=2)
            print(f"rect c={c} 200MeV: gamma2={g2:.4f}+/-{sd:.4f}  Var_ratio={ratio:.3f}")
        return 0
    out = build_curves()
    out["numerical_floor"] = numerical_floor_check()
    json.dump(out, open(os.path.join(OUT, "homog_boundary.json"), "w"), indent=1)
    print("=== gamma2(cell) at the printable 2.5 mm cell and fine cells (rect/gyroid) ===")
    for name in TOPOS:
        for p in MOMENTA:
            cur = out["gamma2_curves"][name][str(p)]
            print(f"  {name:11s} p={p:4d}: g2@2.5mm={cur.get('2.5',float('nan')):.2f} "
                  f"@0.1mm={cur.get('0.1',float('nan')):.2f} @0.01mm={cur.get('0.01',float('nan')):.2f}")
    print(f"\n=== cell_homog(p) [mm] at gamma2 in {GAMMA2_THRESHS}  vs  cell_break (ratio {BREAK_RATIO}) ===")
    for name in TOPOS:
        for p in MOMENTA:
            cb = out["cell_break"][name][str(p)]
            chs = out["cell_homog"][name][str(p)]
            s = " ".join(f"g2<{t}:{('%.4f'%chs[str(t)]) if chs[str(t)] else 'none'}" for t in GAMMA2_THRESHS)
            print(f"  {name:11s} p={p:4d}: cell_break={('%.3f'%cb) if cb else 'none':>6}  {s}")
    print("\n=== C3 numerical-floor check (rect @ 100 MeV; gamma2 stable vs spc & N -> physical) ===")
    for cell, row in out["numerical_floor"].items():
        print(f"  {cell}: " + "  ".join(f"{k}={v}" for k, v in row.items()))
    print("wrote data/analysis/homog_boundary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
