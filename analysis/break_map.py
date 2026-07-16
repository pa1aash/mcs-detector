#!/usr/bin/env python3
"""break_map.py -- S5(rebuilt) Phase C: map the homogeneous-approximation break
with the validated transport-aware ray-tracer on ANALYTIC geometry (no Geant4).

Per (topology, cell, momentum): run the transport SDE -> Var(t_actual) (wandering)
and Var(t_straight) (fixed transverse). The break observable is
    ratio(cell,p) = Var(t_actual)/Var(t_straight) = Dk4_TA/Dk4_LI
(the S3 a(p) and D4 factor cancel in the ratio). ratio->1 where straight-chord
holds; the break = the cell where ratio departs from 1 beyond THRESHOLD.
Also records the absolute Dk4_TA, Dk4_LI (3 a^2 Var (1+d4)) for the Geant4 cross
-check, and the analytic lateral scale y_rms(p).

CPU-light, analytic geometry -> no resolution/step/MSC artifact.
"""
from __future__ import annotations
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "geom"))
import transport_raytrace as tr  # noqa: E402

OUT = os.path.join(ROOT, "data", "analysis")
L = 10.0
F = 0.40
THRESHOLD = 0.10          # break = |1 - ratio| > 10%
PRINT_FDM = 0.5
PRINT_SLA = 0.05
CELLS = [3.0, 2.0, 1.0, 0.7, 0.5, 0.3, 0.2, 0.15, 0.1, 0.07, 0.05, 0.03]
MOMENTA = [100, 200, 500, 1000]
TOPOS = ["rectilinear", "gyroid"]
# S3 D4 amplitude renorm per momentum (STATE locked): kappa_M''/(6a^2)
D4 = {100: -0.209, 200: -0.167, 500: -0.193, 1000: -0.196}
N_PROTON = 20000
STEPS_PER_CELL = 20
N_SEED = 2


def ratio_at(name, cell, p):
    a = tr.A_OF_P[p]
    rs, vas, vss = [], [], []
    for s in range(N_SEED):
        r = tr.transport_trace(name, cell, F, a, L=L, n_proton=N_PROTON,
                               steps_per_cell=STEPS_PER_CELL,
                               rng=np.random.default_rng(100 * p + s))
        va = np.var(r["t_actual"], ddof=1)
        vs = np.var(r["t_straight"], ddof=1)
        rs.append(va / vs); vas.append(va); vss.append(vs)
    return (float(np.mean(rs)), float(np.std(rs)),
            float(np.mean(vas)), float(np.mean(vss)),
            float(np.std(r["y"], ddof=1)))


def break_cell(curve):
    """curve: list of (cell, ratio) sorted DESC in cell. Return the cell where
    ratio first crosses (1-THRESHOLD) going to smaller cells (linear interp in
    log-cell)."""
    c = sorted(curve, key=lambda t: -t[0])
    target = 1.0 - THRESHOLD
    for i in range(1, len(c)):
        (c0, r0), (c1, r1) = c[i - 1], c[i]
        if r0 >= target and r1 < target:
            lr = (np.log(c0) + (target - r0) * (np.log(c1) - np.log(c0)) / (r1 - r0))
            return float(np.exp(lr))
    if c[-1][1] >= target:
        return None             # never breaks down to smallest cell
    return float(c[0][0])       # already broken at largest cell (shouldn't happen)


def main():
    grid = {}
    for name in TOPOS:
        for p in MOMENTA:
            for cell in CELLS:
                rmean, rsd, va, vs, yrms = ratio_at(name, cell, p)
                d4 = D4[p]; a = tr.A_OF_P[p]
                grid[(name, p, cell)] = dict(
                    ratio=rmean, ratio_sd=rsd, var_act=va, var_str=vs,
                    Dk4_TA=3 * a ** 2 * va * (1 + d4),
                    Dk4_LI=3 * a ** 2 * vs * (1 + d4), y_rms_mm=yrms)
            print(f"{name:11s} p={p:4d}: " + " ".join(
                f"{c:g}:{grid[(name,p,c)]['ratio']:.2f}" for c in CELLS))

    breaks = {}
    for name in TOPOS:
        breaks[name] = {}
        for p in MOMENTA:
            curve = [(c, grid[(name, p, c)]["ratio"]) for c in CELLS]
            bc = break_cell(curve)
            breaks[name][p] = bc

    # y_rms(p) analytic (full-solid-path); lattice y_rms reported from the runs
    yrms_solid = {p: float(np.sqrt(tr.A_OF_P[p] * L ** 3 / 3.0)) for p in MOMENTA}

    out = dict(
        method="transport_aware_raytrace_analytic_geometry",
        threshold=THRESHOLD, L_mm=L, f=F, print_FDM_mm=PRINT_FDM,
        print_SLA_mm=PRINT_SLA, cells_mm=CELLS, momenta_MeV=MOMENTA,
        n_proton=N_PROTON, steps_per_cell=STEPS_PER_CELL, n_seed=N_SEED,
        cell_break_mm={n: {str(p): breaks[n][p] for p in MOMENTA} for n in TOPOS},
        y_rms_solid_mm={str(p): yrms_solid[p] for p in MOMENTA},
        grid={f"{n}|{p}|{c}": grid[(n, p, c)] for (n, p, c) in grid},
    )
    with open(os.path.join(OUT, "e0_break.json"), "w") as f:
        json.dump(out, f, indent=2)

    print("\n=== cell_break(momentum) [mm], threshold 10% ===")
    for name in TOPOS:
        print(f"  {name}: " + "  ".join(
            f"{p}MeV:{(f'{breaks[name][p]*1e3:.0f}um' if breaks[name][p] else 'sub-0.01mm')}"
            for p in MOMENTA))
    return out, breaks


if __name__ == "__main__":
    main()
