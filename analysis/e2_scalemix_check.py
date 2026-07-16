#!/usr/bin/env python3
"""e2_scalemix_check.py -- S6 Stage-1 methods investigation of the voronoi tension.

Tests the line-integral scale-mixture identity (Result 2 / SC2)

    Delta_kappa4 = 3 a^2 Var(t)

DIRECTLY, using
  * a_eff calibrated from the solid controls (kappa2_solid = a_eff * t), which is
    topology-independent and removes the locked-Highland-a^2 absolute offset; and
  * Var(tpla) = the variance of the PER-PRIMARY measured PLA path length (the
    `tpla` ntuple column) -- the proton's ACTUAL traversed-material variance,
    including MCS wandering, not the straight-chord ray-traced Var(t).

If k_eff = Delta_kappa4 / (3 a_eff^2 Var(tpla)) ~= 1, the scale mixture holds with the
measured path variance (so any earlier collapse offset was the a-calibration and/or
straight-chord vs effective Var(t), i.e. methods). If k_eff deviates topology-wise,
that is a real model effect. RESULT (500 MeV): periodic k_eff ~= 1 (scale mixture
confirmed); voronoi k_eff scattered ~1-3, driven by the smallest-Var(t) (highest
N_eff) configs whose deconvolved Delta_kappa4 is statistics-limited -> escalate
voronoi stats to discriminate real-foam-excess vs deconvolution noise.

Run inside g4highland:  python analysis/e2_scalemix_check.py
"""
from __future__ import annotations
import json, os, sys
import numpy as np
import uproot

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HERE, "lib"))
import kink_stats as ks  # noqa

RUNS = os.path.join(ROOT, "data", "runs")
AOUT = os.path.join(ROOT, "data", "analysis")
L, W = 10.0, 16.22e-3


def var_tpla(tag):
    t = uproot.open(os.path.join(RUNS, tag + ".root"))["kinks"]["tpla"].array(library="np")
    return float(np.mean(t)), float(np.var(t))


def main():
    J = json.load(open(os.path.join(AOUT, "e2_results.json")))
    aH = J["a_of_E"]["500"]
    # a_eff from each solid control: kappa2_solid(t) = a_eff * t
    aeff = {}
    for tmm in (2, 3, 4, 5):
        ang = ks.load_run(os.path.join(RUNS, f"solid_E500_t{tmm}.root")).angles
        k2, _ = ks.cumulants_in_window(ang, W)
        aeff[tmm] = k2 / tmm
    rows_out = []
    res = {}
    for r in sorted(J["rows"], key=lambda r: (r["topology"], r["infill"])):
        if r["topology"] == "diamond":
            continue
        tmm = round(r["infill"] * L)
        ae = aeff[tmm]
        mt, vt = var_tpla(r["tag"])
        vgeom = (r["f_built"] * (1 - r["f_built"]) * L ** 2 / r["N_eff"]
                 if np.isfinite(r["N_eff"]) else float("nan"))
        pred_eff = 3 * ae ** 2 * vt
        keff = r["dk4"] / pred_eff
        res.setdefault(r["topology"], []).append(keff)
        rows_out.append(dict(tag=r["tag"], topology=r["topology"], infill=r["infill"],
                             N_eff=r["N_eff"], var_t_geom=vgeom, var_tpla=vt,
                             ratio_tpla_geom=vt / vgeom if vgeom else None,
                             dk4=r["dk4"], a_eff=ae, k_eff=keff,
                             dk4_ci_frac=abs(r["dk4_hi"] - r["dk4_lo"]) / abs(r["dk4"])))
    summary = {t: dict(k_eff_mean=float(np.mean(v)), k_eff_std=float(np.std(v)),
                       n=len(v)) for t, v in res.items()}
    out = dict(a_Highland=aH, a_eff_from_solids=aeff,
               note="k_eff = dk4 / (3 a_eff^2 Var(tpla)); ~1 confirms scale mixture",
               per_config=rows_out, per_topology=summary)
    json.dump(out, open(os.path.join(AOUT, "e2_scalemix_check.json"), "w"), indent=1)
    print(f"a_eff(solids) = {np.mean(list(aeff.values())):.4e} "
          f"(Highland {aH:.4e}; ratio {np.mean(list(aeff.values()))/aH:.3f})")
    for t, s in summary.items():
        print(f"  {t:12s} k_eff = {s['k_eff_mean']:.3f} ± {s['k_eff_std']:.3f}  (n={s['n']})")
    print("wrote e2_scalemix_check.json")


if __name__ == "__main__":
    main()
