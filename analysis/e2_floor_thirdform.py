#!/usr/bin/env python3
"""e2_floor_thirdform.py -- floor-model THIRD-FORM robustness check (Fable 4.2 / audit
FLOOR-3FORM).

STATUS: SCAFFOLD, NOT YET RUN (Phase-1 carry-forward). The two floor forms already tried are
(i) a free piecewise-linear interpolation (abandoned; left an interpolation-artifact
systematic) and (ii) the locked constrained quadratic kappa_M(t)=b*t+c*t^2 through the origin
(current, e2_analysis.build_floor). This script adds a THIRD physically-motivated form and
re-derives the periodic k_eff under it, to confirm the periodic result is robust to the floor
functional form (expected: <~1% change, since the periodic geometry signal is 14-27% of the
raw kappa_4 -- well above the floor-model-sensitive regime).

WHY IT COULD NOT BE RUN IN THE PHASE-1 (compute-free) PASS: the deconvolution needs the raw
solid-control kappa_4 values (from data/runs/solid_E*_t*.root) to FIT a new floor form, and
the per-proton tpla distributions (from the campaign .root files) to evaluate <kM_new(tpla)>.
Neither the .root files nor the six control kappa_4 values are committed to the repo (the rows
in e2_results_combined.json store tpla_mean and tpla_var but not the full tpla arrays, and no
file records the control kappa_4). So this must run in the g4highland environment on the box
where the .root files live. It exits cleanly with a CARRY-FORWARD message if they are absent.

Run (on the compute box, inside g4highland):
    python analysis/e2_floor_thirdform.py [--energies 200 500]

Third form implemented: a monotone constrained CUBIC kappa_M(t)=b*t+c*t^2+d*t^3 through the
origin, least-squares fit to the same controls t={2,3,4,5,8,16} mm as the locked quadratic,
plus (if scipy is available) a shape-preserving monotone PCHIP interpolant through the same
six points as an independent, non-parametric third form. Both are compared, per periodic
config, against the committed quadratic-floor k_eff. A >1% periodic shift is a STOP: it would
mean the floor model, not the geometry, is carrying the signal, and must be surfaced honestly
(Fable 4.2) rather than smoothed over.
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HERE, "lib"))
sys.path.insert(0, os.path.join(ROOT, "geom"))

RUNS = os.path.join(ROOT, "data", "runs")
PERIODIC = ("rectilinear", "schwarzp", "gyroid")
INFILLS = (0.20, 0.30, 0.40, 0.50)
CELL, L = 2.5, 10.0
SOLID_THICKS = (2, 3, 4, 5, 8, 16)


def _have_raw():
    """The solid controls at both energies must exist to fit any floor form."""
    need = [os.path.join(RUNS, f"solid_E{E}_t{t}.root")
            for E in (200, 500) for t in SOLID_THICKS]
    return all(os.path.exists(p) for p in need)


def fit_forms(kt, kv):
    """Return {name: callable kM(t)} for the three floor forms, all through the origin
    and fit to the SAME control points (kt, kv) with kM(0)=0 prepended by the caller."""
    forms = {}
    # (ii) locked quadratic  b t + c t^2   (the current e2_analysis form; reproduced here)
    Aq = np.vstack([kt, kt ** 2]).T
    bq, cq = np.linalg.lstsq(Aq, kv, rcond=None)[0]
    tmax = kt.max()
    forms["quadratic"] = lambda tt: (bq * np.clip(tt, 0, tmax) + cq * np.clip(tt, 0, tmax) ** 2)
    # (iii-a) constrained cubic  b t + c t^2 + d t^3
    Ac = np.vstack([kt, kt ** 2, kt ** 3]).T
    bc, cc, dc = np.linalg.lstsq(Ac, kv, rcond=None)[0]
    forms["cubic"] = lambda tt: (bc * np.clip(tt, 0, tmax) + cc * np.clip(tt, 0, tmax) ** 2
                                 + dc * np.clip(tt, 0, tmax) ** 3)
    # (iii-b) monotone shape-preserving PCHIP through the six points (non-parametric third form)
    try:
        from scipy.interpolate import PchipInterpolator
        pk = PchipInterpolator(kt, kv, extrapolate=True)
        forms["pchip"] = lambda tt: np.clip(pk(np.clip(tt, 0, tmax)), 0, None)
    except Exception as e:
        print(f"  (scipy PCHIP unavailable: {e}; cubic-only third form)")
    return forms


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--energies", nargs="+", type=int, default=[200, 500])
    args = ap.parse_args()

    if not _have_raw():
        print("CARRY-FORWARD: solid-control and campaign .root files are not present in this "
              "checkout (data/runs/*.root). The floor third-form refit needs the raw control "
              "kappa_4 values and per-proton tpla arrays, which are regenerable only on the "
              "compute box. Re-run there. No robustness number is asserted until this runs "
              "(Fable 4.2 / audit FLOOR-3FORM stays OPEN).")
        return 0

    import kink_stats as ks           # noqa
    import raytrace as rt             # noqa
    import topologies as topo         # noqa
    from e2_analysis import W, a_eff_of_E, geom_asbuilt, load_tpla  # reuse locked machinery

    stop = False
    for E in args.energies:
        aeff = a_eff_of_E(E)
        # fit floor forms from the solid controls at this energy
        kt, kv = [0.0], [0.0]
        for t in SOLID_THICKS:
            p = os.path.join(RUNS, f"solid_E{E}_t{t}.root")
            _, k4 = ks.cumulants_in_window(ks.load_run(p).angles, W[E])
            kt.append(float(t)); kv.append(k4)
        forms = fit_forms(np.array(kt), np.array(kv))
        print(f"\n=== E={E} MeV: periodic k_eff under each floor form ===")
        per = {name: [] for name in forms}
        for topo_name in PERIODIC:
            for infill in INFILLS:
                tag = f"camp_{topo_name}_f{int(round(infill*100)):02d}_E{E}"
                rp = os.path.join(RUNS, tag + ".root")
                if not os.path.exists(rp):
                    continue
                ang = ks.load_run(rp).angles
                _, k4s = ks.cumulants_in_window(ang, W[E])
                tpla = load_tpla(tag)
                tvar = float(np.var(tpla))
                for name, kM in forms.items():
                    dk4 = k4s - float(np.mean(kM(tpla)))
                    per[name].append(dk4 / (3 * aeff ** 2 * tvar) if tvar > 0 else np.nan)
        base = float(np.nanmean(per["quadratic"]))
        for name in forms:
            m = float(np.nanmean(per[name]))
            d = 100 * (m - base) / base if base else float("nan")
            flag = "  <-- STOP (>1%)" if abs(d) > 1.0 and name != "quadratic" else ""
            print(f"  {name:10s} periodic mean k_eff = {m:.4f}  (Delta vs quadratic {d:+.2f}%){flag}")
            if abs(d) > 1.0 and name != "quadratic":
                stop = True

    if stop:
        print("\nSTOP: a third floor form shifts the periodic k_eff by >1%. Surface this "
              "honestly (Fable 4.2) -- the floor model may be carrying the signal.")
        return 2
    print("\nOK: periodic k_eff is robust to the floor functional form (<1% across three "
          "forms). Add the comparison to Section 8.2 and close FLOOR-3FORM.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
