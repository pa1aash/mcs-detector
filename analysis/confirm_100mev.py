#!/usr/bin/env python3
"""confirm_100mev.py -- S5-verify Phase 1c/1b. Geant4-confirm the 100 MeV break.

Primary observable: Var(t_actual) -- the Geant4 proton PLA path length `tpla` IS
t_actual (the wandering solid path), compared directly to the transport SDE. The
break ratio Var(t_act)/Var(t_str) needs no floor subtraction. Secondary: Delta
kappa4 (G4 vs LI vs TA) with the locked W(100), solid t=4 control floor-scaled to
the lattice <t>. Step-regime (1b): rect 0.56 at voxel res20 (28um, above the S3
0.02mm floor) vs res40 (14um, below) -- if Var(t_act)/Dk4 are stable, sub-floor
steps do not corrupt the observable.
"""
from __future__ import annotations
import os, sys, json
import numpy as np, uproot

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HERE, "lib"))
sys.path.insert(0, os.path.join(ROOT, "geom"))
import kink_stats as ks            # noqa
import transport_raytrace as tr    # noqa
RUNS = os.path.join(ROOT, "data", "runs")
W = 79.85e-3; A = tr.A_OF_P[100]; D4 = -0.209; L = 10.0
CELLS = [("rectilinear", 1.0, "c1p0"), ("rectilinear", 0.56, "c0p56"),
         ("rectilinear", 0.3, "c0p3"), ("gyroid", 0.56, "c0p56"),
         ("gyroid", 0.38, "c0p38"), ("gyroid", 0.2, "c0p2")]


def g4(path):
    try:
        f = uproot.open(path)["kinks"]
    except Exception:
        return None, None      # incomplete / still-writing ROOT
    tpla = f["tpla"].array(library="np")
    ang = np.concatenate([f["thetax"].array(library="np"),
                          f["thetay"].array(library="np")])
    return tpla, ang[np.isfinite(ang)]


def kappa4(ang):
    _, k4 = ks.cumulants_in_window(ang, W)
    k4m, lo, hi = ks.bootstrap_kappa4(ang, W)
    return k4m, lo, hi


def main():
    # solid floor at t=4 (100 MeV)
    sol = ks.load_run(os.path.join(RUNS, "solid_E100_t4.root")).angles
    k4_sol4, _, _ = kappa4(sol)

    rows = []
    for name, cell, tag in CELLS:
        rp = os.path.join(RUNS, f"lat_{name}_{tag}_E100.root")
        if not os.path.exists(rp):
            print("MISSING", rp); continue
        tpla, ang = g4(rp)
        if tpla is None:
            print("INCOMPLETE", rp); continue
        var_g4 = float(tpla.var(ddof=1)); tmean = float(tpla.mean())
        k4s, lo, hi = kappa4(ang)
        floor = k4_sol4 * (tmean / 4.0)           # scale floor to lattice <t>
        dk4_g4 = k4s - floor
        # transport at this cell @100
        r = tr.transport_trace(name, cell, 0.40, A, n_proton=60000,
                               steps_per_cell=20, rng=np.random.default_rng(9))
        var_ta = float(np.var(r["t_actual"], ddof=1))
        var_str = float(np.var(r["t_straight"], ddof=1))
        dk4_ta = 3 * A**2 * var_ta * (1 + D4)
        dk4_li = 3 * A**2 * var_str * (1 + D4)
        rows.append(dict(name=name, cell=cell, tmean=tmean,
                         var_g4=var_g4, var_ta=var_ta, var_str=var_str,
                         dvar=var_g4 / var_ta - 1,
                         ratio_g4=var_g4 / var_str, ratio_ta=var_ta / var_str,
                         dk4_g4=dk4_g4, dk4_g4_lo=lo - floor, dk4_g4_hi=hi - floor,
                         dk4_ta=dk4_ta, dk4_li=dk4_li,
                         ta_in_ci=(lo - floor) <= dk4_ta <= (hi - floor)))

    # step-regime check (1b): rect 0.56 res20 vs res40
    step = {}
    for tag in ("res20", "res40"):
        rp = os.path.join(RUNS, f"lat_rectilinear_c0p56_{tag}_E100.root")
        if os.path.exists(rp):
            tpla, ang = g4(rp)
            if tpla is None:
                continue
            k4s, lo, hi = kappa4(ang)
            floor = k4_sol4 * (float(tpla.mean()) / 4.0)
            step[tag] = dict(var=float(tpla.var(ddof=1)), dk4=k4s - floor,
                             lo=lo - floor, hi=hi - floor,
                             voxel_um={"res20": 28, "res40": 14}[tag])

    out = dict(momentum=100, W_mrad=W*1e3, solid_floor_t4=k4_sol4, rows=rows, step=step)
    json.dump(out, open(os.path.join(ROOT, "data", "analysis",
                                     "confirm_100mev.json"), "w"), indent=1)

    print("=== Var(t_actual): Geant4 vs transport @100 MeV ===")
    for r in rows:
        print(f"  {r['name']:11s} c{r['cell']}: f={r['tmean']/L:.3f} "
              f"Var_G4={r['var_g4']:.2f} Var_TA={r['var_ta']:.2f} ({r['dvar']*100:+.1f}%) "
              f"ratio_G4={r['ratio_g4']:.3f} ratio_TA={r['ratio_ta']:.3f}  "
              f"Dk4: G4={r['dk4_g4']:.2e} TA={r['dk4_ta']:.2e} LI={r['dk4_li']:.2e} "
              f"TAinCI={r['ta_in_ci']}")
    if step:
        print("\n=== step-regime (rect 0.56: res20 28um vs res40 14um) ===")
        for tag, s in step.items():
            print(f"  {tag} (voxel {s['voxel_um']}um): Var={s['var']:.2f} "
                  f"Dk4={s['dk4']:.3e} CI[{s['lo']:.2e},{s['hi']:.2e}]")
        if "res20" in step and "res40" in step:
            dv = step["res40"]["var"]/step["res20"]["var"]-1
            ov = not (step["res40"]["hi"] < step["res20"]["lo"] or
                      step["res20"]["hi"] < step["res40"]["lo"])
            print(f"  ΔVar(res40 vs res20)={dv*100:+.1f}%  Δκ4 CIs overlap: {ov}  "
                  f"-> {'STABLE' if abs(dv)<0.03 and ov else 'UNSTABLE (flag)'}")
    return out


if __name__ == "__main__":
    main()
