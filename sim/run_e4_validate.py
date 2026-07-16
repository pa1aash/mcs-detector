#!/usr/bin/env python3
"""run_e4_validate.py -- S7 / E4 Geant4 spot-check of the DOMINANT print artifact.

Runs the perturbed rectilinear-f40 fields (dimensional tolerance under/over-extrusion)
under LOCKED physics (option4, no env override), locked W(E), and the locked all-order
kappa_M(tpla) floor, then compares the MEASURED Delta_kappa4 shift to the ray-traced
prediction in e4_print_realism.json. The nominal is the on-disk locked campaign run
camp_rectilinear_f40_E200. Confirms the (S5-validated) ray-tracer's artifact shift is
physical.

Usage:  python sim/run_e4_validate.py [--nevt 3e6] [--artifacts tolerance_under tolerance_over]
"""
from __future__ import annotations
import argparse, json, os, subprocess, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BIN = os.path.join(HERE, "build", "mcs_sim")
VOX_E4 = os.path.join(ROOT, "data", "geom_stats", "voxel_e4")
RUNS = os.path.join(ROOT, "data", "runs")
OUTDIR = os.path.join(ROOT, "data", "runs_e4")
sys.path.insert(0, os.path.join(ROOT, "analysis"))
sys.path.insert(0, os.path.join(ROOT, "analysis", "lib"))
import kink_stats as ks   # noqa
import e2_analysis as e2  # noqa

E, F, CELL = 200, 0.40, 2.5
MAXSTEP, CUT = 0.1, 0.05


def macro(stem, out, nevt):
    return "\n".join([
        "/mcs/det/geom voxel", f"/mcs/det/voxel {stem}.raw",
        "/mcs/det/topology rectilinear", f"/mcs/det/infill {F}",
        f"/mcs/det/maxStep {MAXSTEP} mm", f"/run/setCut {CUT} mm",
        "/run/initialize", "/gun/particle proton", f"/gun/energy {E} MeV",
        f"/mcs/gun/spotXY {CELL} mm", f"/analysis/setFileName {out}",
        f"/run/printProgress {max(1,nevt//4)}", f"/run/beamOn {nevt}"]) + "\n"


def dk4_of(root_path):
    import uproot
    W = e2.W[E]; kM = e2.build_floor(E)
    f = uproot.open(root_path)["kinks"]
    ang = np.concatenate([f["thetax"].array(library="np"), f["thetay"].array(library="np")])
    tpla = f["tpla"].array(library="np")
    k2, _ = ks.cumulants_in_window(ang, W)
    k4s, se = e2._k4_pt_se(ang, W, np.random.default_rng(9))
    dk4 = k4s - float(np.mean(kM(tpla)))
    return dict(dk4=float(dk4), dk4_se=float(se), k2=float(k2),
                gamma2=float(dk4 / k2 ** 2), tpla_mean=float(np.mean(tpla)),
                f_built=float(np.mean(tpla)) / e2.L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nevt", default="3e6")
    ap.add_argument("--artifacts", nargs="+", default=["tolerance_under", "tolerance_over"])
    ap.add_argument("--threads", type=int, default=6)
    args = ap.parse_args()
    nevt = int(float(args.nevt))
    os.makedirs(OUTDIR, exist_ok=True)

    rt_pred = json.load(open(os.path.join(ROOT, "data", "analysis", "e4_print_realism.json")))
    nom_pred = rt_pred["configs"]["rectilinear_f40"]

    nominal = dk4_of(os.path.join(RUNS, "camp_rectilinear_f40_E200.root"))
    out = {"E": E, "nominal_g4": nominal, "artifacts": {}}
    print(f"nominal (locked camp rect f40 @200): dk4={nominal['dk4']:.3e} g2={nominal['gamma2']:.2f}")
    for aname in args.artifacts:
        stem = os.path.join(VOX_E4, f"rectilinear_f40_c{CELL:g}_{aname}_vox")
        if not os.path.exists(stem + ".raw"):
            print(f"  [skip {aname}] no perturbed field"); continue
        tag = f"e4_rect_f40_{aname}_E{E}"
        outp = os.path.join(OUTDIR, tag)
        mac = os.path.join(OUTDIR, f"_{tag}.mac")
        open(mac, "w").write(macro(stem, outp, nevt))
        t0 = time.time()
        rc = subprocess.run([BIN, mac, str(args.threads)],
                            stdout=open(outp + ".run.log", "w"),
                            stderr=subprocess.STDOUT).returncode
        os.remove(mac)
        g4 = dk4_of(outp + ".root")
        meas_shift = g4["dk4"] / nominal["dk4"] - 1.0
        pred_shift = nom_pred["perturbed"][aname]["dDelta_kappa4_pct"] / 100.0
        out["artifacts"][aname] = dict(
            g4=g4, meas_dk4_shift_pct=meas_shift * 100, raytrace_pred_pct=pred_shift * 100,
            wall_s=time.time() - t0, rc=rc)
        print(f"  {aname:16s}: G4 dk4={g4['dk4']:.3e}  shift MEAS {meas_shift*100:+.1f}% "
              f"vs RAYTRACE {pred_shift*100:+.1f}%  (wall {time.time()-t0:.0f}s)")
    json.dump(out, open(os.path.join(ROOT, "data", "analysis", "e4_validate.json"), "w"), indent=1)
    print("wrote data/analysis/e4_validate.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
