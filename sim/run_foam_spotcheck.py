#!/usr/bin/env python3
"""run_foam_spotcheck.py -- S6 Stage-2 Phase-1: full Geant4 gamma2 at the 0.2 mm foam-scale
cell (rect+gyroid, 200 MeV) at two coarse voxel resolutions (step-regime guard), then
gamma2 = Delta_kappa4 / kappa2^2 with the LOCKED W(200) window and the matched t=f*L solid
control (all-order floor), compared to the transport-tool prediction (~2.2 rect / ~1.2
gyroid). Reuses solid_E200_t4 (f*L=4 mm) as the matched control.

Usage:  python sim/run_foam_spotcheck.py [--nevt 1e7]
"""
from __future__ import annotations
import argparse, json, os, subprocess, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BIN = os.path.join(HERE, "build", "mcs_sim")
RUNS = os.path.join(ROOT, "data", "runs")
VOX = os.path.join(ROOT, "data", "geom_stats", "voxel")
sys.path.insert(0, os.path.join(ROOT, "analysis"))
sys.path.insert(0, os.path.join(ROOT, "analysis", "lib"))
sys.path.insert(0, os.path.join(ROOT, "geom"))
import kink_stats as ks       # noqa
import e2_analysis as e2      # noqa
import transport_raytrace as tr  # noqa

CELL, INFILL, E, F = 0.2, 0.40, 200, 0.40
TOPOS, RESES = ("rectilinear", "gyroid"), (4, 8)
MAXSTEP, CUT = 0.1, 0.05


def voxstem(name, res):
    return os.path.join(VOX, f"{name}_f{int(INFILL*100):02d}_c{CELL:g}_res{res}_foam_vox")


def macro(name, res, out, nevt):
    return "\n".join([
        "/mcs/det/geom voxel", f"/mcs/det/voxel {voxstem(name,res)}.raw",
        f"/mcs/det/topology {name}", f"/mcs/det/infill {INFILL}",
        f"/mcs/det/maxStep {MAXSTEP} mm", f"/run/setCut {CUT} mm",
        "/run/initialize", "/gun/particle proton", f"/gun/energy {E} MeV",
        f"/mcs/gun/spotXY {CELL} mm", f"/analysis/setFileName {out}",
        f"/run/printProgress {max(1,nevt//4)}", f"/run/beamOn {nevt}"]) + "\n"


def run_one(name, res, nevt, threads=6):
    tag = f"foam_{name}_res{res}_E{E}"
    out = os.path.join(RUNS, tag)
    mac = os.path.join(RUNS, f"_{tag}.mac")
    open(mac, "w").write(macro(name, res, out, nevt))
    t0 = time.time()
    rc = subprocess.run([BIN, mac, str(threads)], stdout=open(out + ".run.log", "w"),
                        stderr=subprocess.STDOUT).returncode
    try:
        os.remove(mac)
    except FileNotFoundError:
        pass
    stuck = open(out + ".run.log").read().count("Stuck Track")
    return tag, rc, time.time() - t0, stuck


def gamma2_of(tag):
    """gamma2 = Delta_kappa4 / kappa2^2, all-order floor subtraction, W(200)."""
    import uproot
    W = e2.W[E]; aeff = e2.a_eff_of_E(E); kM = e2.build_floor(E)
    f = uproot.open(os.path.join(RUNS, tag + ".root"))["kinks"]
    ang = np.concatenate([f["thetax"].array(library="np"), f["thetay"].array(library="np")])
    tpla = f["tpla"].array(library="np")
    k2, _ = ks.cumulants_in_window(ang, W)
    k4s, se = e2._k4_pt_se(ang, W, np.random.default_rng(5))
    dk4 = k4s - float(np.mean(kM(tpla)))
    g2 = dk4 / k2 ** 2
    g2_ci = 1.96 * se / k2 ** 2
    return dict(gamma2=float(g2), gamma2_ci=float(g2_ci), k2=float(k2), dk4=float(dk4),
                f_built=float(np.mean(tpla)) / e2.L, tpla_mean=float(np.mean(tpla)),
                n=int(ang.size))


def gamma2_straight_raw(name, res):
    """Straight-chord gamma2 of the SAME voxel field Geant4 ran (geometry-matched 'no-
    averaging' reference; f-mismatch-immune). gamma2_G4/gamma2_straight = Var(t_act)/
    Var(t_str) -> 1 at foam scale (wandering too small to average) is the claim to test."""
    stem = voxstem(name, res)
    m = open(stem + ".raw.meta").read().split()
    Nx, Ny, Nz, voxel = int(m[0]), int(m[1]), int(m[2]), float(m[3])
    chi = np.fromfile(stem + ".raw", dtype=np.uint8).reshape(Nx, Ny, Nz)
    t = chi.reshape(Nx * Ny, Nz).sum(axis=1) * voxel
    return float(3 * np.var(t, ddof=1) / t.mean() ** 2)


def transport_gamma2(name):
    a = tr.A_OF_P[E]
    r = tr.transport_trace(name, CELL, F, a, n_proton=40000, steps_per_cell=20,
                           rng=np.random.default_rng(7))
    ta = r["t_actual"]
    return float(3 * np.var(ta, ddof=1) / ta.mean() ** 2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nevt", default="1e7")
    args = ap.parse_args()
    nevt = int(float(args.nevt))
    res_out = {}
    for name in TOPOS:
        res_out[name] = {"transport_gamma2": transport_gamma2(name), "by_res": {}}
        for res in RESES:
            tag, rc, wall, stuck = run_one(name, res, nevt)
            g = gamma2_of(tag)
            g["gamma2_straight"] = gamma2_straight_raw(name, res)   # geometry-matched ref
            g["ratio_g4_straight"] = g["gamma2"] / g["gamma2_straight"]
            g.update(rc=rc, wall=round(wall, 1), stuck=stuck, tag=tag)
            res_out[name]["by_res"][str(res)] = g
            print(f"{name:11s} res{res}: g2_G4={g['gamma2']:.3f}+/-{g['gamma2_ci']:.3f} "
                  f"g2_straight(raw)={g['gamma2_straight']:.3f}  ratio_G4/straight="
                  f"{g['ratio_g4_straight']:.3f}  f={g['f_built']:.3f} stuck={stuck} "
                  f"wall={wall:.0f}s")
        g4 = res_out[name]["by_res"]
        # CONFIRMED if ratio g2_G4/g2_straight ~ 1 (wandering does NOT average at foam scale)
        # on BOTH resolutions; voxel-stable if the ratio agrees across res.
        r4, r8 = g4["4"]["ratio_g4_straight"], g4["8"]["ratio_g4_straight"]
        res_out[name]["ratio_voxel_stable"] = bool(abs(r4 - r8) < 0.10)
        print(f"  -> ratio_G4/straight res4={r4:.3f} res8={r8:.3f}  voxel-stable="
              f"{res_out[name]['ratio_voxel_stable']}  (transport analytic g2="
              f"{res_out[name]['transport_gamma2']:.2f})")
    out = dict(cell_mm=CELL, energy=E, infill=INFILL, nevt=nevt, by_topo=res_out)
    json.dump(out, open(os.path.join(ROOT, "data", "analysis", "foam_spotcheck.json"), "w"),
              indent=1)
    print("wrote data/analysis/foam_spotcheck.json")


if __name__ == "__main__":
    main()
