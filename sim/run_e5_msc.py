#!/usr/bin/env python3
"""run_e5_msc.py -- S7 / E5 MSC-model systematic driver.

Re-runs a SUBSET of the 200+500 MeV campaign at the printable 2.5 mm cell under an
ALTERNATIVE multiple-scattering model, plus the matched solid controls (which also
re-derive the Highland `a` and the kappa_M(t) floor for the E5(b) boundary). Outputs
go to data/runs_e5/<variant>/ with the SAME tag names the locked pipeline uses, so
analysis/e2_analysis.py --runs-dir data/runs_e5/<variant> re-fits the collapse with
no other change. Geometry (voxel fields) is REUSED from the locked campaign -> the
ray-traced N_eff is bit-identical, isolating the MSC effect on Delta_kappa4.

The alternative MSC is selected by env vars read in PhysicsList (see PhysicsList.cc):
  opt3   : MCS_EM_PHYSICS=option3  (UrbanMsc-based EM, vs the locked option4
           WentzelVI + single-Coulomb split)
  thlim<X>: MCS_MSC_THETALIMIT_DEG=X (WentzelVI MSC/single-Coulomb boundary)

Usage:
  python sim/run_e5_msc.py --variant opt3 --energies 500 [--threads 6] [--list]
"""
from __future__ import annotations
import argparse, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BIN = os.path.join(HERE, "build", "mcs_sim")
VOX = os.path.join(ROOT, "data", "geom_stats", "voxel")
CELL, DEPTH, MAXSTEP, CUT = 2.5, 10.0, 0.1, 0.05

# Variant -> env vars set for mcs_sim.
VARIANTS = {
    "opt3":     {"MCS_EM_PHYSICS": "option3"},
    "thlim2":   {"MCS_MSC_THETALIMIT_DEG": "2.0"},   # tighten WentzelVI->single split
}

# Collapse subset: anchors (N_eff~2) + mid + high-N_eff voronoi, spanning N_eff 2->13.
# (topo, infill, nevt) -- nevt matches the locked campaign for rect/schwarzp/gyroid;
# voronoi at 1e7 (vs locked 3e7) still resolves the deconvolved Delta_kappa4 to <30%.
LATTICE_PROFILES = {
    # full: the original subset; voronoi escalated to 3e7 (the alt-MSC signal is
    # ~40-70% smaller than locked, so the high-N_eff voronoi needs the locked 3e7
    # stats to resolve the deconvolved Delta_kappa4 to <30%).
    "full": [
        ("rectilinear", 0.40, 3_000_000),   # N_eff ~ 2.07 (anchor)
        ("schwarzp",    0.40, 3_000_000),   # N_eff ~ 2.18 (anchor)
        ("gyroid",      0.40, 6_000_000),   # N_eff ~ 3.88 (mid)
        ("gyroid",      0.20, 6_000_000),   # N_eff ~ 5.55 (mid)
        ("voronoi",     0.40, 30_000_000),  # N_eff ~ 9.32 (high)
        ("voronoi",     0.50, 30_000_000),  # N_eff ~ 12.62 (high)
    ],
    # esc_voronoi: only the two high-N_eff voronoi at 3e7 (escalate to recover the
    # lever arm; the periodic/gyroid points are already done at their targets).
    "esc_voronoi": [
        ("voronoi", 0.40, 30_000_000),
        ("voronoi", 0.50, 30_000_000),
    ],
    # cheap: the resolvable periodic+gyroid anchors (no voronoi); for a second energy
    # where only the N_eff 2-5.5 lever arm is needed to confirm the systematic.
    "cheap": [
        ("rectilinear", 0.40, 3_000_000),
        ("schwarzp",    0.40, 3_000_000),
        ("gyroid",      0.40, 6_000_000),
        ("gyroid",      0.20, 6_000_000),
    ],
}
# Solid controls: a_eff (t={2,3,4,5}) + kappa_M(t) floor curvature (t={8,16}).
SOLID_THICKS = (2, 3, 4, 5, 8, 16)
SOLID_NEVT = 10_000_000


def voxel_path(topo, infill):
    return os.path.join(VOX, f"{topo}_f{int(round(infill*100)):02d}_c{CELL:g}_camp_vox.raw")


def macro_lattice(topo, infill, E, out, nevt):
    return "\n".join([
        "/mcs/det/geom voxel", f"/mcs/det/voxel {voxel_path(topo,infill)}",
        f"/mcs/det/topology {topo}", f"/mcs/det/infill {infill}",
        f"/mcs/det/maxStep {MAXSTEP} mm", f"/run/setCut {CUT} mm",
        "/run/initialize", "/gun/particle proton", f"/gun/energy {E} MeV",
        f"/mcs/gun/spotXY {CELL} mm",
        f"/analysis/setFileName {out}",
        f"/run/printProgress {max(1, nevt//4)}", f"/run/beamOn {nevt}"]) + "\n"


def foam_voxstem(name, res):
    return os.path.join(VOX, f"{name}_f40_c0.2_res{res}_foam_vox")


def macro_foam(name, res, E, out, nevt):
    # 0.2 mm foam-scale cell, f=0.40, matches the locked run_foam_spotcheck config.
    return "\n".join([
        "/mcs/det/geom voxel", f"/mcs/det/voxel {foam_voxstem(name,res)}.raw",
        f"/mcs/det/topology {name}", "/mcs/det/infill 0.4",
        f"/mcs/det/maxStep {MAXSTEP} mm", f"/run/setCut {CUT} mm",
        "/run/initialize", "/gun/particle proton", f"/gun/energy {E} MeV",
        "/mcs/gun/spotXY 0.2 mm", f"/analysis/setFileName {out}",
        f"/run/printProgress {max(1,nevt//4)}", f"/run/beamOn {nevt}"]) + "\n"


def macro_solid(thick, E, out, nevt):
    return "\n".join([
        "/mcs/det/geom slab", "/mcs/det/material PLA",
        f"/mcs/det/thickness {thick} mm",
        f"/mcs/det/maxStep {MAXSTEP} mm", f"/run/setCut {CUT} mm",
        "/run/initialize", "/gun/particle proton", f"/gun/energy {E} MeV",
        f"/analysis/setFileName {out}",
        f"/run/printProgress {max(1, nevt//4)}", f"/run/beamOn {nevt}"]) + "\n"


def done(outdir, tag, nevt):
    import json
    meta = os.path.join(outdir, tag + ".root.meta.json")
    root = os.path.join(outdir, tag + ".root")
    if not (os.path.exists(meta) and os.path.exists(root)):
        return False
    try:
        return json.load(open(meta)).get("n_events", 0) >= nevt
    except Exception:
        return False


def run_one(env, outdir, tag, mac_text, threads):
    out = os.path.join(outdir, tag)
    mac = os.path.join(outdir, f"_{tag}.mac")
    open(mac, "w").write(mac_text.replace("__OUT__", out))
    full_env = dict(os.environ); full_env.update(env)
    t0 = time.time()
    rc = subprocess.run([BIN, mac, str(threads)], env=full_env,
                        stdout=open(out + ".run.log", "w"),
                        stderr=subprocess.STDOUT).returncode
    wall = time.time() - t0
    os.remove(mac)
    stuck = 0
    try:
        stuck = open(out + ".run.log").read().count("Stuck Track")
    except Exception:
        pass
    return rc, wall, stuck


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", default="opt3", choices=list(VARIANTS))
    ap.add_argument("--profile", default="full", choices=list(LATTICE_PROFILES))
    ap.add_argument("--energies", nargs="+", type=int, default=[500])
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--solids-only", action="store_true")
    ap.add_argument("--lattice-only", action="store_true")
    ap.add_argument("--foam", action="store_true",
                    help="also run the 0.2 mm foam-scale rect res8 point (E5b boundary)")
    args = ap.parse_args()

    env = VARIANTS[args.variant]
    outdir = os.path.join(ROOT, "data", "runs_e5", args.variant)
    os.makedirs(outdir, exist_ok=True)

    jobs = []  # (kind, tag, macro_text, nevt)
    for E in args.energies:
        if not args.lattice_only:
            for t in SOLID_THICKS:
                tag = f"solid_E{E}_t{t}"
                jobs.append(("solid", tag, macro_solid(t, E, "__OUT__", SOLID_NEVT), SOLID_NEVT))
        if not args.solids_only:
            for topo, infill, nevt in LATTICE_PROFILES[args.profile]:
                tag = f"camp_{topo}_f{int(round(infill*100)):02d}_E{E}"
                jobs.append(("lattice", tag,
                             macro_lattice(topo, infill, E, "__OUT__", nevt), nevt))
        if args.foam and E == 200:
            jobs.append(("foam", f"foam_rectilinear_res8_E{E}",
                         macro_foam("rectilinear", 8, E, "__OUT__", 10_000_000), 10_000_000))

    if args.list:
        for kind, tag, _, nevt in jobs:
            print(f"{'DONE' if done(outdir,tag,nevt) else 'todo'}  {tag:32s} nevt={nevt:.0e}")
        print(f"\n{len(jobs)} configs -> {outdir}  (variant={args.variant} env={env})")
        return 0

    print(f"[E5/{args.variant}] env={env}  out={outdir}  {len(jobs)} configs")
    t_start = time.time()
    for kind, tag, mac, nevt in jobs:
        if done(outdir, tag, nevt):
            print(f"  SKIP {tag} (done)"); continue
        rc, wall, stuck = run_one(env, outdir, tag, mac, args.threads)
        flag = "OK" if rc == 0 else f"RC={rc}"
        print(f"  {tag:32s} {flag:5s} wall={wall:7.1f}s stuck={stuck}", flush=True)
    print(f"\n[E5/{args.variant}] total wall {(time.time()-t_start):.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
