#!/usr/bin/env python3
"""run_m1.py -- M1 reconstruction-impact telescope campaign driver.

Frozen matrix in configs/M1.yaml. Runs the 6-plane telescope (mcs_sim /mcs/det/telescope 1)
over {rect, gyroid, slab, empty(target-out)} x {200,1000 MeV} x {massless, Si} x {pencil,
divergence}, writing the "telescope" ntuple to data/runs_m1/<tag>.root. Resumable (skips a
tag whose .root already has >= n_events). Seeds + git_commit recorded in each meta.

Usage: python sim/run_m1.py [--nevt 5e6] [--threads 8] [--list]
"""
from __future__ import annotations
import argparse, hashlib, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BIN = os.path.join(HERE, "build", "mcs_sim")
RUNS = os.path.join(ROOT, "data", "runs_m1")
VOX = os.path.join(ROOT, "data", "geom_stats", "voxel")
GENV = os.path.join(ROOT, "geom", "make_campaign_voxels.py")
RUNLOG = os.path.join(ROOT, "RUNLOG.csv")
CELL, DEPTH, INFILL, SLAB_T = 2.5, 10.0, 0.30, 3.0
ENERGIES = (200, 1000)

# (targets, planeSi, beam) waves -- configs/M1.yaml
WAVES = [
    (("rectilinear", "gyroid", "slab", "empty"), False, "pencil"),
    (("rectilinear", "slab", "empty"), True, "pencil"),
    (("rectilinear", "slab"), False, "divergence"),
]


def seed_for(tag, base=4001):
    return (base + int(hashlib.sha256(tag.encode()).hexdigest()[:8], 16)) % (2 ** 31 - 1)


def voxel_path(topo):
    return os.path.join(VOX, f"{topo}_f{int(round(INFILL*100)):02d}_c{CELL:g}_camp_vox.raw")


def ensure_voxel(topo):
    if not os.path.exists(voxel_path(topo)):
        subprocess.run([sys.executable, GENV, topo, str(INFILL)], check=True,
                       stdout=subprocess.DEVNULL)


def jobs():
    out = []
    for targets, si, beam in WAVES:
        for tgt in targets:
            for E in ENERGIES:
                planes = "si" if si else "massless"
                out.append((tgt, si, beam, E,
                            f"m1_{tgt}_{planes}_{beam}_E{E}"))
    return out


def macro(tgt, si, beam, E, out, seed, nevt):
    L = []
    if tgt in ("rectilinear", "gyroid"):
        L += [f"/mcs/det/geom voxel", f"/mcs/det/voxel {voxel_path(tgt)}",
              f"/mcs/det/topology {tgt}", f"/mcs/det/infill {INFILL}"]
    elif tgt == "slab":
        L += ["/mcs/det/geom slab", "/mcs/det/material PLA",
              f"/mcs/det/thickness {SLAB_T} mm"]
    else:  # empty = target-out (vacuum slab)
        L += ["/mcs/det/geom slab", "/mcs/det/material G4_Galactic",
              f"/mcs/det/thickness {SLAB_T} mm"]
    L += ["/mcs/det/telescope 1", f"/mcs/det/planeSi {1 if si else 0}",
          "/mcs/det/maxStep 0.1 mm", "/run/setCut 0.05 mm", "/run/initialize",
          "/gun/particle proton", f"/gun/energy {E} MeV", "/mcs/gun/zStart -350 mm"]
    if beam == "divergence":
        L += ["/mcs/gun/divergence 1e-3", "/mcs/gun/spotXY 1.0 mm"]
    else:
        L += [f"/mcs/gun/spotXY {CELL if tgt in ('rectilinear','gyroid') else 2.5} mm"]
    L += [f"/random/setSeeds {seed} 0", f"/analysis/setFileName {out}",
          f"/run/printProgress {max(1, nevt//4)}", f"/run/beamOn {nevt}"]
    return "\n".join(L) + "\n"


def done(tag, nevt):
    import json
    meta = os.path.join(RUNS, tag + ".root.meta.json")
    root = os.path.join(RUNS, tag + ".root")
    if not (os.path.exists(meta) and os.path.exists(root)):
        return False
    try:
        return json.load(open(meta)).get("n_events", 0) >= nevt
    except Exception:
        return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nevt", default="5e6")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()
    nevt = int(float(args.nevt))
    os.makedirs(RUNS, exist_ok=True)
    try:
        os.environ.setdefault("MCS_GIT_COMMIT", subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT).decode().strip())
    except Exception:
        pass
    commit = os.environ.get("MCS_GIT_COMMIT", "")[:7]

    todo = jobs()
    if args.list:
        for tgt, si, beam, E, tag in todo:
            print(f"{'DONE' if done(tag, nevt) else 'todo'}  {tag}")
        print(f"\n{len(todo)} configs @ {nevt:.0e}")
        return 0

    for tgt, si, beam, E, tag in todo:
        if done(tag, nevt):
            print(f"  skip {tag} (done)"); continue
        if tgt in ("rectilinear", "gyroid"):
            ensure_voxel(tgt)
        sd = seed_for(tag)
        out = os.path.join(RUNS, tag)
        mac = os.path.join(RUNS, f"_{tag}.mac")
        open(mac, "w").write(macro(tgt, si, beam, E, out, sd, nevt))
        t0 = time.time()
        rc = subprocess.run([BIN, mac, str(args.threads)],
                            stdout=open(out + ".run.log", "w"),
                            stderr=subprocess.STDOUT).returncode
        wall = time.time() - t0
        os.remove(mac)
        with open(RUNLOG, "a") as f:
            f.write(f"M1_impact,configs/M1.yaml,{tag},{sd},{nevt},{wall:.1f},"
                    f"{args.threads},{rc},data/runs_m1/{tag}.root,{commit},M1 telescope\n")
        print(f"  {tag:42s} {'OK' if rc==0 else f'RC={rc}'} wall={wall:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
