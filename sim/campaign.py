#!/usr/bin/env python3
"""campaign.py -- S6 resumable campaign driver.

Matrix: 5 topologies x 4 infills x 3 energies = 60 lattice configs, plus baselines
(empty-frame + solid-control + equal-areal-density solid slab per energy). Builds
voxel geometry on demand (geom/make_campaign_voxels.py), generates a macro, runs
mcs_sim (locked physics), and writes data/runs/<tag>.root + provenance. RESUMABLE:
skips any config whose .root already has >= the requested n_events (from the meta
sidecar). Per-run line appended to data/runs/CAMPAIGN_LOG.tsv.

Usage:
  python sim/campaign.py --nevt 1e6 [--energies 200 500 1000] [--threads 6]
                         [--lattice-only | --baselines-only] [--dry 1e4]
                         [--only <topo>:<infill>:<E>] [--list]
Tiered budget: pass --nevt per run; escalate specific configs by re-running with a
larger --nevt (the driver re-runs only configs below the new target).
"""
from __future__ import annotations
import argparse, hashlib, json, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BIN = os.path.join(HERE, "build", "mcs_sim")
RUNS = os.path.join(ROOT, "data", "runs")
VOX = os.path.join(ROOT, "data", "geom_stats", "voxel")
GENV = os.path.join(ROOT, "geom", "make_campaign_voxels.py")
LOG = os.path.join(RUNS, "CAMPAIGN_LOG.tsv")

TOPOS = ("rectilinear", "gyroid", "schwarzp", "diamond", "voronoi")
INFILLS = (0.20, 0.30, 0.40, 0.50)
ENERGIES = (200, 500, 1000)
CELL, DEPTH = 2.5, 10.0
MAXSTEP, CUT = 0.1, 0.05
# Areal-density-matched solid controls: a lattice at infill f has mean traversed
# material f*L, so its matched homogeneous baseline is a t = f*L solid slab. The
# kappa4 deconvolution Delta_kappa4 = kappa4(struct) - kappa4(solid@f*L) requires
# the f-MATCHED thickness (a single t=4mm baseline gives the wrong, even negative,
# Delta_kappa4 for f != 0.40 -- caught in Phase-0c calibration). One control per
# infill: f={0.2,0.3,0.4,0.5} -> t={2,3,4,5} mm.
SOLID_THICKS = sorted({round(f * DEPTH) for f in INFILLS})   # {2,3,4,5} mm


def seed_for(tag, base=12345):
    """Deterministic per-config master RNG seed (emitted via /random/setSeeds so
    every run is reproducible and the seed is recorded in the meta sidecar)."""
    return (base + int(hashlib.sha256(tag.encode()).hexdigest()[:8], 16)) % (2 ** 31 - 1)


def tag_lattice(topo, infill, E):
    return f"camp_{topo}_f{int(round(infill*100)):02d}_E{E}"


def voxel_path(topo, infill):
    return os.path.join(VOX, f"{topo}_f{int(round(infill*100)):02d}_c{CELL:g}_camp_vox.raw")


def done(tag, nevt):
    meta = os.path.join(RUNS, tag + ".root.meta.json")
    root = os.path.join(RUNS, tag + ".root")
    if not (os.path.exists(meta) and os.path.exists(root)):
        return False
    try:
        return json.load(open(meta)).get("n_events", 0) >= nevt
    except Exception:
        return False


def ensure_voxel(topo, infill):
    if not os.path.exists(voxel_path(topo, infill)):
        subprocess.run([sys.executable, GENV, topo, str(infill)], check=True,
                       stdout=subprocess.DEVNULL)


def write_macro(kind, **kw):
    L = []
    if kind == "lattice":
        L += [f"/mcs/det/geom voxel", f"/mcs/det/voxel {kw['vox']}",
              f"/mcs/det/topology {kw['topo']}", f"/mcs/det/infill {kw['infill']}"]
    else:  # slab (solid control / areal) or empty
        mat = "G4_Galactic" if kind == "empty" else "PLA"
        L += [f"/mcs/det/geom slab", f"/mcs/det/material {mat}",
              f"/mcs/det/thickness {kw['thick']} mm"]
    L += [f"/mcs/det/maxStep {MAXSTEP} mm", f"/run/setCut {CUT} mm",
          "/run/initialize", "/gun/particle proton",
          f"/gun/energy {kw['E']} MeV"]
    if kind == "lattice":
        L.append(f"/mcs/gun/spotXY {CELL} mm")
    L += [f"/random/setSeeds {kw['seed']} 0",
          f"/analysis/setFileName {kw['out']}",
          f"/run/printProgress {max(1, kw['nevt']//4)}",
          f"/run/beamOn {kw['nevt']}"]
    return "\n".join(L) + "\n"


def run_one(tag, mac_text, threads, log_extra):
    out = os.path.join(RUNS, tag)
    mac = os.path.join(RUNS, f"_{tag}.mac")
    open(mac, "w").write(mac_text.replace("__OUT__", out))
    t0 = time.time()
    rc = subprocess.run([BIN, mac, str(threads)],
                        stdout=open(out + ".run.log", "w"),
                        stderr=subprocess.STDOUT).returncode
    wall = time.time() - t0
    os.remove(mac)
    stuck = 0
    try:
        stuck = open(out + ".run.log").read().count("Stuck Track")
    except Exception:
        pass
    with open(LOG, "a") as f:
        f.write(f"{tag}\t{rc}\t{wall:.1f}\t{stuck}\t{log_extra}\n")
    return rc, wall, stuck


def jobs(args):
    out = []
    if not args.baselines_only:
        for topo in TOPOS:
            for infill in INFILLS:
                for E in args.energies:
                    out.append(("lattice", topo, infill, E))
    if not args.lattice_only:
        for E in args.energies:
            out.append(("empty", None, None, E))
            for t in SOLID_THICKS:                 # areal-matched controls {2,3,4,5} mm
                out.append(("solid", None, t, E))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nevt", default="1e6")
    ap.add_argument("--energies", nargs="+", type=int, default=list(ENERGIES))
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--lattice-only", action="store_true")
    ap.add_argument("--baselines-only", action="store_true")
    ap.add_argument("--only", default=None, help="topo:infill:E")
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()
    nevt = int(float(args.nevt))
    os.makedirs(RUNS, exist_ok=True)
    # Record the code version in every run's meta sidecar (read by RunAction).
    try:
        os.environ.setdefault("MCS_GIT_COMMIT", subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT).decode().strip())
    except Exception:
        pass

    todo = jobs(args)
    if args.only:
        t, f, e = args.only.split(":")
        todo = [j for j in todo if j[0] == "lattice" and j[1] == t
                and abs(j[2] - float(f)) < 1e-6 and j[3] == int(e)]

    planned, skipped, ran = [], 0, 0
    for kind, topo, infill, E in todo:
        if kind == "lattice":
            tag = tag_lattice(topo, infill, E)
        elif kind == "empty":
            tag = f"empty_E{E}"
        else:
            tag = f"solid_E{E}_t{infill}"
        planned.append((kind, topo, infill, E, tag))

    if args.list:
        for kind, topo, infill, E, tag in planned:
            print(f"{'DONE' if done(tag, nevt) else 'todo'}  {tag}")
        print(f"\n{len(planned)} configs; "
              f"{sum(done(t[4], nevt) for t in planned)} already done at >= {nevt:.0e}")
        return 0

    t_start = time.time()
    for kind, topo, infill, E, tag in planned:
        if done(tag, nevt):
            skipped += 1; continue
        sd = seed_for(tag)
        if kind == "lattice":
            ensure_voxel(topo, infill)
            mac = write_macro("lattice", vox=voxel_path(topo, infill), topo=topo,
                              infill=infill, E=E, out="__OUT__", nevt=nevt, seed=sd)
            extra = f"{topo}\t{infill}\t{E}\tnevt={nevt}\tseed={sd}"
        elif kind == "empty":
            mac = write_macro("empty", thick=16, E=E, out="__OUT__", nevt=nevt, seed=sd)
            extra = f"empty\t-\t{E}\tnevt={nevt}\tseed={sd}"
        else:
            mac = write_macro("solid", thick=infill, E=E, out="__OUT__", nevt=nevt, seed=sd)
            extra = f"solid\t-\t{E}\tnevt={nevt}\tseed={sd}"
        rc, wall, stuck = run_one(tag, mac, args.threads, extra)
        ran += 1
        flag = "OK" if rc == 0 else f"RC={rc}"
        print(f"  {tag:32s} {flag:5s} wall={wall:6.1f}s stuck={stuck}")
    print(f"\nran {ran}, skipped {skipped} (already done), "
          f"total wall {(time.time()-t_start):.0f}s. Log: {LOG}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
