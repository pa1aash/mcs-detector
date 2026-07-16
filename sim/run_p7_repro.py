#!/usr/bin/env python3
"""run_p7_repro.py -- P7 provenance re-run: the full 200+500 MeV headline campaign,
fresh, seeded, at the original tiered statistics, so every figure input regenerates
from tracked, reproducible artifacts.

Tiers (mirroring the original campaign design, data/runs/CAMPAIGN_QC_200.md):
  rectilinear / schwarzp 3e6; gyroid 6e6; diamond 1e6; voronoi 3e7 per config;
  solid controls t in {2,3,4,5,8,16} mm at 1e7; empty frame at 3e6.
Resumable: skips any run whose meta already records >= the target n_events.
Uses sim/campaign.py machinery (tags, deterministic seeds, macro, meta, log).
"""
from __future__ import annotations
import json, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import campaign as cp  # noqa

ENERGIES = (200, 500)
NEVT_LAT = {"diamond": 1_000_000, "rectilinear": 3_000_000, "schwarzp": 3_000_000,
            "gyroid": 6_000_000, "voronoi": 30_000_000}
NEVT_SOLID, NEVT_EMPTY = 10_000_000, 3_000_000
SOLID_THICKS = (2, 3, 4, 5, 8, 16)
THREADS = 10
RUNLOG = os.path.join(ROOT, "RUNLOG.csv")


def log_runlog(tag, seed, nevt, wall, rc, note):
    commit = os.environ.get("MCS_GIT_COMMIT", "working")
    with open(RUNLOG, "a") as f:
        f.write(f"P7_repro,{tag},{tag},{seed},{nevt},{wall:.1f},{THREADS},{rc},"
                f"data/runs/{tag}.root,{commit},{note}\n")


def run(kind, tag, nevt, mac_kw):
    if cp.done(tag, nevt):
        print(f"skip {tag} (done)")
        return 0
    seed = cp.seed_for(tag)
    mac = cp.write_macro(kind, seed=seed, out="__OUT__", nevt=nevt, **mac_kw)
    print(f"run  {tag}  nevt={nevt:.0e} seed={seed}", flush=True)
    rc, wall, stuck = cp.run_one(tag, mac, THREADS, f"P7_repro nevt={nevt}")
    log_runlog(tag, seed, nevt, wall, rc, f"stuck={stuck}")
    if rc != 0:
        print(f"FAIL {tag} rc={rc}", flush=True)
    return rc


def main():
    t0 = time.time()
    for E in ENERGIES:
        # baselines first (floor calibration)
        run("empty", f"empty_E{E}", NEVT_EMPTY, dict(thick=16, E=E))
        for t in SOLID_THICKS:
            run("solid", f"solid_E{E}_t{t}", NEVT_SOLID, dict(thick=t, E=E))
        # lattices, cheap tiers first
        for topo in ("diamond", "rectilinear", "schwarzp", "gyroid", "voronoi"):
            for infill in cp.INFILLS:
                cp.ensure_voxel(topo, infill)
                tag = cp.tag_lattice(topo, infill, E)
                run("lattice", tag, NEVT_LAT[topo],
                    dict(vox=cp.voxel_path(topo, infill), topo=topo, infill=infill, E=E))
    print(f"P7 repro complete in {(time.time()-t0)/3600:.2f} h wall", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
