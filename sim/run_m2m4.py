#!/usr/bin/env python3
"""run_m2m4.py -- M2 (tilt scan) + M4 (carbon) campaign driver.

M2 (data/runs_m2/): rect+gyroid voxel at tilt {0,5,15,30} deg + rect fine {1,2,3} deg,
  200 MeV; PLA solid controls t={2,3,4,5} + empty (floor). Uses /mcs/det/tilt.
M4 (data/runs_m4/): C_amorph (rho 1.7) gyroid+voronoi at f={0.05,0.30}, 200 MeV;
  carbon solid controls t={1,2,3,4} + empty. Uses /mcs/det/material C_amorph.
Resumable; seeds + git_commit in each meta. Usage: python sim/run_m2m4.py [--nevt 3e6].
"""
from __future__ import annotations
import argparse, hashlib, json, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BIN = os.path.join(HERE, "build", "mcs_sim")
VOX = os.path.join(ROOT, "data", "geom_stats", "voxel")
GENV = os.path.join(ROOT, "geom", "make_campaign_voxels.py")
RUNLOG = os.path.join(ROOT, "RUNLOG.csv")
CELL, L, E = 2.5, 10.0, 200


def seed_for(tag, base=6001):
    return (base + int(hashlib.sha256(tag.encode()).hexdigest()[:8], 16)) % (2 ** 31 - 1)


def vox(topo, f):
    return os.path.join(VOX, f"{topo}_f{int(round(f*100)):02d}_c{CELL:g}_camp_vox.raw")


def ensure_vox(topo, f):
    if not os.path.exists(vox(topo, f)):
        subprocess.run([sys.executable, GENV, topo, str(f)], check=True, stdout=subprocess.DEVNULL)


def jobs():
    J = []  # (dir, tag, kind, params)
    # M2 tilt scan. Geant4 only in the block-clean regime: the 7.5 mm (3-cell) block has
    # half-width 3.75 mm, so a straight beam through a target tilted by theta stays inside
    # when the beam drift (L*tan theta) + spot/2 < half-width -> theta <~ 10 deg. The
    # analytic N_eff(theta) rocking curve (m2_tilt.py, infinite lattice) covers 0-30 deg;
    # Geant4 anchors it in the clean regime (15-30 deg = [PREDICTED], needs wider blocks).
    for topo in ("rectilinear", "gyroid"):
        for t in (0, 5, 10):
            J.append(("runs_m2", f"m2_{topo}_tilt{t}", "lattice",
                      dict(topo=topo, f=0.30, mat="PLA", tilt=t)))
    for t in (1, 2, 3):
        J.append(("runs_m2", f"m2_rectilinear_tilt{t}", "lattice",
                  dict(topo="rectilinear", f=0.30, mat="PLA", tilt=t)))
    for t in (2, 3, 4, 5):
        J.append(("runs_m2", f"m2_solid_t{t}", "solid", dict(mat="PLA", thick=t)))
    J.append(("runs_m2", "m2_empty", "empty", dict(thick=16)))
    # M4 carbon
    for topo in ("gyroid", "voronoi"):
        for f in (0.05, 0.30):
            J.append(("runs_m4", f"m4_{topo}_f{int(round(f*100)):02d}", "lattice",
                      dict(topo=topo, f=f, mat="C_amorph", tilt=0)))
    for t in (1, 2, 3, 4):
        J.append(("runs_m4", f"m4_csolid_t{t}", "solid", dict(mat="C_amorph", thick=t)))
    J.append(("runs_m4", "m4_empty", "empty", dict(thick=16)))
    return J


def macro(kind, p, out, seed, nevt):
    L_ = []
    if kind == "lattice":
        L_ += [f"/mcs/det/geom voxel", f"/mcs/det/voxel {vox(p['topo'], p['f'])}",
               f"/mcs/det/material {p['mat']}", f"/mcs/det/topology {p['topo']}",
               f"/mcs/det/infill {p['f']}", f"/mcs/det/tilt {p['tilt']} deg"]
    elif kind == "solid":
        L_ += ["/mcs/det/geom slab", f"/mcs/det/material {p['mat']}",
               f"/mcs/det/thickness {p['thick']} mm"]
    else:  # empty
        L_ += ["/mcs/det/geom slab", "/mcs/det/material G4_Galactic",
               f"/mcs/det/thickness {p['thick']} mm"]
    L_ += ["/mcs/det/maxStep 0.1 mm", "/run/setCut 0.05 mm", "/run/initialize",
           "/gun/particle proton", f"/gun/energy {E} MeV"]
    if kind == "lattice":
        L_.append(f"/mcs/gun/spotXY {CELL} mm")
    L_ += [f"/random/setSeeds {seed} 0", f"/analysis/setFileName {out}",
           f"/run/printProgress {max(1, nevt//4)}", f"/run/beamOn {nevt}"]
    return "\n".join(L_) + "\n"


def done(d, tag, nevt):
    runs = os.path.join(ROOT, "data", d)
    meta, root = os.path.join(runs, tag + ".root.meta.json"), os.path.join(runs, tag + ".root")
    if not (os.path.exists(meta) and os.path.exists(root)):
        return False
    try:
        return json.load(open(meta)).get("n_events", 0) >= nevt
    except Exception:
        return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nevt", default="3e6")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()
    nevt = int(float(args.nevt))
    try:
        os.environ.setdefault("MCS_GIT_COMMIT", subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT).decode().strip())
    except Exception:
        pass
    commit = os.environ.get("MCS_GIT_COMMIT", "")[:7]
    for d, tag, kind, p in jobs():
        runs = os.path.join(ROOT, "data", d)
        os.makedirs(runs, exist_ok=True)
        if args.list:
            print(f"{'DONE' if done(d,tag,nevt) else 'todo'}  {d}/{tag}"); continue
        if done(d, tag, nevt):
            print(f"  skip {tag}"); continue
        if kind == "lattice":
            ensure_vox(p["topo"], p["f"])
        sd = seed_for(tag)
        out = os.path.join(runs, tag)
        mac = os.path.join(runs, f"_{tag}.mac")
        open(mac, "w").write(macro(kind, p, out, sd, nevt))
        t0 = time.time()
        rc = subprocess.run([BIN, mac, str(args.threads)], stdout=open(out + ".run.log", "w"),
                            stderr=subprocess.STDOUT).returncode
        wall = time.time() - t0
        os.remove(mac)
        with open(RUNLOG, "a") as f:
            f.write(f"{d.replace('runs_','').upper()},sim/run_m2m4.py,{tag},{sd},{nevt},"
                    f"{wall:.1f},{args.threads},{rc},data/{d}/{tag}.root,{commit},M2/M4\n")
        print(f"  {tag:28s} {'OK' if rc==0 else f'RC={rc}'} wall={wall:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
