#!/usr/bin/env python3
"""run_m3.py -- M3 campaign driver: third energy (1000 MeV) + wide-N_eff cells.

Two blocks, pre-registered in ANALYSIS_PLAN.md (section M3), frozen configs/M3.yaml.
Ordered so the Definition-of-Done-critical parts complete first (an early stop still
leaves the collapse's third energy + >=2-decade lever arm intact):

  M3-A  third energy @ 1000 MeV, standard 2.5 mm cell -> data/runs/ with NATIVE
        campaign tags (camp_<topo>_f<ff>_E1000, solid_E1000_t<t>, empty_E1000) so
        analysis/e2_analysis.py --energies 1000 reads them with no changes. Lattice
        rect/gyroid/schwarzp x f{0.20,0.30,0.40,0.50}; PLA solids t{2,3,4,5,8,16}
        (the frozen floor set); empty frame. 1e7 primaries (MCS width 8.95 mrad at
        1 GeV is small -> kappa4 needs the stats).

  M3-B  wide-N_eff cells @ f=0.30, rect/gyroid/schwarzp x cell{0.5,1.0,5.0} mm ->
        data/runs_m3/ with m3_<topo>_f30_c<cell>_E<E>. All cells >= c_break at both
        energies (0.5 mm = 500 um >= c_break(200)~291 um >= c_break(1000)~72 um), so
        the straight-chord Var(t) model is valid (no wandering-regime contamination).
        200 MeV first (3e6, big signal; reuses committed a_eff(200) + the M2 PLA-solid
        floor), then 1000 MeV (1e7; reuses the M3-A 1000 solids/floor).

Physics + tags are byte-identical to campaign.py for the M3-A block (same macro, same
seed_for base=12345), so those runs are indistinguishable from a campaign.py 1000 MeV
run. Resumable (skips any tag whose .root already has >= its target n_events). Seeds +
git_commit recorded in every meta. Usage:
  python sim/run_m3.py [--only A|B200|B1000] [--list] [--threads 8]
"""
from __future__ import annotations
import argparse, hashlib, json, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BIN = os.path.join(HERE, "build", "mcs_sim")
VOX = os.path.join(ROOT, "data", "geom_stats", "voxel")
GENV = os.path.join(ROOT, "geom", "make_campaign_voxels.py")
RUNLOG = os.path.join(ROOT, "RUNLOG.csv")
L, INFILL = 10.0, 0.30
TOPOS = ("rectilinear", "gyroid", "schwarzp")
FILLS = (0.20, 0.30, 0.40, 0.50)
SOLID_THICKS = (2, 3, 4, 5, 8, 16)            # frozen floor set (ANALYSIS_PLAN sec 0)
WIDE_CELLS = (0.5, 1.0, 5.0)
NEVT_1000, NEVT_WIDE200 = int(1e7), int(3e6)  # per-config targets


def seed_native(tag, base=12345):              # matches campaign.py::seed_for
    return (base + int(hashlib.sha256(tag.encode()).hexdigest()[:8], 16)) % (2 ** 31 - 1)


def seed_wide(tag, base=8001):
    return (base + int(hashlib.sha256(tag.encode()).hexdigest()[:8], 16)) % (2 ** 31 - 1)


def vox_path(topo, f, cell):
    return os.path.join(VOX, f"{topo}_f{int(round(f*100)):02d}_c{cell:g}_camp_vox.raw")


def ensure_vox(topo, f, cell):
    if not os.path.exists(vox_path(topo, f, cell)):
        cmd = [sys.executable, GENV, topo, str(f)]
        if abs(cell - 2.5) > 1e-9:
            cmd.append(str(cell))
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)


def jobs():
    """(runs_dir, tag, kind, params, E, nevt) tuples in DoD-priority order."""
    J = []
    # --- M3-A: third energy @ 1000 MeV, standard 2.5 mm cell, native campaign tags ---
    for topo in TOPOS:
        for f in FILLS:
            tag = f"camp_{topo}_f{int(round(f*100)):02d}_E1000"
            J.append(("runs", tag, "lattice",
                      dict(topo=topo, f=f, cell=2.5, spot=2.5), 1000, NEVT_1000))
    for t in SOLID_THICKS:
        J.append(("runs", f"solid_E1000_t{t}", "solid", dict(thick=t), 1000, NEVT_1000))
    J.append(("runs", "empty_E1000", "empty", dict(thick=16), 1000, NEVT_1000))
    # --- M3-B: wide-N_eff cells @ f=0.30 ---
    for E, nevt in ((200, NEVT_WIDE200), (1000, NEVT_1000)):
        for topo in TOPOS:
            for cell in WIDE_CELLS:
                # encode cell in dot-free micrometres: G4 /analysis/setFileName reads a
                # "." in the output name as a file extension (a dotted tag segfaults on write).
                tag = f"m3_{topo}_f30_c{int(round(cell*1000))}um_E{E}"
                J.append(("runs_m3", tag, "lattice",
                          dict(topo=topo, f=INFILL, cell=cell, spot=cell), E, nevt))
    return J


def macro(kind, p, E, out, seed, nevt):
    Lm = []
    if kind == "lattice":
        Lm += ["/mcs/det/geom voxel", f"/mcs/det/voxel {vox_path(p['topo'], p['f'], p['cell'])}",
               f"/mcs/det/topology {p['topo']}", f"/mcs/det/infill {p['f']}"]
    elif kind == "solid":
        Lm += ["/mcs/det/geom slab", "/mcs/det/material PLA",
               f"/mcs/det/thickness {p['thick']} mm"]
    else:  # empty frame -> floor(0)=0
        Lm += ["/mcs/det/geom slab", "/mcs/det/material G4_Galactic",
               f"/mcs/det/thickness {p['thick']} mm"]
    Lm += ["/mcs/det/maxStep 0.1 mm", "/run/setCut 0.05 mm", "/run/initialize",
           "/gun/particle proton", f"/gun/energy {E} MeV"]
    if kind == "lattice":
        Lm.append(f"/mcs/gun/spotXY {p['spot']} mm")
    Lm += [f"/random/setSeeds {seed} 0", f"/analysis/setFileName {out}",
           f"/run/printProgress {max(1, nevt//4)}", f"/run/beamOn {nevt}"]
    return "\n".join(Lm) + "\n"


def done(d, tag, nevt):
    runs = os.path.join(ROOT, "data", d)
    meta = os.path.join(runs, tag + ".root.meta.json")
    root = os.path.join(runs, tag + ".root")
    if not (os.path.exists(meta) and os.path.exists(root)):
        return False
    try:
        return json.load(open(meta)).get("n_events", 0) >= nevt
    except Exception:
        return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", choices=["A", "B200", "B1000"], default=None,
                    help="run only one block (A=1000 std cell, B200/B1000=wide cells)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--nevt-a", default=None, help="override M3-A n_events (pilot)")
    ap.add_argument("--nevt-b", default=None, help="override M3-B n_events (pilot)")
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()
    try:
        os.environ.setdefault("MCS_GIT_COMMIT", subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT).decode().strip())
    except Exception:
        pass
    commit = os.environ.get("MCS_GIT_COMMIT", "")[:7]

    def block_of(d, tag, E):
        if d == "runs":
            return "A"
        return "B200" if E == 200 else "B1000"

    for d, tag, kind, p, E, nevt in jobs():
        blk = block_of(d, tag, E)
        if args.only and blk != args.only:
            continue
        if blk == "A" and args.nevt_a:
            nevt = int(float(args.nevt_a))
        if blk in ("B200", "B1000") and args.nevt_b:
            nevt = int(float(args.nevt_b))
        if args.list:
            print(f"{'DONE' if done(d,tag,nevt) else 'todo'}  [{blk}] {d}/{tag} @ {nevt:.0e}")
            continue
        if done(d, tag, nevt):
            print(f"  skip {tag} (done)"); continue
        runs = os.path.join(ROOT, "data", d)
        os.makedirs(runs, exist_ok=True)
        if kind == "lattice":
            ensure_vox(p["topo"], p["f"], p["cell"])
        sd = seed_native(tag) if d == "runs" else seed_wide(tag)
        out = os.path.join(runs, tag)
        mac = os.path.join(runs, f"_{tag}.mac")
        open(mac, "w").write(macro(kind, p, E, out, sd, nevt))
        t0 = time.time()
        rc = subprocess.run([BIN, mac, str(args.threads)],
                            stdout=open(out + ".run.log", "w"),
                            stderr=subprocess.STDOUT).returncode
        wall = time.time() - t0
        os.remove(mac)
        with open(RUNLOG, "a") as f:
            f.write(f"M3_{blk},sim/run_m3.py,{tag},{sd},{nevt},{wall:.1f},"
                    f"{args.threads},{rc},data/{d}/{tag}.root,{commit},M3 3rd-energy+wide-Neff\n")
        print(f"  [{blk}] {tag:34s} {'OK' if rc==0 else f'RC={rc}'} wall={wall:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
