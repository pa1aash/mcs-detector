#!/usr/bin/env python3
"""run_m4_solids_ext.py -- two extra carbon solid controls (t=8,16 mm; 200 MeV; 3e6)
so the M4 floor fit covers the full per-proton path range (fix for the floor
under-coverage found by the math-referee pass; see DEVIATIONS D8).
Reuses sim/run_m2m4.py machinery (tags, seeds, macro, RUNLOG line).
"""
import os, sys, time, subprocess, json
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import run_m2m4 as m

THREADS = 4
for t in (8, 16):
    tag = f"m4_csolid_t{t}"
    if m.done("runs_m4", tag, 3_000_000):
        print("skip", tag); continue
    sd = m.seed_for(tag)
    runs = os.path.join(ROOT, "data", "runs_m4")
    out = os.path.join(runs, tag)
    mac = os.path.join(runs, f"_{tag}.mac")
    open(mac, "w").write(m.macro("solid", dict(mat="C_amorph", thick=t), out, sd, 3_000_000))
    t0 = time.time()
    rc = subprocess.run([m.BIN, mac, str(THREADS)], stdout=open(out + ".run.log", "w"),
                        stderr=subprocess.STDOUT).returncode
    wall = time.time() - t0
    os.remove(mac)
    with open(m.RUNLOG, "a") as f:
        f.write(f"M4,sim/run_m4_solids_ext.py,{tag},{sd},3000000,{wall:.1f},{THREADS},{rc},"
                f"data/runs_m4/{tag}.root,{os.environ.get('MCS_GIT_COMMIT','working')[:7]},"
                f"D8 floor-range extension\n")
    print(tag, "OK" if rc == 0 else f"RC={rc}", f"wall={wall:.0f}s")
