#!/usr/bin/env python3
"""break_sensitivity.py -- S5-verify Phase 2. Re-process the transport break map
at thresholds {5,10,20,30}% -> c(p, topology, threshold), and tabulate cell_break
vs FDM nozzle sizes {0.25,0.4,0.5,0.8} mm and SLA (~0.05 mm). Shows how the break
location and its printability classification move with the (arbitrary) threshold
and the (vendor-dependent) nozzle floor. CPU-light: reads e0_break.json.
"""
from __future__ import annotations
import csv, json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
A = os.path.join(ROOT, "data", "analysis")
J = json.load(open(os.path.join(A, "e0_break.json")))
cells = J["cells_mm"]; momenta = J["momenta_MeV"]; topos = list(J["cell_break_mm"])
THRESH = [0.05, 0.10, 0.20, 0.30]
FDM = [0.25, 0.4, 0.5, 0.8]; SLA = 0.05


def break_at(name, p, thr):
    curve = sorted([(c, J["grid"][f"{name}|{p}|{c}"]["ratio"]) for c in cells],
                   key=lambda t: -t[0])
    tgt = 1.0 - thr
    for i in range(1, len(curve)):
        (c0, r0), (c1, r1) = curve[i - 1], curve[i]
        if r0 >= tgt and r1 < tgt:
            return float(np.exp(np.log(c0) + (tgt - r0) *
                                (np.log(c1) - np.log(c0)) / (r1 - r0)))
    return None if curve[-1][1] >= tgt else curve[0][0]


# c(p, topology, threshold)
rows = []
for name in topos:
    for p in momenta:
        for thr in THRESH:
            b = break_at(name, p, thr)
            rows.append(dict(topology=name, p=p, threshold=thr,
                             cell_break_um=(b * 1e3 if b else None)))
with open(os.path.join(A, "break_sensitivity.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["topology", "p", "threshold", "cell_break_um"])
    w.writeheader(); w.writerows(rows)

# markdown
L = ["# Break threshold + nozzle sensitivity (Phase 2)\n",
     "## c(p, topology, threshold) [µm] — how cell_break moves with the threshold\n"]
for name in topos:
    L.append(f"**{name}**\n")
    L.append("| p [MeV] | 5% | 10% | 20% | 30% |")
    L.append("|--:|--:|--:|--:|--:|")
    for p in momenta:
        vals = []
        for thr in THRESH:
            b = break_at(name, p, thr)
            vals.append(f"{b*1e3:.0f}" if b else "<30")
        L.append(f"| {p} | " + " | ".join(vals) + " |")
    L.append("")

# printability classification vs nozzle (at 10% and the 5-20% band)
L.append("## cell_break vs nozzle floor — does the break land in printable range?\n")
L.append("Entry = ✓ if cell_break ≥ nozzle (the break sits at/above a printable "
         "feature, i.e. the criterion bites for that printer); — otherwise "
         "(break is finer than the printer can make → design-rule for that printer).\n")
L.append("| topo | p | c_break(10%) µm | FDM 0.25 | FDM 0.4 | FDM 0.5 | FDM 0.8 | SLA 0.05 |")
L.append("|--|--:|--:|:--:|:--:|:--:|:--:|:--:|")
for name in topos:
    for p in momenta:
        b = break_at(name, p, 0.10)
        if not b:
            continue
        bu = b * 1e3
        cells_cls = [("✓" if b >= n else "—") for n in FDM] + \
                    ["✓" if b >= SLA else "—"]
        L.append(f"| {name} | {p} | {bu:.0f} | " + " | ".join(cells_cls) + " |")

# robustness across the 5-20% band (Phase-4 anti-fishing input)
L.append("\n## Robustness band (5–20% threshold) for the headline rect@100 MeV\n")
for name in ["rectilinear", "gyroid"]:
    b5, b10, b20 = (break_at(name, 100, t) for t in (0.05, 0.10, 0.20))
    L.append(f"- {name} @100 MeV: c_break = {b5*1e3:.0f} / {b10*1e3:.0f} / "
             f"{b20*1e3:.0f} µm (5/10/20%). FDM-0.5 overlap "
             f"(c_break ≥ 0.5 mm): "
             f"{'/'.join('Y' if b>=0.5 else 'N' for b in (b5,b10,b20))}.")
open(os.path.join(A, "break_sensitivity.md"), "w").write("\n".join(L) + "\n")

print("=== c(p,topology,threshold) [µm] ===")
for name in topos:
    for p in momenta:
        vals = [f"{(break_at(name,p,t) or 0)*1e3:.0f}" for t in THRESH]
        print(f"  {name:11s} {p:4d}: 5%={vals[0]} 10%={vals[1]} 20%={vals[2]} 30%={vals[3]}")
print("\nrect@100 FDM-0.5 overlap across 5/10/20%:",
      "/".join("Y" if (break_at("rectilinear",100,t) or 0) >= 0.5 else "N"
               for t in (0.05, 0.10, 0.20)))
