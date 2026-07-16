#!/usr/bin/env python3
"""write_e0_break_md.py -- render data/analysis/e0_break.md from e0_break.json
(Phase C output) + the Phase-D Geant4 confirmation, with the break table,
topology-contrast result, and the sub-printable extrapolation scope."""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
A = os.path.join(ROOT, "data", "analysis")
J = json.load(open(os.path.join(A, "e0_break.json")))
cells, momenta, topos = J["cells_mm"], J["momenta_MeV"], list(J["cell_break_mm"])
fdm, sla = J["print_FDM_mm"], J["print_SLA_mm"]


def br(name, p):
    b = J["cell_break_mm"][name][str(p)]
    return b


L = ["# E0 — the homogeneous-approximation break (transport-aware map)\n",
     f"Method: transport-aware ray-tracer on ANALYTIC geometry (no facets/voxels/"
     f"steps -> no resolution/MSC artifact), validated against full Geant4 at "
     f"printable cells (`transport_vs_g4.md`). Break = cell where "
     f"Var(t_actual)/Var(t_straight)=Δκ4_TA/Δκ4_LI departs from 1 by "
     f">{J['threshold']:.0%}. L={J['L_mm']} mm, f={J['f']}, "
     f"{J['n_proton']} protons/point, {J['steps_per_cell']} steps/cell, "
     f"{J['n_seed']} seeds.\n",
     "## Break table — cell_break(momentum) [µm]\n",
     "| p [MeV] | rectilinear | gyroid | y_rms(p) solid [µm] | vs FDM 0.5 mm |",
     "|--:|--:|--:|--:|--:|"]
for p in momenta:
    br_r, br_g = br("rectilinear", p), br("gyroid", p)
    yr = J["y_rms_solid_mm"][str(p)] * 1e3
    sr = f"{br_r*1e3:.0f}" if br_r else "<10"
    sg = f"{br_g*1e3:.0f}" if br_g else "<10"
    worst = max(br_r or 0, br_g or 0)
    L.append(f"| {p} | {sr} | {sg} | {yr:.0f} | {worst/fdm:.2f}× |")

L.append("\n## Ratio map (Var(t_act)/Var(t_str))\n")
for name in topos:
    L.append(f"**{name}**\n")
    L.append("| p\\cell | " + " | ".join(f"{c:g}" for c in cells) + " |")
    L.append("|--:|" + "|".join(["--:"] * len(cells)) + "|")
    for p in momenta:
        L.append(f"| {p} | " + " | ".join(
            f"{J['grid'][f'{name}|{p}|{c}']['ratio']:.2f}" for c in cells) + " |")
    L.append("")

# topology contrast
L.append("## Topology-contrast test\n")
diffs = []
for p in momenta:
    a, b = br("rectilinear", p), br("gyroid", p)
    if a and b:
        diffs.append(abs(a - b) / max(a, b))
        L.append(f"- {p} MeV: rect break {a*1e3:.0f} µm vs gyroid {b*1e3:.0f} µm "
                 f"(differ {abs(a-b)/max(a,b)*100:.0f}%)")
maxd = max(diffs) if diffs else 0
verdict = ("SAME cell within ~%.0f%% -> the break is a near-universal length-scale "
           "criterion set by the lateral-scattering scale, weakly topology-dependent"
           % (maxd * 100)) if maxd < 0.35 else (
    "DIFFERENT -> the break is topology-dependent")
L.append(f"\n**Topology contrast: {verdict}.**\n")

L.append("## Scope (Phase D honesty)\n")
L.append(f"The break map uses the transport-aware tool on analytic geometry. It is "
         f"validated against full Geant4 at printable cells (≥0.5 mm, see "
         f"`transport_vs_g4.md`): Var(t_actual) and Δκ4 agree within bootstrap error "
         f"there. Below Geant4's validated step regime (sub-printable, fine cells) "
         f"NO Geant4 measurement is claimed — the break there rests on the "
         f"transport-aware map, justified by the printable-cell agreement where the "
         f"regimes overlap. The FDM print limit is {fdm} mm; fine SLA ~{sla} mm.\n")
open(os.path.join(A, "e0_break.md"), "w").write("\n".join(L) + "\n")
print("wrote e0_break.md")
print("break table:")
for p in momenta:
    print(f"  {p} MeV: rect={br('rectilinear',p)}, gyroid={br('gyroid',p)} mm")
