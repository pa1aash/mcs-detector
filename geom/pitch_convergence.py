"""pitch_convergence.py -- Phase 2 Task 2: geometry convergence at the REAL campaign
operating pitches (audit CONV-PITCH).

The campaign builds voxel fields at a TOPOLOGY-DEPENDENT pitch (geom/make_campaign_voxels.py):
rectilinear at spc=48 (52 um), the four others (gyroid, Schwarz-P, diamond, Voronoi) at
spc=32 (78 um). The original resolution_convergence.py only tested spc={40,60,80,120}, so
neither real operating pitch was ever on the grid, and 78 um is COARSER than the whole tested
range. This script traces f and Var(t) at spc = {32,40,48,60,80,120} for all five topologies
and asks, at each topology's real operating pitch, whether f and Var(t) are converged to <1%
against the finest (spc=120) reference.

Geometry-only (ray-tracer + analytic/Voronoi generator); no Geant4, no .root. Runs anywhere.
Voronoi uses the SAME seed as the committed campaign f=0.40 field (make_campaign_voxels SEED
6006 + round(100*f) = 6046) so the test is on the realization the paper actually used.
Writes data/geom_stats/pitch_convergence.md.
"""
from __future__ import annotations
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "analysis", "lib"))
import raytrace as rt          # noqa: E402
import topologies as topo      # noqa: E402

L, CELL, F = 10.0, 2.5, 0.40
SPCS = [32, 40, 48, 60, 80, 120]
OPPITCH = {"rectilinear": 48, "gyroid": 32, "schwarzp": 32, "diamond": 32, "voronoi": 32}
CAMP_SEED = 6006 + int(round(100 * F))     # 6046, matches make_campaign_voxels f=0.40
OUT = os.path.join(os.path.dirname(HERE), "data", "geom_stats")


def trace_analytic(name, spc):
    p, _ = topo.tune_analytic(name, CELL, F, n=max(80, spc), tol=8e-4)
    chi, dz = topo.ray_chi_analytic(name, CELL, p, L, nxy=spc, dz=CELL / spc)
    s = rt.stats_from_chi(chi, dz, L, corr_frac=0.9)
    return s.f, s.var_t, s.N_eff_exact


def trace_voronoi(spc):
    box = (3 * CELL, 3 * CELL, L)
    vx = CELL / spc
    gs = (int(round(box[0] / vx)), int(round(box[1] / vx)), int(round(L / vx)))
    chi, _, _ = topo.voronoi_field(CELL, box, F, gs, np.random.default_rng(CAMP_SEED))
    nx, ny, nz = chi.shape
    s = rt.stats_from_chi(chi.reshape(nx * ny, nz).astype(np.uint8), L / nz, L, corr_frac=0.9)
    return s.f, s.var_t, s.N_eff_exact


def main():
    series = {}
    for name in ("rectilinear", "gyroid", "schwarzp", "diamond", "voronoi"):
        series[name] = []
        for spc in SPCS:
            f, var, neff = trace_voronoi(spc) if name == "voronoi" else trace_analytic(name, spc)
            series[name].append((spc, f, var, neff))

    lines = ["# Voxel-pitch convergence at the REAL campaign operating pitches (Phase 2 Task 2)\n",
             f"Campaign cell {CELL} mm, f={F}, L={L} mm. Pitch = {CELL} mm / spc. Operating pitches: "
             "rectilinear spc=48 (52 um); gyroid/Schwarz-P/diamond/Voronoi spc=32 (78 um). "
             "Voronoi = committed campaign realization (seed 6046). Convergence reference = spc=120; "
             "PASS = |dev| < 1% in BOTH f and Var(t) at the operating pitch.\n"]
    verdict = {}
    for name, ser in series.items():
        op = OPPITCH[name]
        f_ref, v_ref = ser[-1][1], ser[-1][2]                 # spc=120 reference
        row_op = next(r for r in ser if r[0] == op)
        df = row_op[1] / f_ref - 1.0
        dv = (row_op[2] / v_ref - 1.0) if v_ref != 0 else float("nan")
        ok = (abs(df) < 0.01) and (abs(dv) < 0.01)
        verdict[name] = dict(op_pitch_spc=op, op_pitch_um=1000 * CELL / op,
                             f_op=row_op[1], var_op=row_op[2], f_ref=f_ref, var_ref=v_ref,
                             df_pct=100 * df, dv_pct=100 * dv, pass_=bool(ok))
        lines.append(f"## {name}  (operating pitch spc={op}, {1000*CELL/op:.0f} um)\n")
        lines.append("| spc | pitch [um] | f | Var(t) [mm^2] | N_eff | df vs 120 | dVar vs 120 |")
        lines.append("|--:|--:|--:|--:|--:|--:|--:|")
        for spc, f, var, neff in ser:
            mark = "  <-- OPERATING" if spc == op else ""
            lines.append(f"| {spc}{mark} | {1000*CELL/spc:.0f} | {f:.4f} | {var:.4f} | {neff:.3f} "
                         f"| {100*(f/f_ref-1):+.2f}% | {100*(var/v_ref-1) if v_ref else float('nan'):+.2f}% |")
        lines.append(f"\n**Operating pitch: f dev {100*df:+.2f}%, Var(t) dev {100*dv:+.2f}% vs spc=120 "
                     f"-> {'PASS (<1%)' if ok else 'FAIL (>=1%)'}**")
        if name == "diamond":
            lines.append(f"\n_Note: diamond Var(t) is near-zero ({v_ref:.4f} mm^2 at spc=120), the "
                         "N_eff->inf corner; relative Var deviations are amplified by the tiny "
                         "denominator and are the known numerical-floor issue, distinct from a "
                         "genuine coarse-pitch geometric failure -- read the absolute Var(t) column._")
        lines.append("")

    open(os.path.join(OUT, "pitch_convergence.md"), "w").write("\n".join(lines) + "\n")
    json.dump(verdict, open(os.path.join(OUT, "pitch_convergence.json"), "w"), indent=1)

    print("=== operating-pitch convergence (vs spc=120 reference) ===")
    for name, v in verdict.items():
        print(f"  {name:12s} spc={v['op_pitch_spc']:>3} ({v['op_pitch_um']:.0f}um): "
              f"f {v['df_pct']:+.2f}%  Var {v['dv_pct']:+.2f}%  -> {'PASS' if v['pass_'] else 'FAIL'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
