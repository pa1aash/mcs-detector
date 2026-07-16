"""resolution_convergence.py -- S5(rebuilt) Phase A1.

For each topology at the campaign cell (2.5 mm, f=0.40), ray-trace f and Var(t)
at increasing sampling resolution spc in {40,60,80,120} (spc = transverse rays
per cell AND longitudinal steps per cell, dz=cell/spc). Find the converged spc
where BOTH f and Var(t) change <1% between successive resolutions. Each topology
may converge at a different spc.

Resolution here is the geometry SAMPLING the ray-tracer (and the campaign voxel
build) uses; the analytic field is exact, so this measures how finely it must be
sampled to pin f and Var(t). Writes data/geom_stats/resolution_convergence.md.
"""
from __future__ import annotations
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "analysis", "lib"))
import raytrace as rt          # noqa: E402
import topologies as topo      # noqa: E402

L = 10.0
CELL = 2.5
F = 0.40
SPCS = [40, 60, 80, 120]
TOL = 0.01
OUT = os.path.join(os.path.dirname(HERE), "data", "geom_stats")


def trace_analytic(name, spc):
    p, _ = topo.tune_analytic(name, CELL, F, n=max(80, spc), tol=8e-4)
    chi, dz = topo.ray_chi_analytic(name, CELL, p, L, nxy=spc, dz=CELL / spc)
    s = rt.stats_from_chi(chi, dz, L, corr_frac=0.9)
    return s.f, s.var_t, s.l_int, s.N_eff_exact


def trace_voronoi(spc, rng):
    # voxel grid at resolution spc per cell over a 3x3xL block; rays = z-columns
    box = (3 * CELL, 3 * CELL, L)
    gs = (int(round(box[0] / (CELL / spc))), int(round(box[1] / (CELL / spc))),
          int(round(L / (CELL / spc))))
    chi, _, _ = topo.voronoi_field(CELL, box, F, gs, rng)
    nx, ny, nz = chi.shape
    rays = chi.reshape(nx * ny, nz).astype(np.uint8)
    dz = L / nz
    s = rt.stats_from_chi(rays, dz, L, corr_frac=0.9)
    return s.f, s.var_t, s.l_int, s.N_eff_exact


def converged_spc(series):
    """series: list of (spc,f,var). Return first spc where |Δf|,|Δvar| <TOL vs prev."""
    for i in range(1, len(series)):
        _, f0, v0 = series[i - 1]
        spc, f1, v1 = series[i]
        if abs(f1 / f0 - 1) < TOL and abs(v1 / v0 - 1) < TOL:
            return spc
    return None


def main():
    rng = np.random.default_rng(4040)
    rows = {}
    for name in ("rectilinear", "gyroid", "schwarzp", "diamond", "voronoi"):
        series = []
        for spc in SPCS:
            if name == "voronoi":
                f, var, lint, neff = trace_voronoi(spc, np.random.default_rng(4040))
            else:
                f, var, lint, neff = trace_analytic(name, spc)
            series.append((spc, f, var, lint, neff))
        rows[name] = series

    lines = ["# Resolution convergence (Phase A1)\n",
             f"Campaign cell {CELL} mm, f={F}, L={L} mm. spc = transverse rays per "
             f"cell and longitudinal steps per cell (dz=cell/spc). Converged = both "
             f"f and Var(t) change <{TOL:.0%} between successive spc.\n"]
    conv = {}
    for name, series in rows.items():
        lines.append(f"## {name}\n")
        lines.append("| spc | f | Var(t) | l_int | N_eff_exact | Δf | ΔVar |")
        lines.append("|--:|--:|--:|--:|--:|--:|--:|")
        for i, (spc, f, var, lint, neff) in enumerate(series):
            if i == 0:
                df = dv = "—"
            else:
                df = f"{f/series[i-1][1]-1:+.2%}"
                dv = f"{var/series[i-1][2]-1:+.2%}"
            lines.append(f"| {spc} | {f:.4f} | {var:.4f} | {lint:.4f} | {neff:.4f} "
                         f"| {df} | {dv} |")
        c = converged_spc([(s[0], s[1], s[2]) for s in series])
        conv[name] = c
        lines.append(f"\n**converged spc = {c if c else 'NOT by 120 (FLAG)'}**\n")

    lines.append("## Converged spc per topology\n")
    lines.append("| topology | converged spc |\n|--|--:|")
    for name, c in conv.items():
        lines.append(f"| {name} | {c if c else '>120 FAIL'} |")
    with open(os.path.join(OUT, "resolution_convergence.md"), "w") as f:
        f.write("\n".join(lines) + "\n")

    print("converged spc:", conv)
    for name, series in rows.items():
        print(f"  {name:12s}: " + "  ".join(
            f"spc{s[0]}:f={s[1]:.4f},Var={s[2]:.3f}" for s in series))
    # persist converged spc for A2
    import json
    with open(os.path.join(OUT, "converged_spc.json"), "w") as f:
        json.dump({k: (v if v else 120) for k, v in conv.items()}, f)
    return 0 if all(conv.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
