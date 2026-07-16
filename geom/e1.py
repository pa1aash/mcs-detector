"""e1.py -- S4 / E1 geometry statistics and the N_eff collapse inputs.

For every topology x infill x cell size: ray-trace the beam-direction line
integral, compute p(t) cumulants, the autocovariance C(u), l_int, Var(t), and
BOTH N_eff forms (exact = f(1-f)L^2/Var(t); asymptotic = L/(2 l_int)). Writes
data/geom_stats/summary.csv and per-geometry raw .npz (campaign cell).

Cross-checks:
  * theory.py B7 shape-law identity  gamma2 = 3(1-f)/(f N_eff_exact)  == 3 Var/<t>^2
  * theory.py B8 (asymptotic) -> exact in the many-cell (small-cell) limit
  * theory.py B9 extruded limit -> N_eff = 1 for a z-constant field
  * D2: Voronoi foam solid-chord distribution vs an exponential (how Markov-like)

Run inside g4highland:  python geom/e1.py
"""
from __future__ import annotations

import csv
import os
import sys

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "analysis", "lib"))
import raytrace as rt          # noqa: E402
import theory as th            # noqa: E402
import topologies as topo      # noqa: E402

L = 10.0
CAMPAIGN_CELL = 2.5
INFILLS = (0.20, 0.30, 0.40, 0.50)
SWEEP_CELLS = (0.5, 1.0, 1.6, 2.5, 3.3, 4.0)     # 10/cell = 20..2.5 cells deep
OUT = os.path.join(os.path.dirname(HERE), "data", "geom_stats")


SPC = 80   # S5(rebuilt) Phase A: converged sampling resolution (resolution_convergence.md)


def trace_analytic(name, cell, f, nxy=SPC):
    p, _ = topo.tune_analytic(name, cell, f, n=max(80, SPC), tol=8e-4)
    chi, dz = topo.ray_chi_analytic(name, cell, p, L, nxy=nxy, dz=cell / SPC)
    return rt.stats_from_chi(chi, dz, L, corr_frac=0.9)


def trace_voronoi(cell, f, rng):
    vox = max(cell / 60.0, 0.02)
    box = (4 * cell, 4 * cell, L)
    gs = tuple(int(round(b / vox)) for b in box)
    delta = topo.voronoi_delta(cell, box, gs, rng)
    wall = float(np.quantile(delta, f))
    chi3 = (delta < wall)
    nx, ny, nz = chi3.shape
    rays = chi3.reshape(nx * ny, nz).astype(np.uint8)
    return rt.stats_from_chi(rays, box[2] / nz, L, corr_frac=0.9), rays, box[2] / nz


def foam_chord_fit(rays, dz):
    """solid-phase chord lengths along z (uncensored runs), exponential KS fit."""
    lengths = []
    for row in rays:
        # runs of 1s not touching either end (uncensored)
        idx = np.flatnonzero(np.diff(np.r_[0, row, 0]))
        starts, ends = idx[0::2], idx[1::2]
        for s, e in zip(starts, ends):
            if s > 0 and e < row.size:            # drop edge-censored chords
                lengths.append((e - s) * dz)
    lengths = np.asarray(lengths)
    if lengths.size < 50:
        return dict(n=int(lengths.size), mean=float("nan"), ks=float("nan"), p=float("nan"))
    scale = lengths.mean()
    ks, pval = stats.kstest(lengths, "expon", args=(0.0, scale))
    return dict(n=int(lengths.size), mean=float(scale), ks=float(ks), p=float(pval))


def main():
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(909)
    rows = []
    print(f"{'geometry':>11} {'f':>5} {'cell':>5} {'l_int':>7} {'Neff_ex':>8} "
          f"{'Neff_as':>8} {'Var_t':>9}")

    cells = sorted(set(SWEEP_CELLS) | {CAMPAIGN_CELL})
    for name in topo.ALL_TOPOS:
        for f in INFILLS:
            for cell in cells:
                if name == "voronoi":
                    s, rays, dz = trace_voronoi(cell, f, rng)
                else:
                    s = trace_analytic(name, cell, f)
                rows.append((name, s.f, cell, s.l_int, s.N_eff_exact,
                             s.N_eff_asymp, s.var_t, *s.kappa))
                if cell == CAMPAIGN_CELL:
                    np.savez(os.path.join(OUT, f"{name}_f{int(f*100):02d}.npz"),
                             u=s.u, C=s.C, t_mean=s.t_mean, var_t=s.var_t,
                             f=s.f, l_int=s.l_int, N_eff_exact=s.N_eff_exact,
                             N_eff_asymp=s.N_eff_asymp, kappa=np.array(s.kappa))
                    print(f"{name:>11} {s.f:>5.2f} {cell:>5.2f} {s.l_int:>7.4f} "
                          f"{s.N_eff_exact:>8.2f} {s.N_eff_asymp:>8.2f} {s.var_t:>9.3e}")

    with open(os.path.join(OUT, "summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["geometry", "f_designed", "cell_size", "l_int",
                    "N_eff_exact", "N_eff_asymp", "Var_t", "k1", "k2", "k3", "k4"])
        for r in rows:
            w.writerow([r[0], f"{r[1]:.4f}", r[2], f"{r[3]:.5f}", f"{r[4]:.4f}",
                        f"{r[5]:.4f}", f"{r[6]:.6e}", f"{r[7]:.6e}", f"{r[8]:.6e}",
                        f"{r[9]:.6e}", f"{r[10]:.6e}"])

    # ---- cross-checks ----
    print("\n--- theory.py cross-checks ---")
    # B7 identity on the campaign points
    b7max = 0.0
    for r in rows:
        f, ne, var = r[1], r[4], r[6]
        if not np.isfinite(ne) or var <= 0:
            continue                      # projection-homogeneous (Var(t)->0): N_eff->inf
        g2_a = 3.0 * var / (f * L) ** 2
        g2_b = th.gamma2(f, ne)
        if g2_b > 0:
            b7max = max(b7max, abs(g2_a / g2_b - 1.0))
    print(f"B7 shape-law identity gamma2=3Var/<t>^2 == 3(1-f)/(f Neff): "
          f"max rel diff = {b7max:.2e}  ({'OK' if b7max < 1e-6 else 'CHECK'})")

    # B9 extruded: z-constant field -> N_eff = 1
    chi_ext = np.tile((rng.random((4000, 1)) < 0.4).astype(np.uint8), (1, 400))
    s_ext = rt.stats_from_chi(chi_ext, L / 400, L)
    print(f"B9 extruded (z-constant, f~0.4): N_eff_exact = {s_ext.N_eff_exact:.3f} "
          f"(theory 1.0)  ({'OK' if abs(s_ext.N_eff_exact-1) < 0.05 else 'CHECK'})")

    # B8 asymptotic: exact vs asymp converge as cell shrinks (gyroid f0.4)
    gy = [(r[2], r[4], r[5]) for r in rows if r[0] == "gyroid" and abs(r[1]-0.4) < 0.03]
    gy.sort()
    print("B8 asymptotic convergence (gyroid f0.4): cell -> Neff_exact / Neff_asymp")
    for cell, ne, na in gy:
        print(f"    cell={cell:>4} Neff_ex={ne:>7.2f} Neff_as={na:>7.2f} ratio={ne/na:.3f}")

    # D2 foam Markov-ness (campaign cell, f0.4)
    _, rays, dz = trace_voronoi(CAMPAIGN_CELL, 0.40, np.random.default_rng(11))
    fit = foam_chord_fit(rays, dz)
    print(f"\nD2 Voronoi-foam solid-chord exponential fit (cell {CAMPAIGN_CELL}, f0.4): "
          f"n={fit['n']} mean={fit['mean']:.3f} mm KS={fit['ks']:.3f} p={fit['p']:.3f}")
    with open(os.path.join(OUT, "foam_markov_fit.txt"), "w") as fh:
        fh.write(f"voronoi solid-chord exponential KS fit (cell={CAMPAIGN_CELL}, f=0.40)\n")
        fh.write(f"n={fit['n']} mean_chord_mm={fit['mean']:.4f} "
                 f"KS_stat={fit['ks']:.4f} p_value={fit['p']:.4f}\n")

    nex = [r[4] for r in rows]
    print(f"\nsummary.csv: {len(rows)} rows; N_eff_exact range "
          f"[{min(nex):.2f}, {max(nex):.2f}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
