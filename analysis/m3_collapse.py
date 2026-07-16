#!/usr/bin/env python3
"""m3_collapse.py -- M3 wide-N_eff collapse analysis (block B) + 3-energy summary.

Block A (1000 MeV, standard 2.5 mm cell) is analysed by the canonical frozen path
`e2_analysis.py --energies 1000 --out-tag m3_1000` (native tags in data/runs/). This script
handles the wide-cell block B (data/runs_m3/) and assembles the multi-energy statement.

For each wide-cell run m3_<topo>_f30_c<cell>um_E<E>.root it applies the FROZEN procedures
(ANALYSIS_PLAN sec 0 + M3): N_eff from the as-built field ray-traced at that cell (same
estimator as e2_analysis.geom_asbuilt), Delta_kappa4 by the all-order floor subtraction
(200 MeV: committed a_eff(200) + the M2 PLA-solid floor; 1000 MeV: the block-A 1000 solids
via e2_analysis), and C = Delta_kappa4/(3 a_eff^2 f(1-f) L^2). The wide cells span N_eff from
~0.5 (5 mm cell) to ~50 (0.5 mm cell) = >=2 decades, giving the tightest single-energy
exponent. Every cell is checked against the c_break guard (cell >= 5 y_rms) so no
wandering-regime point enters the fit.

Run inside g4highland (after M3 finishes): python analysis/m3_collapse.py
"""
from __future__ import annotations
import glob, json, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HERE, "lib"))
sys.path.insert(0, os.path.join(ROOT, "geom"))
import kink_stats as ks   # noqa
import raytrace as rt     # noqa
import topologies as topo  # noqa
import e2_analysis as e2  # noqa  (a_eff_of_E, build_floor for 1000 MeV)

RUNS_M3 = os.path.join(ROOT, "data", "runs_m3")
RUNS_M2 = os.path.join(ROOT, "data", "runs_m2")
OUT = os.path.join(ROOT, "results", "M3")
L, F = 10.0, 0.30
WIN = {200: 37.84e-3, 1000: 8.95e-3}     # locked absolute windows
D4 = {200: -0.167, 1000: -0.196}         # amplitude renorm (exponent-preserving)
AEFF200 = 3.8743e-6                       # committed (e2_results_combined)
A_OF_P = {100: 1.7197e-05, 200: 3.9146e-06, 500: 7.2058e-07, 1000: 2.1943e-07}  # solid-path a
# committed per-energy OLS collapse exponents (200+500 pooled campaign; e2_results_combined)
COMMITTED = {"200": (-0.951, -1.040, -0.866), "500": (-1.001, -1.095, -0.908),
             "200+500 pooled": (-0.983, -1.052, -0.914)}


def floor200():
    """kappa_M(t)=b t + c t^2 from the M2 PLA solids (the 200 MeV floor); returns kM(t)."""
    kt, kv = [0.0], [0.0]
    for t in (2, 3, 4, 5):
        p = os.path.join(RUNS_M2, f"m2_solid_t{t}.root")
        if os.path.exists(p):
            _, k4 = ks.cumulants_in_window(ks.load_run(p).angles, WIN[200])
            kt.append(float(t)); kv.append(k4)
    kt, kv = np.array(kt), np.array(kv)
    b, c = np.linalg.lstsq(np.vstack([kt, kt ** 2]).T, kv, rcond=None)[0]
    tmax = kt.max()
    return lambda tt: b * np.clip(tt, 0, tmax) + c * np.clip(tt, 0, tmax) ** 2


def aeff_floor(E):
    if E == 200:
        return AEFF200, floor200()
    return e2.a_eff_of_E(E), e2.build_floor(E)         # 1000: reads data/runs/solid_E1000_t*


def wide_geom(topo_name, cell):
    """(f, Var(t), N_eff) of the as-built field ray-traced at this cell (e2 estimator)."""
    p, _ = topo.tune_analytic(topo_name, cell, F, n=80, tol=8e-4)
    chi, dz = topo.ray_chi_analytic(topo_name, cell, p, L, nxy=80, dz=cell / 40.0)
    s = rt.stats_from_chi(chi, dz, L, corr_frac=0.9)
    return s.f, s.var_t, s.N_eff_exact


def y_rms(E):
    """Transverse MCS wander over L for the homogenised medium (a_g = a f)."""
    return float(np.sqrt(A_OF_P[E] * F * L ** 3 / 3.0))


def main():
    os.makedirs(OUT, exist_ok=True)
    files = sorted(glob.glob(os.path.join(RUNS_M3, "m3_*_f30_c*um_E*.root")))
    if not files:
        print("no wide-cell runs found in data/runs_m3 yet (M3 block B not finished).")
        return 1
    rows = []
    for fp in files:
        m = re.match(r"m3_(\w+?)_f30_c(\d+)um_E(\d+)", os.path.basename(fp))
        if not m:
            continue
        topo_name, cell_um, E = m.group(1), int(m.group(2)), int(m.group(3))
        cell = cell_um / 1000.0
        # skip incomplete / mid-write runs (meta below the per-energy production target)
        target = 3e6 if E == 200 else 1e7
        meta_p = fp + ".meta.json"
        try:
            nev = json.load(open(meta_p)).get("n_events", 0) if os.path.exists(meta_p) else 0
        except Exception:
            nev = 0
        if nev < 0.9 * target:
            print(f"  skip {os.path.basename(fp)} (n_events {nev:.0e} < target)")
            continue
        aeff, kM = aeff_floor(E)
        rd = ks.load_run(fp)
        import uproot
        with uproot.open(fp) as fr:
            tpla = fr["kinks"]["tpla"].array(library="np")
        k2, k4 = ks.cumulants_in_window(rd.angles, WIN[E])
        floor_mean = float(np.mean(kM(tpla)))
        dk4 = (k4 - floor_mean) / (1.0 + D4[E])
        fb, var_t, neff = wide_geom(topo_name, cell)
        C = dk4 / (3.0 * aeff ** 2 * fb * (1 - fb) * L ** 2) if dk4 > 0 else float("nan")
        yr = y_rms(E)
        rows.append(dict(topo=topo_name, E=E, cell_mm=cell, f_built=fb, N_eff=neff,
                         dk4=float(dk4), C=float(C), CxNeff=float(C * neff),
                         y_rms_mm=yr, cell_over_5yrms=cell / (5 * yr),
                         straight_chord_ok=bool(cell >= 5.0 * yr)))
    # SCALE INVARIANCE (D6): for periodic lattices N_eff is cell-independent, so the wide-cell
    # points are replicates at one N_eff per topology. The test is that (i) N_eff does not move
    # with cell and (ii) the collapse point C (and C*N_eff) is stationary across the 10x cell range.
    scale_inv = {}
    for E in sorted({r["E"] for r in rows}):
        for topo_name in sorted({r["topo"] for r in rows if r["E"] == E}):
            g = [r for r in rows if r["E"] == E and r["topo"] == topo_name
                 and r["dk4"] > 0 and np.isfinite(r["C"])]
            if len(g) < 2:
                continue
            ne = np.array([r["N_eff"] for r in g]); Cv = np.array([r["C"] for r in g])
            scale_inv[f"{topo_name}@{E}"] = dict(
                cells_mm=[r["cell_mm"] for r in g],
                N_eff_mean=float(ne.mean()), N_eff_spread_pct=float(100 * ne.std() / ne.mean()),
                C_mean=float(Cv.mean()), C_cv_pct=float(100 * Cv.std() / Cv.mean()),
                CxNeff_mean=float(np.mean(ne * Cv)))
    out = dict(block="B (cell scale-invariance; N_eff cell-independent for periodic lattices, D6)",
               L=L, f=F, window_mrad=WIN, committed_per_energy_exponents=COMMITTED,
               wide_cell_rows=rows, scale_invariance=scale_inv,
               note="Block A (1000 MeV std cell, the third energy) -> e2_analysis.py "
                    "--energies 1000 --out-tag m3_1000; the collapse EXPONENT comes from the "
                    "standard-cell topology span, not these replicate wide-cell points")
    json.dump(out, open(os.path.join(OUT, "collapse_wide.json"), "w"), indent=1)
    print(f"M3 block B (scale invariance): {len(rows)} runs")
    for k, s in scale_inv.items():
        print(f"  {k:18s} N_eff={s['N_eff_mean']:.2f} (spread {s['N_eff_spread_pct']:.1f}% over "
              f"cells {s['cells_mm']}), C={s['C_mean']:.3f} (CV {s['C_cv_pct']:.0f}%), "
              f"C*N_eff={s['CxNeff_mean']:.2f}")
    print(f"wrote {OUT}/collapse_wide.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
