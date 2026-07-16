# USAGE

How to reproduce the figures, tables, and manuscript, and a reference to the pipeline scripts.
For the project overview, installation, and computational requirements, see `README.md`.

## Quickstart

```bash
conda env create -f environment.yml     # one-time; provides Geant4 11.3.2 + the Python stack
conda activate g4highland

make all      # fast, no-new-compute: figures from committed data + unit tests + ledger audit
make test     # analytic theory-identity unit tests (21 checks)
make audit    # re-verify every headline result against committed data (claims ledger)
make help     # list all targets
```

None of the above runs a simulation. The full Monte-Carlo campaign (Geant4, ~310 core-hours) is a
separate path — `make campaign` prints its entry points.

## Worked example: reproduce the N_eff collapse figure

The central Result-2 figure is `figs/fig_neff_collapse_3energy.pdf` (the parameter-free
`N_eff^-1` collapse). It regenerates from committed data with **no simulation**:

```bash
conda activate g4highland
cd figs/src
python fig_neff_collapse_3energy.py
```

- **Reads (committed):** `data/analysis/e2_results_combined.json` (the 200 + 500 MeV per-config
  deconvolved `Δκ4`, `N_eff`, width ratios) and `data/analysis/e2_results_m3_1000.json` (the
  1000 MeV periodic points), via the shared `figs/src/figstyle.py` loader.
- **Writes:** `figs/fig_neff_collapse_3energy.pdf` — a vector PDF with embedded fonts.
- **Expected content:** all topologies at 200/500 MeV plus the periodic families at 1000 MeV
  collapse onto the derived `C = N_eff^-1` line; the pooled periodic OLS slope is
  `-0.932 [-0.992, -0.872]`; QC-excluded points (floor-dominated Voronoi, diamond consistency
  corner) are shown but daggered.

To rebuild **all** figures at once: `make figures` (from the repo root). The one figure that needs
geometry rather than analysis data — the geometry panel `fig_geometry.pdf` — triggers a
deterministic voxel regeneration (`geom/make_campaign_voxels.py`, ~seconds–minutes, still no
Geant4).

## Pipeline script reference

The pipeline runs geometry → simulation → cumulant analysis → figures. Fast-path (no-Geant4)
stages are marked ●; simulation stages ○.

| Script | Reads | Writes | Notes |
|---|---|---|---|
| ● `geom/make_campaign_voxels.py <topo> <fill>` | analytic topology definitions | `data/geom_stats/voxel/*.raw`, `*.meta` | Voxelises the 5 topologies × 4 infills at the 2.5 mm campaign cell (res=32) |
| ● `geom/transport_raytrace.py` | analytic material field `χ(x,y,z)` | (library) | Transport-aware ray-tracer in the small-angle MCS diffusion limit; no voxels/facets |
| ○ `sim/campaign.py` | config matrix; builds voxels + macros on demand | `data/runs/*.root` | Resumable driver: 5 topologies × 4 infills × 3 energies + baselines; runs the `mcs_sim` engine |
| ○ `sim/run_campaign.sh` | (locked physics) | `data/runs/*.root` | Highland-validation + `κ_M`-linearity campaign (energy × thickness on solid PLA + empty-frame) |
| ● `analysis/lib/kink_stats.py` | a run's `.root` ntuple (uproot) | (library) | Robust width estimators + acceptance-defined cumulants; Highland reference |
| ● `analysis/e2_analysis.py [--energies 200 500] [--nevt-min N]` | `data/runs/*.root` | `data/analysis/e2_results_*.json` | Cumulant analysis + all-order floor deconvolution + per-config QC + collapse fit |
| ● `analysis/e3_affordability.py` | `data/analysis/e2_results_*.json`, `data/runs/*.root` | `data/analysis/e3_affordability.{json,md}` | Shape-channel proton cost / beam-time (calibration-anchored) |
| ● `analysis/e5_msc_systematic.py` | `.root` (WentzelVI vs Urban) | `data/analysis/e5_msc_systematic.{json,md}` | Multiple-scattering-model systematic band |
| ● `analysis/lib/theory.py` | — | (library) | Closed-form cumulant / `N_eff` identities (unit-tested by `test_theory.py`) |
| ● `figs/src/fig_*.py` | `data/analysis/*.json`, `results/**` | `figs/fig_*.pdf` | One script per publication figure; shared style in `figstyle.py` |
| ● `figs/src/fig_mechanism.tex` | — | `figs/fig_mechanism.pdf` | Native TikZ/PGF scale-mixture schematic |
| ● `tools/check_ledger.py` | `claims_ledger.csv` + committed data | (exit 0 / non-zero) | Independent re-verification of headline results (each row maps a number to its source file) |

Options common to the analysis scripts: `--energies` selects the energy subset; `--nevt-min`
(default `500000`) is the per-config statistics floor used in QC; `--out-tag` labels intermediate
outputs. The figure scripts take no arguments (paths are resolved through `figstyle.py`).

## Expected behaviour (verifiable invariants)

Running the fast path against the committed data reproduces these project-specific invariants;
`make audit` checks them mechanically:

- **Width invariance (Result 1):** the mean-path-corrected width ratio `R_width` lies in
  `[0.981, 1.006]` — all 52 campaign configurations close within **1.9 %** of unity.
- **Scale-mixture closure (Result 2):** on the measured per-proton path variance, the excess
  fourth cumulant closes as `k_eff = Δκ4 / [3 a_eff² Var(t_pla)] = 1.00 ± 0.01` — the
  scale-mixture identity to the **percent level**.
- **N_eff collapse:** the periodic-family pooled log-log exponent is `-0.932 [-0.992, -0.872]`,
  energy-universal across 200–1000 MeV.
- **Deterministic reproduction:** given the random seeds recorded in `RUNLOG.csv`, the seeded
  campaign reproduces bit-for-bit; the fresh-seed re-run in `results/phase0_repro/` reproduces the
  500 MeV headline within its stated confidence intervals.
- **Ledger:** `tools/check_ledger.py` exits 0 — every reported number matches its committed
  source (`make audit`).
