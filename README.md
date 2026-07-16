# Multiple Coulomb Scattering in Lattice and Foam Detector Structures Beyond the Homogeneous-Slab Approximation

Simulation, analysis, and manuscript source for a Geant4 study of how multiple Coulomb
scattering (MCS) in porous low-mass detector structures — 3D-printed lattices and reticulated
carbon foams — departs from the homogeneous equal-areal-density slab that is routinely assumed.
Using a validated Geant4 11.3.2 transport campaign over **five topologies** (rectilinear,
Schwarz-P, gyroid, diamond, and stochastic Voronoi) across a range of infills and proton
energies, the paper shows the slab approximation is only *half* valid: the kink-angle
distribution **width** is geometry-blind (topology-independent to 1.9 % after the realised-path
correction), but the distribution **shape** is not — porosity imprints an excess fourth cumulant
that follows a parameter-free `N_eff^-1` scaling in an effective cell count set by the structure's
chord autocovariance. The repository contains the Geant4 application, the analytic and voxel
geometry generators, a transport-aware ray-tracer, the cumulant-analysis and deconvolution
pipeline, every figure-generation script, and the LaTeX manuscript.

## Authors

| Author | ORCID |
|---|---|
| Palaash Gang (corresponding) | [0009-0006-7448-9488](https://orcid.org/0009-0006-7448-9488) |
| Raunav Mendiratta | [0009-0003-4426-2088](https://orcid.org/0009-0003-4426-2088) |
| Palash Bharti | [0009-0000-8926-5520](https://orcid.org/0009-0000-8926-5520) |
| Dhruv Punjabi | [0009-0002-6678-6192](https://orcid.org/0009-0002-6678-6192) |

## Installation

Everything runs inside a single Conda environment (name of record `g4highland`):

```bash
conda env create -f environment.yml
conda activate g4highland
```

This provides Geant4 11.3.2 (from conda-forge, with its data files), the C++17 toolchain to
build the simulation, and the full Python analysis stack. No manual or out-of-band installations
are required. The **fast reproduction path below needs only the Python subset** (NumPy, SciPy,
Matplotlib, scikit-image, pytest) — no Geant4 — and was verified from a clean environment.

## Reproduce

Run from the repository root (`make help` lists everything):

| Command | What it does | Needs Geant4 / new compute? |
|---|---|---|
| `make figures` | Regenerate all 11 publication figures from the committed `data/analysis` outputs (data plots) plus deterministic voxel geometry (the geometry panel) | No |
| `make test` | Run the analytic unit tests (`analysis/lib/test_theory.py`, 21 checks of the theory identities) | No |
| `make audit` | Independently re-verify every headline result against committed data via the claims ledger (`tools/check_ledger.py` over `claims_ledger.csv`) | No |
| `make all` | `figures` + `test` + `audit` — the complete no-new-compute reproduction | No |
| `make campaign` | Prints the **full Monte-Carlo campaign** entry points (`sim/run_campaign.sh`, `sim/campaign.py`, the high-statistics escalations) that produce the raw `.root` output | **Yes** |

## Computational requirements

- **Fast / cached path** (`make figures test audit all`): no special hardware. Runs in **minutes**
  on any laptop from the committed reduced data; produces every figure and re-verifies every
  headline number without any simulation.
- **Full simulation campaign** (`make campaign`): the Geant4 transport is single-threaded per
  event, **embarrassingly parallel** across dozens of independent (topology × fill × energy)
  jobs, and totals **approximately 310 core-hours** in aggregate (the figure recorded in the
  manuscript's computing acknowledgment: ~63 core-hours for the original seeded campaign plus
  ~245 for the high-statistics diamond/Voronoi and pitch/floor reruns). A compute-optimized cloud
  instance (dozens of cores) is recommended for full reproduction — this is **not** a laptop-scale
  job. The raw per-event `.root` output is regenerable and is **not** committed (see below).

## Results artifacts

All reduced/derived outputs needed to regenerate the figures and tables are committed. Key files:

| File | Contents |
|---|---|
| `data/analysis/e2_results_combined.json` | Headline 200 + 500 MeV per-configuration campaign: deconvolved excess `Δκ4` with bootstrap CIs, `N_eff`, `k_eff`, and width ratios (drives the collapse and width figures) |
| `data/analysis/e2_results_m3_1000.json` | 1000 MeV periodic-family collapse points |
| `data/analysis/abmax_results.json` | AB-MAX high-statistics (3×10⁷) diamond rerun and four-seed Voronoi realisation scatter |
| `data/analysis/transport_vs_g4.json` | Transport-aware ray-tracer vs full-Geant4 `Δκ4`/`Var(t)` closure |
| `data/analysis/e5_msc_systematic.json` | WentzelVI-vs-Urban multiple-scattering-model systematic |
| `data/analysis/e3_affordability.json` | Shape-channel proton budget and beam-time cost |
| `data/analysis/pt_hist.json` | Per-configuration `p(t)` line-integral distributions and mean-path ratios |
| `data/analysis/homog_boundary.json` | Homogenisation-boundary (`c_break`) data |
| `data/analysis/w_sensitivity.json` | Acceptance-window sensitivity of `Δκ4`/`γ2` |
| `data/analysis/confirm_100mev.json`, `e0_break.json` | 100 MeV break onset and the cell-break map |
| `data/analysis/e4_print_realism.json`, `e4_validate.json` | Print-realism (roughness, microporosity) stress tests |
| `data/analysis/foam_spotcheck.json`, `e5a_rawk4_check.json`, `geom_consistency.json` | Foam-scale replication, raw-`κ4` floor-sensitivity, voxel-vs-tessellation consistency |
| `results/M1/` | Reconstruction-impact study (track loss, large-angle-tail and impact-parameter excesses) |
| `results/M2/` | Orientation/tilt scan (`neff_theta.json`, `collapse_theta.json`) |
| `results/M3/` | Wide-cell scale-invariance (`collapse_wide.json`) |
| `results/M4/` | Carbon material-transfer check (`carbon.json`) |
| `results/N2/` | Wandering-path-length variance |
| `results/N3/` | Bayesian `N_eff`/`f` recovery posterior |
| `results/N4/` | Beam-optics/momentum-spread systematics |
| `results/phase0_repro/` | Fresh-seed reproduction of the 500 MeV headline campaign |
| `claims_ledger.csv` | Machine-checked ledger tying every reported number to the script/data file that produced it (re-verified by `make audit`) |
| `RUNLOG.csv` | Per-run provenance (seeds, configuration) for the seeded campaign |

## Data and provenance

This project generates its own geometry and simulation data — **no external datasets are fetched**.

- **Physics:** Geant4 11.3.2, physics list `FTFP_BERT` + `G4EmStandardPhysics_option4` with the
  WentzelVI + single-Coulomb multiple-scattering model (`sim/src/PhysicsList.cc`), locked and
  recorded per run.
- **Geometry:** the lattice/foam/TPMS structures are generated analytically and voxelised by the
  `geom/` scripts; nothing is downloaded.
- **Raw output:** per-event `.root` ntuples land in `data/runs/` and are **gitignored** (large and
  fully regenerable from the recorded seeds). All **reduced/derived** analysis outputs (the
  `data/analysis/*.json` and `results/**` above) **are** committed and are sufficient to regenerate
  every figure and table via the fast path.
- **Reproducibility record:** the 200 + 500 MeV headline campaign was re-run fresh with
  deterministic seeds (`sim/run_p7_repro.py`; see `results/phase0_repro/` and `DEVIATIONS.md`), so
  every figure input regenerates from tracked artifacts.

## Pinned environment

Exact pins from `environment.yml` (conda-forge):

| Package | Version |
|---|---|
| python | 3.12 |
| geant4 | 11.3.2 |
| numpy | 2.4.* |
| scipy | 1.17.* |
| matplotlib | 3.10.* |
| uproot | 5.7.* |
| awkward | 2.9.* |
| scikit-image | 0.26.* |
| trimesh | (latest on solve) |
| pytest | (latest on solve) |
| cmake | ≥ 3.16 |
| clang | C++17 (macOS arm64: clang 18) |

## How to cite

- **Paper:** P. Gang, R. Mendiratta, Palash Bharti, Dhruv Punjabi, *Multiple Coulomb Scattering in
  Lattice and Foam Detector Structures Beyond the Homogeneous-Slab Approximation*, Journal of
  Instrumentation (JINST), forthcoming — **[VOLUME/PAGE — pending acceptance]**.
- **Software and data:** archived on Zenodo — **[ZENODO CONCEPT DOI — pending deposit]**.

See `CITATION.cff` for machine-readable citation metadata.

## License

- **Code** — MIT License (`LICENSE`).
- **Data products** — Creative Commons Attribution 4.0 International, CC BY 4.0 (`LICENSE-data`).
