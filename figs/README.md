# Publication figures

`figs/` contains the eleven publication figures accompanying the paper.  Every
publication PDF is generated from a tracked source in `figs/src/`; the graph
scripts share `figstyle.py`, while the mechanism is native TikZ/PGFPlots.  Regenerate
them all from the committed data with, from the repository root,

```bash
make figures
```

The Python environment is recorded in [`environment.yml`](../environment.yml).

## Active figure set

| Order | File | Role and evidence | Primary source |
|---:|---|---|---|
| 1 | `fig_mechanism.pdf` | Straight-chord Gaussian scale-mixture mechanism; explanatory vector schematic | `src/fig_mechanism.tex` |
| 2 | `fig_geometry.pdf` | Five matched as-built topology cutaways; vector surfaces from the raw voxel fields | `src/fig_geometry.py`, `data/geom_stats/voxel/` |
| 3 | `fig_pt.pdf` | Literal straight-chord line-integral probabilities at a common designed mean | `data/analysis/pt_hist.json` |
| 4 | `fig_validation.pdf` | Solid-control Highland trend and Geant4/transport closure | `data/analysis/transport_vs_g4.json` plus committed solid-control values |
| 5 | `fig_width_invariance.pdf` | Main “half-valid” result: realised-path width closure and direct shape-law closure | `data/analysis/e2_results_combined.json`, `e2_results_m3_1000.json` |
| 6 | `fig_neff_collapse_3energy.pdf` | Three-energy parameter-free effective-cell-count collapse and residual closure | same seeded E2 products |
| 7 | `fig_boundary.pdf` | Practical cell-size evidence and explicitly extrapolated shape-recovery region | `homog_boundary.json`, `e0_break.json`, `e5_msc_systematic.json` |
| 8 | `fig_rocking.pdf` | Orientation prediction and Geant4 closure through the simulated tilt range | `results/M2/collapse_theta.json` |
| 9 | `fig_impact.pdf` | Kink-angle and track-fit consequences of using homogeneous process noise | `results/M1/f6_hist.json` |
| 10 | `fig_affordability.pdf` | Fourth-moment statistical cost and parameter-recovery precision | `e3_affordability.json`, `results/N3/posterior.json` |
| 11 | `fig_systematics.pdf` | Physics-model and fabrication-systematic stress tests | `e4_print_realism.json`, `e4_validate.json`, `e5_msc_systematic.json` |

Each publication PDF above is built by `make figures`.

## Visual conventions

- topology is redundant in colour, marker shape, and line pattern;
- energy is encoded by marker fill;
- physics-list changes use edge/line style;
- all visible text is at least 8.5 pt at final JINST display size;
- evidence, transport prediction, and extrapolation use distinct visual treatments;
- PDFs have opaque white backgrounds and are checked for raster soft masks.

The table above records each figure's role and its committed data source; the authoritative
captions appear in the accompanying paper.
