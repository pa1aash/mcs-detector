# configs/ — version-controlled campaign configs (work order Ground Rule 9)

One YAML per campaign (`M1.yaml`, `M2.yaml`, …), frozen before the campaign runs
(cross-referenced from `ANALYSIS_PLAN.md`). A config fully specifies a reproducible
run: geometry matrix, energies, event counts, seed base, physics overrides, and the
output `campaign_id`. The driver (`sim/campaign.py`, and per-campaign wrappers) reads
it, and each run records `seed` + `git_commit` in its `.root.meta.json`.

Schema (minimal):
```yaml
campaign_id: M1_impact          # -> results/M1_impact/
description: downstream reconstruction impact
geometry:
  topologies: [rectilinear, gyroid]
  infills: [0.30]
  cell_mm: 2.5
  depth_mm: 10.0
energies_mev: [200, 1000]
n_events: 5.0e6
seed_base: 12345
physics: {em: option4, cut_mm: 0.05, max_step_mm: 0.1}   # locked unless a systematic
outputs: results/M1_impact
```
Phase 0 used `sim/campaign.py` directly (hardcoded matrix); the YAML layer is introduced
for M1 onward.
