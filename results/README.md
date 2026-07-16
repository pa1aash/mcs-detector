# results/ — per-campaign outputs + provenance (work order Ground Rule 9)

Each campaign writes to `results/<campaign_id>/` with a `MANIFEST.json`:

```json
{
  "campaign_id": "M1_impact",
  "config": "configs/M1.yaml",
  "config_sha256": "<hash of the frozen YAML>",
  "git_commit": "<HEAD at run time>",
  "geant4_version": "11.3.2",
  "n_events_total": 40000000,
  "wall_time_s": 3600,
  "seeds": {"<tag>": <seed>, ...},
  "runs": ["data/runs/<tag>.root", ...],
  "created_utc": "<timestamp>"
}
```

Event-level `.root` files stay under `data/runs/` (gitignored, regenerable); `results/`
holds the derived tables/figures/manifests that ARE committed. `results/convergence/`
holds the step×0.5 / cut×0.5 / voxel×2 convergence re-runs required at every new regime
(new energy, density, tilt, fine cell) — Ground Rule 15.

Phase 0 (reduced-stat reproduction) predates this layer; it ran via `sim/campaign.py`
into `data/runs/` and is logged in `RUNLOG.csv`. M1 onward uses `results/<id>/`.
