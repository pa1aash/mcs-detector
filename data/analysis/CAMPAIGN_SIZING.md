# CAMPAIGN_SIZING.md — S6 GATE-0 report (STOP point before full-budget spend)

Phase-0 pipeline proof + cost calibration. The full campaign (Phase 1) is sized
here and must be run on a cluster — see "Allocation" below. **No full-budget spend
yet, per GATE 0.**

## Pipeline proof (Phase 0a/0b)

- **Resumable driver** `sim/campaign.py`: config matrix → voxel build → macro →
  `mcs_sim` (locked physics) → ROOT + provenance sidecar; skips any config whose
  meta `n_events ≥` target; per-run TSV log; `--only/--list/--lattice-only/
  --baselines-only`; tiered via `--nevt`. Resume verified (re-run skips completed).
- **Dry run, all 66 configs @ 2e4, 655 s, exit 0.** Full chain proven end-to-end
  (run → collect → window → κ4-deconvolve → N_eff fit → QC → figure).
- **0 stuck tracks across ALL 5 topologies** (incl. the diamond/voronoi
  navigator-stress cases) — voxel geometry is navigation-clean.
- At 2e4 all 60 lattice configs are correctly flagged **thin** (deconvolved Δκ4
  unresolved) and the slope fit is suppressed — confirming dry-run-then-scale.

## Timing calibration (Phase 0c) — Voronoi 0.30 @ 500 MeV @ 1e6

| quantity | value |
|--|--|
| wall time (6 threads) | **190 s** @ 1e6 primaries |
| peak resident memory | **0.69 GB** (res=32 voxel, 1.2 Mvox) |
| stuck tracks | 0 |
| raw κ4 fractional CI @ 1e6 | 6.4% (29.7% @ 50k, 14.2% @ 200k, 8.9% @ 500k) |

Per-config wall scales ~linearly with primaries and with voxel count:
- **res=32** topologies (gyroid, schwarzp, diamond, voronoi; 1.2 Mvox):
  ~190 s @ 1e6 → **~32 min @ 1e7**.
- **res=48** rectilinear (4.0 Mvox, needed for thin-strut f-convergence):
  ~3.4× → **~11 min @ 1e6, ~108 min @ 1e7**.

## Two correctness findings the calibration caught (fixed before Phase 1)

1. **Areal-density-matched baseline (critical).** Δκ4 = κ4(struct) − κ4(solid)
   requires the solid control matched to the lattice's mean traversed material
   t = f·L, NOT a single t=4mm slab. With the wrong t=4mm baseline, voronoi f30 @
   500 gave **Δκ4 = −8.2e-12 (unphysical negative)**; with the matched t=3mm
   baseline, **Δκ4 = +2.4e-12 (correct positive)**. *Fixed:* driver now emits solid
   controls at t={2,3,4,5} mm per energy; E2 subtracts the f-matched one.
2. **κ4-stability threshold.** A raw-κ4 CI < 100% is not enough to fit a slope; the
   QC now flags a config thin unless its **deconvolved Δκ4** CI < 30%.

## Statistics budget (the deconvolved-Δκ4 CI, with the corrected baseline)

The trend quantity is the deconvolved Δκ4, whose fractional CI is set by the SIGNAL
size (small at high N_eff). The calibration config is a HARD case:

| config | N_eff | Δκ4 @1e6 | deconv. CI @1e6 | for 30% CI | for 20% CI |
|--|--:|--:|--:|--:|--:|
| voronoi f30 @500 (high-N_eff, small signal) | ~26 | +2.4e-12 | 113% | ~1.4e7 | ~3.2e7 |
| low-N_eff configs (rect/coarse, N_eff~2) | ~2 | ~5–10× larger | (extrapolated) ~15–35% | ~1e6 | ~2e6 |

⇒ **TIERED budget (confirmed necessary):**
- **1e6 baseline** resolves low-to-mid-N_eff configs (large Δκ4) to ≤30%.
- **1e7 escalation** for high-N_eff small-signal configs (the Voronoi f≥0.30 /
  fine end + the collapse anchors). Even 1e7 gives only ~30% CI for the very
  highest N_eff — an honest limitation: the collapse slope will be anchored by the
  well-resolved low-to-mid N_eff configs, with the high-N_eff points carrying wide
  bars (still consistent with the −1 law, not driving it).

## Allocation (the resolved [FILL]s)

60 lattice + 12 solid controls + 3 empty = **75 configs**. Tiered:
- ~48 configs @ 1e6 (the resolvable ones) + ~24 @ 1e7 (high-N_eff + anchors).
- Wall estimate (6-thread-equiv): 48×(0.05–0.18 h) + 24×(0.5–1.8 h) ≈ **15–45
  core-hours-equivalent at 6 threads ≈ 90–270 core-hours**. Dominated by the
  rect res=48 @1e7 configs (~1.8 h each).
- **Phase-1 cost optimization (recommended):** rect needs res=48 only at f≤0.30
  (thin struts); f≥0.40 converges at res=32. Applying res=32 to rect f≥0.40 cuts
  the rect @1e7 cost ~halve → total ≈ **70–200 core-hours**.
- Peak memory ≤ ~1 GB/process (res=32) and ~2.5 GB (res=48 @4Mvox) → trivially
  fits standard cluster nodes; parallelize across configs (embarrassingly).

RESOLVED: [dry-run-vs-full] → dry-run-then-scale (done). [compute] → ~70–270
core-hours on a cluster (above). [statistics] → tiered 1e6/1e7 (above).

## STOP — GATE 0

Pipeline proven, cost sized, two correctness bugs fixed. **Provision the cluster
(~70–270 core-hours, ≤2.5 GB/proc, embarrassingly parallel), then resume:** the
driver skips the dry-run configs and runs the full tiered campaign. Phases 1–5
(campaign → baselines → QC → E2 → homogeneous boundary) follow with the spine
already locked (corrected-direction HEADLINE; Voronoi-foam as the N_eff workhorse).
