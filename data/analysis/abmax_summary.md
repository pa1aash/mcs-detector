# AB-MAX: diamond high-stats + Voronoi reseed (pod, 2026-07-16, Geant4 11.3.2)

40 runs at 3e7 protons, locked physics/windowing/all-order floor. Raw abmax_results.json.
Honest branches applied: diamond NOT promoted (kept excluded), Voronoi NOT restored (kept
demoted) — both now backed by real high-statistics data instead of a disclosed limitation.

## A-MAX — diamond high-stats (30x the original 1e6; spc=80 field, N_eff = analytic nxy=80)
| f | E | N_eff | Δκ4 [rad^4] (95% CI) | C·N_eff | excludes 0? |
|--|--|--:|--|--:|:--:|
| 0.20 | 200 | 12.47 | 5.01e-11 [4.44,5.57]e-11 | 0.87 | yes |
| 0.30 | 200 | 19.07 | 4.46e-11 [3.81,5.11]e-11 | 0.90 | yes |
| 0.40 | 200 | 45.03 | 1.68e-11 [0.87,2.49]e-11 | 0.70 | yes |
| 0.50 | 200 | ∞ | -5.5e-12 [-14.6,3.6]e-12 | — | no (corner) |
| 0.20 | 500 | 12.47 | 2.19e-12 [1.98,2.39]e-12 | 1.08 | yes |
| 0.30 | 500 | 19.07 | 1.74e-12 [1.51,1.98]e-12 | 1.00 | yes |
| 0.40 | 500 | 45.03 | 9.41e-13 [6.69,12.1]e-13 | 1.12 | yes |
| 0.50 | 500 | ∞ | -2.4e-13 [-5.3,0.4]e-13 | — | no (corner) |

**Verdict: RAN and RESOLVED, but kept excluded.** The 1e6->3e7 rerun retires Fable 3.2 ("a
test that couldn't fail because it couldn't run"): the signal now excludes zero for f<=0.40 at
both energies. C·N_eff brackets the derived unity at 500 MeV (1.00-1.12) but runs 10-30% low at
200 MeV in the deepest small-signal/high-N_eff regime; f=0.50 is the N_eff->inf corner (crosses
zero). Diamond stays a consistency corner, not a fit anchor — the exponent is unchanged
(periodic-only pooled -0.932 remains the headline; NOT promoted, so no exponent recompute).

## B-MAX — Voronoi +4 seeds (78µm converged field; seeds 6006+round(100f)+1000k, k=1..4)
Realization scatter across 4 seeds vs the per-run bootstrap CI half-width:

| f | E | C·N_eff (seed mean ± std) | dk4 seed-std / mean CI-hw |
|--|--|--|--:|
| 0.20 | 200 | 0.80 ± 0.23 | 2.50 |
| 0.30 | 200 | 1.06 ± 0.49 | 5.97 |
| 0.40 | 200 | 0.94 ± 0.21 | 2.87 |
| 0.50 | 200 | 0.80 ± 0.21 | 5.37 |
| 0.20 | 500 | 0.83 ± 0.22 | 2.77 |
| 0.30 | 500 | 1.16 ± 0.46 | 5.96 |
| 0.40 | 500 | 0.94 ± 0.19 | 2.66 |
| 0.50 | 500 | 0.84 ± 0.23 | 5.48 |

**Verdict: large scatter -> demotion CONFIRMED and QUANTIFIED.** The seed-to-seed realization
scatter is 2.5-6x the per-proton bootstrap CI at both energies (±20-45% in C·N_eff), so the
single-seed Voronoi Δκ4 genuinely cannot carry quantitative load. Phase-1's demotion to a
qualitative exemplar stands, now with a measured realization-uncertainty number (Fable 3.5).
