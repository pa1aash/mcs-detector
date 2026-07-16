# Width-envelope derivation (uncorrected nominal-control envelope = 7.9%)

Derivation of the 7.9% figure quoted in Section 5.2 (`width_invariance`, ledger row
`width_invariance`). The uncorrected width envelope against the *designed* mean fL is
set by the largest as-built mean-material mismatch <t>/(fL) across the campaign, which
the width ratio tracks identically (geometry-level identity, Section 5.2).

Source: `data/analysis/pt_hist.json` key `mean_t_over_fL` (as-built <t>/(fL) per run,
straight-chord line integral through the committed campaign voxel fields).

**max <t>/(fL) = 1.07882  (camp_voronoi_f30_E500)  =>  max deviation = 7.88% ~ 7.9%.**

Top-5 as-built mean-material mismatches:

| run | <t>/(fL) | deviation |
|--|--:|--:|
| camp_voronoi_f30_E500 | 1.07882 | +7.88% |
| camp_voronoi_f30_E200 | 1.07872 | +7.87% |
| camp_voronoi_f20_E200 | 1.06378 | +6.38% |
| camp_voronoi_f20_E500 | 1.06365 | +6.36% |
| camp_voronoi_f40_E200 | 1.02616 | +2.62% |

Recomputed 2026-07-16 (Phase-1 hygiene, audit ID WIDTH-79). The value was already
stored inside `pt_hist.json` `mean_t_over_fL`; this note makes the max/derivation
explicit and is the provenance target pinned by claims_ledger.csv row `width_invariance`.
