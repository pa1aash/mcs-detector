# Break threshold + nozzle sensitivity (Phase 2)

## c(p, topology, threshold) [µm] — how cell_break moves with the threshold

**rectilinear**

| p [MeV] | 5% | 10% | 20% | 30% |
|--:|--:|--:|--:|--:|
| 100 | 1174 | 559 | 271 | 164 |
| 200 | 504 | 260 | 131 | 79 |
| 500 | 212 | 111 | 56 | 33 |
| 1000 | 128 | 62 | 30 | <30 |

**gyroid**

| p [MeV] | 5% | 10% | 20% | 30% |
|--:|--:|--:|--:|--:|
| 100 | 630 | 376 | 187 | 109 |
| 200 | 297 | 183 | 91 | 53 |
| 500 | 131 | 76 | 39 | <30 |
| 1000 | 74 | 44 | <30 | <30 |

## cell_break vs nozzle floor — does the break land in printable range?

Entry = ✓ if cell_break ≥ nozzle (the break sits at/above a printable feature, i.e. the criterion bites for that printer); — otherwise (break is finer than the printer can make → design-rule for that printer).

| topo | p | c_break(10%) µm | FDM 0.25 | FDM 0.4 | FDM 0.5 | FDM 0.8 | SLA 0.05 |
|--|--:|--:|:--:|:--:|:--:|:--:|:--:|
| rectilinear | 100 | 559 | ✓ | ✓ | ✓ | — | ✓ |
| rectilinear | 200 | 260 | ✓ | — | — | — | ✓ |
| rectilinear | 500 | 111 | — | — | — | — | ✓ |
| rectilinear | 1000 | 62 | — | — | — | — | ✓ |
| gyroid | 100 | 376 | ✓ | — | — | — | ✓ |
| gyroid | 200 | 183 | — | — | — | — | ✓ |
| gyroid | 500 | 76 | — | — | — | — | ✓ |
| gyroid | 1000 | 44 | — | — | — | — | — |

## Robustness band (5–20% threshold) for the headline rect@100 MeV

- rectilinear @100 MeV: c_break = 1174 / 559 / 271 µm (5/10/20%). FDM-0.5 overlap (c_break ≥ 0.5 mm): Y/Y/N.
- gyroid @100 MeV: c_break = 630 / 376 / 187 µm (5/10/20%). FDM-0.5 overlap (c_break ≥ 0.5 mm): Y/N/N.

## [S5-AUDIT/REFRAME NOTE — see DIRECTION_AUDIT.md]

The "criterion bites if cell_break ≥ nozzle" row above frames printability as the
relevance test. Printability stays **dropped** as a load-bearing claim (genuinely
threshold-fragile: rect@100/FDM-0.5 holds at 5–10% but fails at 20%, one topology
at one momentum). The physics result is the momentum-gated **foam-scale** failure
boundary (homogeneous fails for pore d above p\*(d); see `foam_overlap.md` and
`DIRECTION_AUDIT.md`), which is independent of any printer. Direction note: the
homogeneous approximation fails for cells COARSER than cell_break, not finer.
