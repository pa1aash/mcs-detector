# Foam-overlap momentum (Phase 3)

Carbon-foam pore band ~100-500 um (ATLAS/ALICE carbon-foam staves; refs atlasibl2018, feher2022). cell_break(p) at 10% ENTERS the band when 100 um <= cell_break <= 500 um. cell_break decreases with momentum.

| topology | foam-band momentum window | 100 | 200 | 500 | 1000 MeV |
|--|--|--:|--:|--:|--:|
| rectilinear | ~111–567 MeV | 559  | 260✓ | 111✓ | 62  |
| gyroid | ~<100–376 MeV | 376✓ | 183✓ | 76  | 44  |

**Both topologies' validity boundary sweeps THROUGH the ~100-500 um carbon-foam pore band across feasible momenta (~100-550 MeV)** -- the robust anchor of the criterion, independent of the single rect@100 MeV/FDM point.

---

## [S5-AUDIT/REFRAME CORRECTION — direction fixed; see DIRECTION_AUDIT.md]

The "ENTERS / sweeps THROUGH the band" framing above is directionally imprecise.
The correct relevance condition (confirmed from data, commit 53b3106): the
homogeneous-density approximation FAILS for a foam of pore size **d** when
**d > cell_break(p)** — i.e. the pore is COARSER than the wandering-set break, so
the structure is NOT averaged out. Since cell_break decreases with momentum, this
holds for **p ABOVE p\*(d)** (not "across feasible momenta" symmetrically, and not
"low-momentum"):

| threshold | rect p\*(500/250/100 µm) | gyroid p\*(500/250/100 µm) |
|--|--|--|
| 5%  | 202 / 420 / >1000 MeV | 124 / 243 / 694 MeV |
| 10% | 111 / 209 / 567 MeV | <100 / 148 / 376 MeV |
| 20% | <100 / 108 / 268 MeV | <100 / <100 / 183 MeV |

Reading: a foam of pore d is non-Gaussian (homogeneous fails) for **p > p\*(d)**;
coarse (~500 µm) foams fail from ~111 MeV, fine (~100 µm) from ~567 MeV (rect,
10%). At therapy energies the finest foams remain safely averaged. The failure
regime is robust to the threshold choice (onset shifts, stays feasible) but is
**momentum-gated, not universal**.
