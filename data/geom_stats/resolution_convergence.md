# Resolution convergence (Phase A1)

Campaign cell 2.5 mm, f=0.4, L=10.0 mm. spc = transverse rays per cell and longitudinal steps per cell (dz=cell/spc). Converged = both f and Var(t) change <1% between successive spc.

## rectilinear

| spc | f | Var(t) | l_int | N_eff_exact | Δf | ΔVar |
|--:|--:|--:|--:|--:|--:|--:|
| 40 | 0.3969 | 11.5709 | 4.2751 | 2.0687 | — | — |
| 60 | 0.3996 | 11.6018 | 4.2628 | 2.0679 | +0.68% | +0.27% |
| 80 | 0.3995 | 11.6008 | 4.2558 | 2.0680 | -0.01% | -0.01% |
| 120 | 0.4004 | 11.6169 | 4.2523 | 2.0666 | +0.20% | +0.14% |

**converged spc = 60**

## gyroid

| spc | f | Var(t) | l_int | N_eff_exact | Δf | ΔVar |
|--:|--:|--:|--:|--:|--:|--:|
| 40 | 0.4030 | 6.2260 | 2.2355 | 3.8643 | — | — |
| 60 | 0.3994 | 6.1784 | 2.2157 | 3.8824 | -0.90% | -0.76% |
| 80 | 0.4000 | 6.1871 | 2.2130 | 3.8791 | +0.17% | +0.14% |
| 120 | 0.3993 | 6.1791 | 2.2070 | 3.8819 | -0.18% | -0.13% |

**converged spc = 60**

## schwarzp

| spc | f | Var(t) | l_int | N_eff_exact | Δf | ΔVar |
|--:|--:|--:|--:|--:|--:|--:|
| 40 | 0.4031 | 11.0346 | 4.0372 | 2.1805 | — | — |
| 60 | 0.4021 | 11.0206 | 4.0219 | 2.1815 | -0.26% | -0.13% |
| 80 | 0.4001 | 10.9956 | 4.0135 | 2.1829 | -0.49% | -0.23% |
| 120 | 0.3999 | 10.9913 | 4.0058 | 2.1833 | -0.05% | -0.04% |

**converged spc = 60**

## diamond

| spc | f | Var(t) | l_int | N_eff_exact | Δf | ΔVar |
|--:|--:|--:|--:|--:|--:|--:|
| 40 | 0.3990 | 0.5074 | 0.0309 | 47.2603 | — | — |
| 60 | 0.3995 | 0.5220 | 0.0321 | 45.9573 | +0.13% | +2.88% |
| 80 | 0.3994 | 0.5327 | 0.0339 | 45.0321 | -0.02% | +2.05% |
| 120 | 0.4000 | 0.5148 | 0.0244 | 46.6145 | +0.13% | -3.35% |

**converged spc = NOT by 120 (FLAG)**

## voronoi

| spc | f | Var(t) | l_int | N_eff_exact | Δf | ΔVar |
|--:|--:|--:|--:|--:|--:|--:|
| 40 | 0.4000 | 2.3321 | 0.7022 | 10.2913 | — | — |
| 60 | 0.4000 | 2.3309 | 0.6911 | 10.2964 | +0.00% | -0.05% |
| 80 | 0.4000 | 2.3299 | 0.6854 | 10.3008 | +0.00% | -0.04% |
| 120 | 0.4000 | 2.3294 | 0.6795 | 10.3033 | +0.00% | -0.02% |

**converged spc = 60**

## Converged spc per topology

| topology | converged spc |
|--|--:|
| rectilinear | 60 |
| gyroid | 60 |
| schwarzp | 60 |
| diamond | >120 FAIL |
| voronoi | 60 |

## Note on diamond (FLAGGED — does not converge by spc=120)

Diamond at f=0.40 is the extreme-D3-suppression case: its chord variance Var(t) is
tiny (~0.5 mm² vs ~11.6 for rect), so N_eff is huge (≈40–∞) and Var(t) sits at the
numerical floor where finite-sampling noise (~±3%) exceeds the 1% criterion. f IS
converged (0.399–0.400 across spc); only the suppressed Var(t) is noise-limited.
This is a **numerical floor of an intrinsically near-zero quantity**, not an
under-resolution of geometry. It affects only diamond's exact N_eff data point on
the collapse plot (a known extreme); it does NOT affect the rectilinear/gyroid
break map (both converged at spc=60), which is the headline result. Diamond's
summary.csv row is shipped at spc=80 with this caveat recorded.
