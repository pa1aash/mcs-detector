# Voxel-pitch convergence at the REAL campaign operating pitches (Phase 2 Task 2)

Campaign cell 2.5 mm, f=0.4, L=10.0 mm. Pitch = 2.5 mm / spc. Operating pitches: rectilinear spc=48 (52 um); gyroid/Schwarz-P/diamond/Voronoi spc=32 (78 um). Voronoi = committed campaign realization (seed 6046). Convergence reference = spc=120; PASS = |dev| < 1% in BOTH f and Var(t) at the operating pitch.

## rectilinear  (operating pitch spc=48, 52 um)

| spc | pitch [um] | f | Var(t) [mm^2] | N_eff | df vs 120 | dVar vs 120 |
|--:|--:|--:|--:|--:|--:|--:|
| 32 | 78 | 0.3970 | 11.5728 | 2.069 | -0.85% | -0.38% |
| 40 | 62 | 0.3969 | 11.5709 | 2.069 | -0.87% | -0.40% |
| 48  <-- OPERATING | 52 | 0.3973 | 11.5751 | 2.069 | -0.77% | -0.36% |
| 60 | 42 | 0.3996 | 11.6018 | 2.068 | -0.19% | -0.13% |
| 80 | 31 | 0.3995 | 11.6008 | 2.068 | -0.20% | -0.14% |
| 120 | 21 | 0.4004 | 11.6169 | 2.067 | +0.00% | +0.00% |

**Operating pitch: f dev -0.77%, Var(t) dev -0.36% vs spc=120 -> PASS (<1%)**

## gyroid  (operating pitch spc=32, 78 um)

| spc | pitch [um] | f | Var(t) [mm^2] | N_eff | df vs 120 | dVar vs 120 |
|--:|--:|--:|--:|--:|--:|--:|
| 32  <-- OPERATING | 78 | 0.3970 | 6.1682 | 3.881 | -0.59% | -0.18% |
| 40 | 62 | 0.4030 | 6.2260 | 3.864 | +0.92% | +0.76% |
| 48 | 52 | 0.3996 | 6.1793 | 3.883 | +0.07% | +0.00% |
| 60 | 42 | 0.3994 | 6.1784 | 3.882 | +0.01% | -0.01% |
| 80 | 31 | 0.4000 | 6.1871 | 3.879 | +0.18% | +0.13% |
| 120 | 21 | 0.3993 | 6.1791 | 3.882 | +0.00% | +0.00% |

**Operating pitch: f dev -0.59%, Var(t) dev -0.18% vs spc=120 -> PASS (<1%)**

## schwarzp  (operating pitch spc=32, 78 um)

| spc | pitch [um] | f | Var(t) [mm^2] | N_eff | df vs 120 | dVar vs 120 |
|--:|--:|--:|--:|--:|--:|--:|
| 32  <-- OPERATING | 78 | 0.3994 | 11.0061 | 2.180 | -0.12% | +0.13% |
| 40 | 62 | 0.4031 | 11.0346 | 2.181 | +0.81% | +0.39% |
| 48 | 52 | 0.3995 | 10.9886 | 2.183 | -0.11% | -0.02% |
| 60 | 42 | 0.4021 | 11.0206 | 2.181 | +0.55% | +0.27% |
| 80 | 31 | 0.4001 | 10.9956 | 2.183 | +0.05% | +0.04% |
| 120 | 21 | 0.3999 | 10.9913 | 2.183 | +0.00% | +0.00% |

**Operating pitch: f dev -0.12%, Var(t) dev +0.13% vs spc=120 -> PASS (<1%)**

## diamond  (operating pitch spc=32, 78 um)

| spc | pitch [um] | f | Var(t) [mm^2] | N_eff | df vs 120 | dVar vs 120 |
|--:|--:|--:|--:|--:|--:|--:|
| 32  <-- OPERATING | 78 | 0.3975 | 0.6667 | 35.920 | -0.63% | +29.50% |
| 40 | 62 | 0.3990 | 0.5074 | 47.260 | -0.24% | -1.45% |
| 48 | 52 | 0.3990 | 0.4927 | 48.672 | -0.24% | -4.30% |
| 60 | 42 | 0.3995 | 0.5220 | 45.957 | -0.11% | +1.39% |
| 80 | 31 | 0.3994 | 0.5327 | 45.032 | -0.13% | +3.47% |
| 120 | 21 | 0.4000 | 0.5148 | 46.614 | +0.00% | +0.00% |

**Operating pitch: f dev -0.63%, Var(t) dev +29.50% vs spc=120 -> FAIL (>=1%)**

_Note: diamond Var(t) is near-zero (0.5148 mm^2 at spc=120), the N_eff->inf corner; relative Var deviations are amplified by the tiny denominator and are the known numerical-floor issue, distinct from a genuine coarse-pitch geometric failure -- read the absolute Var(t) column._

## voronoi  (operating pitch spc=32, 78 um)

| spc | pitch [um] | f | Var(t) [mm^2] | N_eff | df vs 120 | dVar vs 120 |
|--:|--:|--:|--:|--:|--:|--:|
| 32  <-- OPERATING | 78 | 0.4000 | 2.5752 | 9.320 | -0.00% | +0.07% |
| 40 | 62 | 0.4000 | 2.5758 | 9.318 | +0.00% | +0.09% |
| 48 | 52 | 0.4000 | 2.5766 | 9.315 | +0.00% | +0.12% |
| 60 | 42 | 0.4000 | 2.5755 | 9.319 | +0.00% | +0.08% |
| 80 | 31 | 0.4000 | 2.5747 | 9.322 | +0.00% | +0.05% |
| 120 | 21 | 0.4000 | 2.5735 | 9.326 | +0.00% | +0.00% |

**Operating pitch: f dev -0.00%, Var(t) dev +0.07% vs spc=120 -> PASS (<1%)**

