> **RE-ANCHOR NOTE (2026-07-16, referee M2).** The headline stochastic case is now the
> FITTED Voronoi points ($f=0.40$--$0.50$), not the excluded floor-dominated $f=0.30$.
> Beam-time from the committed `e2_results_combined.json` via `N@10% = N_have*(cif/0.10)^2`:
> voronoi f40 = 5.0e7/6.1e7, f50 = 7.1e7/7.6e7 protons (200/500 MeV); cost ratio ~2-3e5x
> (using the signal-independent width-channel `cif_theta0`). `e3_affordability.json` and
> `fig_affordability` are re-anchored to f40; the voronoi rows of the RATIO TABLE below still
> show the old f30 `.root`-derived `cif_theta0`/`cif_k4raw` and need a compute-box rerun of
> `analysis/e3_affordability.py` (updated `REP`) to refresh them. The beam-time numbers and
> the figure are correct as shown.

# e3_affordability.md -- S7 / E3: affordability of the shape channel

**Calibration-anchored** (no folklore; the 470M-track material-budget-imaging figure is NOT used). All numbers scale the MEASURED deconvolved-Delta_kappa4 fractional CI on the on-disk campaign runs by the moment-estimator law (SE ~ 1/sqrt(N), so cif ~ 1/sqrt(N)): **N(target p) = N_have x (cif_have/p)^2**.

## Flux anchor (verified, in-regime)

The LLUMC/UCSC **phase-II proton-CT scanner** sustains **>1e6 individually-measured protons/s** (1 MHz) at ~200 MeV through a silicon-strip tracker, completing a full CT scan in <=6 min of beam time [Johnson et al., *Phys. Procedia* **90** (2017) 209, doi:10.1016/j.phpro.2017.09.060]. This is the same metrology class (single-proton tracking telescope) and energy regime as the proposed shape channel, so we quote beam-time at **1 MHz** as a conservative, citable single-proton rate.

## Protons + beam-time to resolve the geometry Delta_kappa4 (representative configs)

| E [MeV] | config | N_eff | measured N (cif) | N@30% | N@20% | N@10% | t@30% | t@20% | t@10% (1 MHz) |
|--:|--|--:|--|--:|--:|--:|--:|--:|--:|
| 200 | rectilinear f40 | 2.1 | 3.0e+06 (10%) | 3.28e+05 | 7.37e+05 | 2.95e+06 | 0.3s | 0.7s | 2.9s |
| 200 | gyroid f40 | 3.9 | 6.0e+06 (12%) | 1.04e+06 | 2.33e+06 | 9.33e+06 | 1.0s | 2.3s | 9.3s |
| 200 | voronoi f30 | 13.5 | 3.0e+07 (24%) | 1.86e+07 | 4.19e+07 | 1.68e+08 | 18.6s | 41.9s | 2.8min |
| 500 | rectilinear f40 | 2.1 | 3.0e+06 (9%) | 2.66e+05 | 5.99e+05 | 2.39e+06 | 0.3s | 0.6s | 2.4s |
| 500 | gyroid f40 | 3.9 | 6.0e+06 (13%) | 1.07e+06 | 2.41e+06 | 9.64e+06 | 1.1s | 2.4s | 9.6s |
| 500 | voronoi f30 | 13.5 | 3.0e+07 (27%) | 2.38e+07 | 5.37e+07 | 2.15e+08 | 23.8s | 53.7s | 3.6min |

**Read-out:** a printable/foam-relevant config (rect/gyroid, N_eff~2-4) is resolvable to 20% in **seconds** and to 10% in **<10 s** of 1-MHz beam; the hardest representative case (high-N_eff stochastic voronoi, a small signal on a large floor) needs **~minutes** to 10%. The shape channel is **feasible at a proton-CT-class facility with seconds-to-minutes of beam**, not 'cheap' and not 10^8-track-expensive.

## Cost RATIO vs the standard WIDTH (Highland) channel

To reach the SAME fractional precision on its respective observable, the shape channel costs more protons than the width channel by **N_shape/N_width = (cif_dk4 / cif_theta0)^2** (both CIs ~ 1/sqrt(N)). Decomposed into the moment-order penalty (raw 4th vs 2nd moment) x the small-signal/deconvolution penalty (Delta_kappa4 is a small geometry excess on a large intrinsic floor):

| E [MeV] | config | cif theta0 | cif kappa4_raw | cif Delta_kappa4 | moment-order x | deconv/small-signal x | **total shape/width x** |
|--:|--|--:|--:|--:|--:|--:|--:|
| 200 | rectilinear f40 | 0.12% | 1.4% | 10% | 28x | 48x | **7055x** |
| 200 | gyroid f40 | 0.08% | 1.1% | 12% | 40x | 126x | **26187x** |
| 200 | voronoi f30 | 0.03% | 0.6% | 24% | 67x | 1780x | **743872x** |
| 500 | rectilinear f40 | 0.12% | 1.4% | 9% | 33x | 39x | **5933x** |
| 500 | gyroid f40 | 0.08% | 1.2% | 13% | 45x | 112x | **28020x** |
| 500 | voronoi f30 | 0.03% | 0.5% | 27% | 47x | 2704x | **1036408x** |

**Honest statement.** The width (Highland) channel measures a large primary observable (theta0) to sub-percent precision at ~1e3-1e4 protons; the shape channel measures a small deconvolved 4th-cumulant difference and so costs **~10^3-10^4x more protons for the printable-relevant configs (rect/gyroid), rising to ~10^6x for the hardest high-N_eff voronoi**, at matched fractional precision -- both because it is a 4th moment (the moment-order penalty, ~30-70x) and because Delta_kappa4 is a small excess on a large intrinsic floor (the deconvolution/small-signal penalty, ~40-2700x). That multiplier is the real price of the shape channel. It nonetheless remains affordable in ABSOLUTE terms (seconds-to-minutes at 1 MHz) because the width channel is so cheap to begin with: ~1e6-1e8 protons is seconds-to-minutes of proton-CT-class beam.

## Provenance

- N and cif are MEASURED on the on-disk locked campaign runs (`e2_results_200.json`, `e2_results.json`); scaling is the 1/sqrt(N) moment law. Cross-check: the Phase-0c calibration extrapolated voronoi f30 @500 to ~1.4e7 protons for 30% -- this script's calibration-anchored scaling reproduces the same ~1e7 scale.
- Flux: Johnson et al. 2017 (verified via the published abstract: '>1e6 protons individually measured per second', full scan '<=6 min').

