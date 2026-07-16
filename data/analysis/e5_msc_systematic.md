# e5_msc_systematic.md -- S7 / E5 MSC-model systematic (opt3)

Alternative MSC model: **G4EmStandardPhysics_option3 (UrbanMsc-based)**, vs the locked `option4` (WentzelVI multiple scattering + single Coulomb). The geometry (voxel fields, ray-traced N_eff) is byte-identical to the locked campaign, so the shift is purely the MSC model. This alternative is a LARGE perturbation (a deliberately conservative stress test): see the re-derived scattering power below.

## E5b inputs: scattering power a_eff and floor under the alt MSC

| E [MeV] | a_eff locked | a_eff alt | a_alt/a_locked | Highland resid locked | Highland resid alt |
|--:|--:|--:|--:|--:|--:|
| 200 | 3.9127e-06 | 3.2736e-06 | **0.837** | -3.07% | -9.06% |
| 500 | 7.2061e-07 | 6.0520e-07 | **0.840** | -3.65% | -9.03% |

## E5a -- THE COLLAPSE (clean regime): exponent + Delta_kappa4 band

| E [MeV] | locked OLS slope (CI) | alt OLS slope (CI) | alt in locked CI? | alt brackets -1? |
|--:|--|--|:--:|:--:|
| 200 | -0.951 [-1.040,-0.866] | -1.030 [-1.422,-0.624] | **True** | True |
| 500 | -1.001 [-1.095,-0.908] | -1.182 [-1.416,-0.946] | **False** | True |

NOTE: 200 MeV alt = PERIODIC-only (no voronoi run at 200) -> -1.030 is the clean periodic
slope. 500 MeV alt = ALL-INCLUSIVE (incl. 2 voronoi) -> -1.182; the **PERIODIC-only 500 MeV
alt slope = -1.002** (inside the locked CI). The headline is the periodic slope (-1.00 at both
energies, MSC-robust); the all-inclusive -1.18 is voronoi-floor-driven (see verdict + diagnostic).

**Per-config Delta_kappa4 fractional shift (alt vs locked), within locked bootstrap CI?**

| E | config | N_eff | dk4 locked | dk4 alt | shift | locked CI | within CI? |
|--:|--|--:|--:|--:|--:|--:|:--:|
| 200 | rectilinear_f40 | 2.1 | 5.117e-10 | 2.608e-10 | -49% | 10% | False |
| 200 | gyroid_f20 | 5.5 | 1.310e-10 | 6.816e-11 | -48% | 19% | False |
| 200 | gyroid_f40 | 3.9 | 3.010e-10 | 1.654e-10 | -45% | 12% | False |
| 200 | schwarzp_f40 | 2.2 | 5.310e-10 | 3.127e-10 | -41% | 10% | False |
| 500 | rectilinear_f40 | 2.1 | 1.741e-11 | 8.933e-12 | -49% | 9% | False |
| 500 | gyroid_f20 | 5.5 | 4.286e-12 | 2.415e-12 | -44% | 23% | False |
| 500 | gyroid_f40 | 3.9 | 9.473e-12 | 5.617e-12 | -41% | 13% | False |
| 500 | schwarzp_f40 | 2.2 | 1.791e-11 | 1.067e-11 | -40% | 10% | False |
| 500 | voronoi_f40 | 9.3 | 3.953e-12 | 1.620e-12 | -59% | 15% | False |
| 500 | voronoi_f50 | 12.6 | 4.159e-12 | 1.283e-12 | -69% | 15% | False |

## E5b -- THE BOUNDARY: cell_homog + foam-scale gamma2 under the alt MSC

gamma2 = 3 Var(t_act)/<t>^2 carries the scattering power `a` ONLY through the wandering, so the FOAM-scale gamma2 (cell >> wandering scale) is nearly a-independent, while cell_homog (the fine gamma2->0.1 crossing) moves with `a`.

| E | topo | gamma2@0.2mm locked | gamma2@0.2mm alt | cell_homog locked | cell_homog alt |
|--:|--|--:|--:|--:|--:|
| 200 | rectilinear | 2.13 | 2.13 | 0.06 um | 0.08 um |
| 200 | gyroid | 1.19 | 1.18 | 0.13 um | 0.16 um |
| 500 | rectilinear | 2.12 | 2.12 | 0.00 um | 0.00 um |
| 500 | gyroid | 1.15 | 1.15 | 0.66 um | 0.68 um |

**Foam-scale Geant4 (rect res8, 0.2 mm, 200 MeV) MEASURED under the alt MSC:** gamma2 =
**0.90** (f_built 0.37, geometry = 17% of the raw kappa4 -- a moderate, reliable fraction,
NOT the voronoi small-difference fragility), vs the locked-physics Geant4 gamma2 = 1.97.
**The qualitative broad failure SURVIVES** -- gamma2 = 0.90 is still a large excess kurtosis
(the homogeneous/Highland model predicts gamma2 = 0), Geant4-anchored under BOTH physics
lists, independent of the transport tool. **But the MAGNITUDE carries a large MSC band:
gamma2 ~ 0.9-2.0 (factor ~2)** -- under Urban's thinner tails both the raw structured and the
intrinsic-floor kurtosis shrink (raw gamma2 8.17->5.37, floor 6.20->4.47), and the geometry
excess (their difference) roughly halves.

## Verdict

- **E5a (collapse, clean regime) -- GATE PASS (ruled):** the N_eff^-1 EXPONENT is MSC-robust:
  the **PERIODIC-only OLS slope is -1.00 under both physics lists** (alt: 200 MeV -1.030, 500 MeV
  -1.002; locked periodic -0.987) -- the load-bearing parameter-free result holds. The
  all-inclusive 500 MeV slope steepens to **-1.18** (outside the locked CI) SOLELY via the two
  stochastic-voronoi points (k_eff 1.0->0.5); the raw-kappa4 diagnostic (`e5a_rawk4_check.md`)
  shows voronoi geometry is only ~6% of its raw kappa4 (floor-dominated) and rests on a +/-0.006
  struct-floor gap at floor-model-accuracy level -> **voronoi DEMOTED to qualitative-only**
  (its 2nd documented floor fragility; user-ruled). Exponent MSC band: **-1.00 (periodic,
  load-bearing) -> -1.18 (all-inclusive, voronoi-floor-attributed)**. Prefactor band **k =
  0.80 (alt) <- 0.98 (locked)**, uniform across N_eff -> confirms the slope/prefactor split.
  Absolute Delta_kappa4 band **-40 to -70%** (a_eff -16%; expected a^2 rescaling).

- **E5b (boundary):** TWO distinct, separately-stated results --
  1. **The QUALITATIVE foam-shape failure is MSC-robust (Geant4-anchored):** gamma2 stays O(1)
     and far from the homogeneous gamma2 = 0 under BOTH physics lists (1.97 WentzelVI, 0.90
     Urban). The broad failure does NOT vanish under the MSC change; it does NOT depend on the
     transport tool.
  2. **The PRECISE magnitude / cell_homog carries an MSC systematic that the fine-cell regime
     forbids bounding cleanly (stated limitation):** the foam gamma2 MAGNITUDE has a ~2x MSC
     band (0.9-2.0). The transport tool predicts gamma2@0.2mm is a-independent (2.13 both
     physics), but the Geant4 measurement shows it IS MSC-sensitive -- i.e. the scale-mixture
     closure that the transport tool rests on degrades toward fine cells under Urban (k_eff ~
     0.8 at the 2.5 mm printable cell -> ~0.4 at the 0.2 mm foam cell). So the precise
     cell_homog (~um) CANNOT be cleanly MSC-bounded; it is reported as ~um (sub-foam under both
     physics, so the conclusion -- foams fail the homogeneous kurtosis at all practical momenta
     -- is qualitatively unchanged) WITH this explicit transport-tool MSC limitation, not as a
     clean number.

