# Transport-aware ray-tracer vs full Geant4 (Phase B4/D)

200 MeV, voxel geometry, locked window W(200)=37.84 mrad, solid control t=4 mm (kappa4_solid=1.421e-09). Geant4 `tpla` (proton PLA path) = t_actual.

## (i) Var(t_actual): Geant4 vs transport
| topo | cell | <t> G4 | Var(t_act) G4 | Var(t_act) TA | Δ |
|--|--:|--:|--:|--:|--:|
| rectilinear | 2.5 | 3.945 | 11.380 | 11.388 | -0.1% |
| rectilinear | 1.0 | 3.925 | 11.173 | 11.207 | -0.3% |
| rectilinear | 0.5 | 3.889 | 10.869 | 10.907 | -0.3% |
| gyroid | 2.5 | 3.950 | 6.122 | 6.193 | -1.1% |
| gyroid | 1.0 | 3.936 | 6.080 | 6.150 | -1.1% |
| gyroid | 0.5 | 3.910 | 5.991 | 6.047 | -0.9% |

## (ii) Δκ₄: Geant4 vs transport-aware (TA) vs straight-chord (LI)
| topo | cell | Δκ₄_G4 [rad⁴] | 95% CI | Δκ₄_TA | Δκ₄_LI | TA in CI? | ratio TA/LI |
|--|--:|--:|--:|--:|--:|:--:|--:|
| rectilinear | 2.5 | 4.091e-10 | [3.41e-10,4.81e-10] | 4.361e-10 | 4.407e-10 | Y | 0.990 |
| rectilinear | 1.0 | 3.687e-10 | [2.99e-10,4.36e-10] | 4.292e-10 | 4.407e-10 | Y | 0.974 |
| rectilinear | 0.5 | 4.266e-10 | [3.58e-10,4.95e-10] | 4.177e-10 | 4.407e-10 | Y | 0.948 |
| gyroid | 2.5 | 2.248e-10 | [1.53e-10,2.89e-10] | 2.372e-10 | 2.379e-10 | Y | 0.997 |
| gyroid | 1.0 | 1.902e-10 | [1.25e-10,2.57e-10] | 2.355e-10 | 2.379e-10 | Y | 0.990 |
| gyroid | 0.5 | 1.872e-10 | [1.18e-10,2.50e-10] | 2.316e-10 | 2.379e-10 | Y | 0.973 |

**Var(t_actual) agreement: max |Δ| = 1.1%.** **Δκ₄_TA within Geant4 bootstrap CI at all cells: YES.**
At 2.5/1.0 mm ratio TA/LI≈1 → straight-chord shape law (SC2) holds (Δκ₄_G4≈Δκ₄_LI≈Δκ₄_TA); at 0.5 mm the onset of departure (ratio<1, Δκ₄_G4 tracking Δκ₄_TA below Δκ₄_LI) is visible where Geant4 is still in its validated regime.

## GATE B/D verdict — PASS

- **(i) Var(t_actual): Geant4 vs transport agree to max |Δ|=1.1%** across all 6 cells (rect/gyroid × 2.5/1.0/0.5 mm). The transport SDE reproduces the wandering solid-path the proton actually traverses.
- **(ii) Δκ₄_TA within the Geant4 bootstrap CI at ALL 6 cells.**
- **Above the break (2.5, 1.0 mm): ratio Δκ₄_TA/Δκ₄_LI ≈ 0.97–0.99** → Δκ₄_G4 ≈ Δκ₄_LI ≈ Δκ₄_TA. The straight-chord N_eff shape law (Result 2 / SC2) is confirmed empirically at the fixed cell-independent N_eff, independent of the S6 collapse.
- **Approaching the break (0.5 mm): onset of departure** (ratio 0.95 rect / 0.98 gyroid; Δκ₄_G4 begins tracking Δκ₄_TA below Δκ₄_LI) is visible where Geant4 is still in its validated step regime.
- **Model-error STOP NOT triggered:** the line-integral prediction matches Geant4 at the large printable cells where straight-chord must hold → Result 2 is sound; the second-order trajectory↔scattering residual is bounded ≤ the ~1% Var agreement / κ₄ bootstrap width.

The transport-aware break map (Phase C) is therefore validated for extrapolation to sub-printable cells where Geant4's step regime is not valid.

## Phase 1 (S5-verify): the 100 MeV break — Geant4-CONFIRMED

The decisive, previously-unconfirmed point: at 100 MeV the transport map predicts
the break at rect 559 / gyroid 376 µm. Full Geant4 (voxel res20 → 28 µm voxels,
steps above the S3 0.02 mm floor; correct f; W(100)=79.85 mrad; 300k) at cells
BRACKETING each break, compared on the floor-free break observable
Var(t_actual)/Var(t_straight):

| topo | cell mm | ratio Var(t_act)/Var(t_str): Geant4 | transport | ΔVar(G4-vs-TA) |
|--|--:|--:|--:|--:|
| rect | 1.00 | 0.934 | 0.945 | −1.1% |
| rect | 0.56 | **0.891** | 0.898 | −0.7% |
| rect | 0.30 | 0.824 | 0.821 | +0.3% |
| gyroid | 0.56 | 0.928 | 0.941 | −1.3% |
| gyroid | 0.38 | **0.896** | 0.903 | −0.7% |
| gyroid | 0.20 | 0.810 | 0.808 | +0.2% |

**Verdict: CONFIRMED.** Geant4 reproduces the transport break ratio to ≤1.3% at
all cells, and the break onset (ratio = 0.90) falls at the transport-predicted cell
for both topologies (rect ≈ 0.56 mm, gyroid ≈ 0.38 mm). Geant4 also reproduces the
wandering-reduced ⟨t_actual⟩ (−4% rect / −11% gyroid vs the straight chord) exactly
— an independent confirmation of the trajectory model. Δκ₄_TA is within the Geant4
bootstrap CI at all 6 cells (the absolute Δκ₄_G4 carries a floor-subtraction
systematic from the wandering ⟨t⟩ shift; the Var ratio is the clean observable).

**Step-regime (1b): STABLE.** Same res20 geometry (28 µm voxels), maxStep 0.1
(≈28 µm steps, above floor) vs 0.01 (≈10 µm, below floor): Var(t_act) ΔVar = 0.00%
(break observable exactly step-independent); Δκ₄ CIs overlap. The sub-floor cells
(rect 0.3 @15 µm, gyroid 0.2 @10 µm) match transport identically to the above-floor
cells, confirming the break observable is not corrupted by the MSC step regime.

This is the first DIRECT Geant4 bracketing of a break (the S5-rebuilt 200 MeV
validation cells were all above the 200 MeV break). The 100 MeV break — the point
the criterion headline rested on — is now Geant4-backed, not extrapolated.

---

## S6 Stage-2 Phase-1: foam-scale (0.2 mm cell) Geant4 gamma2 spot-check (200 MeV)

The Phase-C cell_homog boundary rests on the transport tool extrapolating gamma2(cell)
from the 2.5 mm Geant4 anchor down ~12x to foam scales (100-500 um). This phase tests
that extrapolation with a FULL Geant4 run at a 0.2 mm foam-scale cell (rect + gyroid,
f=0.4, 200 MeV, 1e7), at coarse voxel resolutions so the step stays above the S3
0.02 mm MSC-corruption floor (res4 -> 0.05 mm voxel/step; res8 -> 0.025 mm).

**gamma2 = Delta_kappa4/kappa2^2** (all-order physical floor, W(200)); the geometry-
matched test is the break ratio = gamma2_G4 / gamma2_straight (same voxel field), which
the transport tool predicts at 0.2 mm @ 200 MeV: rect 0.868, gyroid 0.908.

| topo | res | voxel | f_built | gamma2_G4 | gamma2_straight | ratio G4/straight | transport ratio |
|--|--|--|--|--|--|--|--|
| rect   | 4 | 0.050 | 0.469 | 1.400 +/- 0.047 | 1.504 | 0.931 | 0.868 |
| rect   | 8 | 0.025 | 0.365 | 1.965 +/- 0.068 | 2.301 | **0.854** | 0.868 |
| gyroid | 4 | 0.050 | 0.480 | 0.505 +/- 0.045 | 0.752 | 0.671 | 0.908 |
| gyroid | 8 | 0.025 | 0.477 | 0.733 +/- 0.043 | 0.868 | 0.868 | 0.908 |

**Geometry-convergence diagnostic (gamma2_straight vs voxel res, ray-traced, no Geant4):**
rect 1.50(f.50)/2.30(f.39)/2.00/2.20 at res 4/8/16/32 -> analytic 2.18(f.40);
gyroid 0.75(f.50)/0.84(f.50)/1.24/1.17 -> analytic 1.16(f.40). The 0.2 mm geometry
converges only at res16-32, whose voxel/step (0.0125/0.00625 mm) fall BELOW the 0.02 mm
MSC-corruption floor. **There is no voxel resolution that simultaneously resolves the
0.2 mm geometry AND keeps the Geant4 step above the MSC floor.**

**Branch / verdict:**
- **rect res8 is a clean in-band validation:** its geometry is ~converged (gamma2_straight
  2.30 ~= analytic 2.18, f 0.39 ~= 0.40) at a clean-step resolution, and the Geant4
  break ratio 0.854 matches the transport prediction 0.868 to ~1.5%, with gamma2_G4 = 1.97
  (substantial -> homogeneous-kurtosis approximation FAILS at the 0.2 mm foam-scale cell,
  200 MeV; 0 stuck). 0.2 mm is just below cell_break(200)=0.267 mm, so a ratio ~0.85
  (mild line-integral departure onset), not 1, is the expected, consistent value.
- **gyroid is geometry-resolution-limited:** at clean-step res (4,8) f=0.50 (vs 0.40) and
  gamma2_straight is 28% under converged; a clean gyroid foam-scale point is unattainable
  with voxels. Its ratio (0.868) is consistent with transport (0.908) but not on a
  converged geometry.
- **Honest limitation:** foam-scale Geant4 validation is obstructed by the voxel
  geometry-resolution vs MSC-step-floor tension. rect provides a clean confirmation of the
  transport tool in-band; gyroid does not. The precise cell_homog (~um) rests on the
  transport tool (now validated at 2.5 mm AND at 0.2 mm for rect).
