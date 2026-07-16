# homog_boundary.md — S6 Stage-2 Phase C: the homogeneous-approximation SHAPE boundary

Replaces the `cell_break` proxy for the foam-failure claim with the direct,
gamma2-based homogeneous-approximation boundary, and validates it with full Geant4
in the foam band. Tooling: `analysis/homog_boundary.py` (transport gamma2 curves +
numerical-floor check), `sim/run_foam_spotcheck.py` + `geom/make_foam_voxels.py`
(0.2 mm Geant4), `figs/src/fig_homog_boundary.py`. Data: `homog_boundary.json`,
`foam_spotcheck.json`, `transport_vs_g4.md`.

## The distinction
`cell_break(p)` (Var(t_act)/Var(t_str) = 0.90) is the **line-integral** boundary —
where the straight-chord approximation to Var(t) starts to fail. It is NOT where the
homogeneous (Gaussian) model is recovered: **at cell_break the excess kurtosis gamma2
is still ~90% of its geometry maximum.** The homogeneous-SHAPE boundary is where
gamma2 -> 0, which is far finer.

gamma2 = 3 Var(t_act)/<t>^2 (the geometry-induced excess kurtosis; momentum-dependent
only through the wandering Var(t_act)).

## C1 — gamma2 MEASURED at the printable 2.5 mm cell (Geant4)
At the printable cell, full Geant4 gives the geometry excess kurtosis directly
(gamma2 = Delta_kappa4/kappa2^2, all-order physical floor, W(E)):

| topology | gamma2 @2.5mm, 200 MeV | @500 MeV | transport tool |
|--|--:|--:|--:|
| rectilinear | 2.22 | 2.19 | 2.12 (ratio 0.97-0.99) |
| gyroid      | 1.28 | 1.16 | 1.15 (ratio 0.91-1.00) |

The homogeneous (Gaussian) model predicts gamma2 = 0; the measured gamma2 ~ 1.2-2.2 is
a large, momentum-independent failure of the homogeneous approximation at the printable
cell, MEASURED at two energies, and matching the transport tool to <=9%.

## C2 — cell_homog(p) via gamma2, vs cell_break(p)
gamma2(cell) stays pinned at its geometry maximum for every cell coarser than the
lateral-wandering scale y_rms (tens of um), and only then falls to zero. At 200 MeV,
e.g. rectilinear: gamma2 = 2.18 @2.5mm, 2.16 @0.1mm, 1.59 @0.01mm. Defining
cell_homog(p) where gamma2 crosses a stated threshold (gamma2 = 0.1) puts it at the
**~um scale — two-to-three orders of magnitude below cell_break** (rect cell_break =
557/267/118/73 um at 100/200/500/1000 MeV) and below every printable and foam feature
size. (At gamma2 = 0.1 the crossing is sub-um and only weakly bracketed by the feasible
transport sweep (>= 5 um); cell_homog is reported as ~um, the conclusion robust to the
exact threshold given the gamma2(cell) shape.)

## C3 — numerical-floor check: the broad failure is PHYSICAL
The fine-cell gamma2 is stable against both the transport step size and the sample
count, so the broad failure is physical, not a tool floor (rect @ 100 MeV):

| cell | gamma2 spread, spc 8->24 | gamma2 spread, N 8k->16k |
|--|--:|--:|
| 0.05 mm | 1.94-1.98 (<3%) | <0.1% |
| 0.02 mm | 1.56-1.61 (<3%) | <1% |
| 0.01 mm | 1.38-1.43 (<4%) | <0.5% |

## Phase 1 — foam-scale Geant4 validation (0.2 mm, 200 MeV)
The C2 boundary rests on extrapolating the transport gamma2 ~12x in cell size below the
2.5 mm Geant4 anchor. A full Geant4 run at a 0.2 mm foam-scale cell tests it, at coarse
voxel resolution so the step stays above the 0.02 mm MSC-corruption floor. Geometry-
matched test: ratio gamma2_G4/gamma2_straight (same voxel field) -> the transport
break-ratio prediction (rect 0.868, gyroid 0.908 at 0.2mm/200MeV).

- **Rectilinear (clean):** at res8 the geometry is ~converged (gamma2_straight 2.30 ~=
  analytic 2.18, f 0.39) at a clean-step voxel (0.025 mm); Geant4 gives **gamma2 = 1.97**
  (substantial -> homogeneous-shape FAILS at the foam-scale cell, 200 MeV; 0 stuck) with
  **ratio 0.854 matching the transport prediction 0.868 to ~1.5%.** 0.2 mm is just below
  cell_break(200)=0.267 mm, so ratio ~0.85 (mild departure onset), not 1, is expected.
- **Gyroid (geometry-resolution-limited):** at clean-step res (4,8) the smooth TPMS gives
  f = 0.50 (vs 0.40) and gamma2_straight 28% under converged; convergence needs res16+
  (voxel/step 0.0125 mm < 0.02 mm floor). **No voxel resolution simultaneously resolves
  the geometry AND keeps the step above the MSC floor**, so a clean gyroid foam-scale
  point is unattainable; its ratio (0.868) is consistent with transport (0.908) but not
  on a converged geometry. (Convergence diagnostic, gamma2_straight vs res: rect
  1.50/2.30/2.00/2.20 -> analytic 2.18; gyroid 0.75/0.84/1.24/1.17 -> 1.16, at res
  4/8/16/32.) See `transport_vs_g4.md`.

## Verdict (foam-failure restated on cell_homog; momentum-gating removed)
- The homogeneous-density (Highland) approximation's **fourth-moment / shape** fails for
  3D-printed lattices AND detector carbon foams (~100-500 um pores) across **all practical
  proton momenta (>~100 MeV)**: gamma2 ~ 1-2.2, Geant4-measured at 2.5 mm (both topologies,
  two energies) and at 0.2 mm (rectilinear).
- The previously-reported momentum gating (500 um > ~111 MeV, 250 um > ~209, 100 um >
  ~567 MeV) and the "therapy-energy fine foams stay averaged" caveat were defined on
  cell_break (the line-integral boundary) and **understated the failure** — they are
  REMOVED. The homogeneous limit returns only below ~um cells or below ~tens of MeV, both
  outside the practical detector regime.
- **Scope (kurtosis, not width):** Result 1 holds — kappa2_struct/kappa2_solid@fL = 1.000
  to <8%; it is the Gaussian tail/kurtosis, not the RMS, that fails.
- **Honest limitation:** the precise cell_homog (~um) rests on the transport tool, now
  Geant4-validated at 2.5 mm (both topologies) and 0.2 mm (rectilinear); the gyroid
  foam-scale point is geometry-resolution-limited by the fundamental voxel
  geometry-vs-MSC-step tension. The qualitative broad failure is measured; the precise
  recovery scale is transport-derived.

Figure: `figs/fig_homog_boundary.pdf` (fig:boundary).
