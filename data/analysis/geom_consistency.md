# Geometry consistency + detector-integration validations (S5.5)

Production `mcs_sim` now macro-selects the target geometry: `slab` (S3 G4Box,
default, unchanged), `stl` (CADMesh ASCII unit cell tiled transversely + along z),
`voxel` (binary field as a G4 parameterisation), and a `reflective` single-cell
option. These are the four gating validations that the bare integration plan
lacked. Physics list / cuts / maxStep unchanged from S3.

## 4. GEOMETRY CONSISTENCY (the STOP gate) — **PASS**

The Geant4 detector must build the SAME geometry the ray-tracer characterised, or
the line-integral prediction and the Geant4 measurement are not comparable. Probe:
straight zero-scattering **geantinos** fired across one transverse cell of the
tiled lattice (printable cell 2.0 mm, 5×5 transverse tiling, depth L=10 mm); the
per-primary PLA path length (`tpla`) is t = ∫χ dz, so ⟨t⟩/L = the Geant4-seen
solid fraction f and Var(t) = the chord variance. Targets = the S4 ray-tracer
(`data/geom_stats/<topology>_f40.npz`). Tolerances: f <1%, Var(t) <3% (the
ray-tracer estimator itself carries ~2%). 200k geantinos.

| topology | f (Geant4) | f (ray-tracer) | Δf | Var(t) G4 | Var(t) ray | ΔVar | verdict |
|--|--:|--:|--:|--:|--:|--:|:--:|
| rectilinear | 0.3987 | 0.4008 | −0.52% | 11.544 | 11.528 | +0.14% | **OK** |
| gyroid      | 0.3988 | 0.4007 | −0.47% | 6.168  | 6.215  | −0.74% | **OK** |

**Both f AND Var(t) match within tolerance → PASS.** The CADMesh load, the
periodic transverse tiling, and the depth tiling reproduce the ray-traced p(t)
the N_eff theory and the S5 line-integral prediction are built on.

*Provenance note:* this required regenerating the E0 unit-cell STLs at higher
marching-cubes resolution (spc=80 vs the S4 spc=40). At spc=40 the rectilinear
**cylindrical struts were under-resolved**, giving a meshed f 0.95% below analytic
→ a 1.78% Geant4-vs-ray f gap that FAILED the <1% gate. At spc=80 the meshed f is
within +0.22% of analytic and the gate passes (gyroid was already exact at spc=40).
See `geom/regen_e0_cells.py`, `geom/stl/e0_pair/manifest.csv`.

**The voxel path (the campaign's actual proton path, see §6) is ALSO geometry-
consistent** — geantino probe through voxelised rect/gyroid blocks (L=10 mm,
voxel 0.0625 mm, `geom/make_lattice_voxel.py`):

| topology | f (voxel G4) | Δf vs ray | Var(t) voxel | ΔVar |
|--|--:|--:|--:|--:|
| rectilinear | 0.3980 | −0.71% | 11.565 | +0.32% |
| gyroid      | 0.3968 | −0.96% | 6.223  | +0.12% |

Both <1%/<3% → the geometry is consistent whether built as tessellated tiling
(geometry probing) or voxels (proton physics).

## 5. VOXEL-VS-TESSELLATED (Voronoi) — **PASS (geometry); κ within bootstrap**

Voronoi is the theory-critical stochastic geometry (anchors the L/2ℓ_int scaling);
its ~461k-facet tessellation is why S4 chose to voxelise it. One Voronoi block
(5 mm cube, cell 2.5 mm, f=0.40) was emitted as BOTH a voxel field (128³,
voxel 0.039 mm) and a tessellated STL from the **same realisation**
(`geom/make_voronoi_voxel.py`).

Geometry (geantino, same beam over the central 2.5 mm column):

| representation | f (⟨t⟩/S) | Var(t) |
|--|--:|--:|
| tessellated (461k facets) | 0.3755 | 0.953 |
| voxel (128³)              | 0.3765 | 0.968 |
| agreement                 | +0.27% | +1.6% |

Kink angles (200 MeV protons, common fixed window 5·σ_core; 80k each):

| representation | σ_core [mrad] | κ₂ [rad²] | κ₄ [rad⁴] | κ₄ 95% CI |
|--|--:|--:|--:|--:|
| tessellated | 2.143 | 5.828e-06 | 5.102e-11 | [4.90e-11, 5.30e-11] |
| voxel       | 2.126 | 5.733e-06 | 4.963e-11 | [4.76e-11, 5.18e-11] |
| agreement   | −0.8% | −1.64% | −2.72% | **CIs overlap → PASS** |

Voxel and tessellated agree on f, Var(t), and the kink κ₂/κ₄ → **the 0.039 mm
voxelisation is artifact-free**; the Voronoi may be voxelised for S6 without
corrupting the dichotomy.

## 6. REFLECTIVE-VS-TRUE-MULTI-CELL (break-specific) — **DECISION: reflective INVALID; build campaign as VOXELS**

The break is about lateral wandering across cells, so a transverse boundary that
folds cross-cell sampling could bias exactly the S5 measurement. Tested at a
near-break cell (0.05 mm, rectilinear, depth 10 mm = 200 cells deep):

| build | f (geantino) | Var(t) | stuck-track warnings | f (200 MeV proton) |
|--|--:|--:|--:|--:|
| reflective single cell (specular momentum-fold) | 0.3899 | 11.55 | 1389/50k | **0.263** (35% low) |
| true multi-cell (7×7 transverse tiling) | 0.3993 | 11.58 | 1388/50k | — |

**The reflective momentum-fold is INVALID.** For straight geantinos it is already
biased (f −2.7%) and throws stuck-track warnings; for the actual **scattering
protons it collapses to f=0.263 (35% low)** — the fold corrupts the very
lateral-transport physics the break depends on. Per the task's pre-registered
contingency, **S5 must use true multi-cell tiling blocks, not reflective.** The
`/mcs/det/reflective` option is retained but flagged INVALID-for-physics in the
code.

**Secondary finding — the PROTON BUILD PATH (carried to S5/S6, important):**
abutting void-padded unit-cell tessellations are **proton-navigation-fragile**.
Measured stuck-track rates (200 MeV protons):

| build | stuck tracks | note |
|--|--:|--|
| tiled unit cells — **rectilinear** c2 (5×5×5) | **1720 / 10000 (17%)** | axis-aligned struts coplanar with cell faces → coincident surfaces |
| tiled unit cells — gyroid c2 (5×5×5) | 0 / 10000 | smooth TPMS crosses faces at angles |
| **voxel** (rect & gyroid & Voronoi) | **0** | regular grid, no inter-cell surfaces |
| single watertight block STL (Voronoi, 461k facets) | 0 | one closed surface |

Straight geantinos rarely cross transverse boundaries, so the §4 geometry-
consistency probe is unaffected; but scattering **protons** repeatedly hit the
coincident inter-cell caps of the tiled rectilinear build. **Decision: the proton
campaign (S5 break sweep + S6) must build lattices as VOXELS** (0 stuck, geometry-
consistent per §4, artifact-free per §5) **or single watertight multi-cell STLs**,
NOT abutted unit-cell tiling. Voxel generators: `geom/make_lattice_voxel.py`
(rect/gyroid), `geom/make_voronoi_voxel.py` (Voronoi). The tiled-STL path is kept
only for geantino geometry probing.

## 7. HIGHLAND REGRESSION (slab path unperturbed) — **PASS**

The solid PLA slab campaign re-run through the modified binary reproduces S3
exactly: σ_core/θ₀−1 at 16 mm = −3.07 / −3.65 / −3.31% (200/500/1000 MeV),
κ₄ linearity R²≈0.987. The z-plane scoring equals the S3 PV-transition scoring on
the slab. Details appended to `VALIDATION.md §5`.

## Summary

| validation | result |
|--|--|
| 4. Geometry consistency (f AND Var match ray-tracer) | **PASS** — tessellated Δf<0.6%/ΔVar<0.8%; voxel Δf<1%/ΔVar<0.4% |
| 5. Voxel vs tessellated (Voronoi) | **PASS** — geometry f+0.27%/Var+1.6%; proton κ₂−1.6%, κ₄ CIs overlap |
| 6. Reflective vs multi-cell | **reflective INVALID**; + tiled-rect proton-fragile → **build campaign as VOXELS** |
| 7. Highland regression (slab) | **PASS** — identical to S3 (−3.07/−3.65/−3.31%) |

**GATE 5.5: PASS.** `mcs_sim` builds; macro-selects slab/stl/voxel/reflective;
geometry consistency holds (Geant4-seen f AND Var(t) match the ray-tracer for both
the tessellated and voxel builds); voxel reproduces the tessellated Voronoi; the
reflective-vs-multi-cell decision is made (reflective invalid → multi-cell, and
specifically **voxels** for the proton campaign because abutted rectilinear unit-
cell tiling is navigation-fragile); Highland regression still passes on the slab;
sub-printable test cells {0.2,0.03,0.01} regenerated; smoke tests run (voxel/
single-block clean, tiled-rect flagged). **No κ₄ campaign was run** — that is the
re-run S5.

Empirical S5/S6 (Δκ₄ sweep) is now UNBLOCKED, building lattices as voxels.
Compute still to be provisioned (cluster; the navigation-heavy deep small-cell
configs at 1e6–1e7 stats are not laptop runs).
