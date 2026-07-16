# `geom/` — topology generators + ray-tracer (built S2)

Generates the five target topologies as STL/GDML for the Geant4 app, and
ray-traces each to produce the geometry statistics behind the `N_eff` collapse
(E1).

> **Status: empty (S0).** Implemented in S2.

## Two jobs

**1. Topology generators** → STL/GDML + parametric placement. The five
topologies (final set fixed in S2):

1. solid control
2. extruded rectilinear lattice
3. cubic / strut lattice
4. gyroid (TPMS)
5. a second TPMS or stochastic foam

Plus the required **equal-areal-density solid slab** baseline (~4 mm solid PLA,
matched to the 40% lattice). Infill variants (~4 fractions including 40%) and a
cell-size sweep hook feed E0.

**2. Ray-tracer** → cast beam-direction rays through each geometry and emit:
- `p(t)` — distribution of material line-integrals along the beam
- `S2(z)` — two-point (autocorrelation) function of the material indicator
- `ℓ_int` — directional integral correlation length
- `N_eff = L / (2 ℓ_int)`

Outputs land in `data/geom_stats/`. CPU-days, not cluster-scale.
