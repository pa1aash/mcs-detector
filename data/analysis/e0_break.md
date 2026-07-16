# E0 — the homogeneous-approximation break (transport-aware map)

Method: transport-aware ray-tracer on ANALYTIC geometry (no facets/voxels/steps -> no resolution/MSC artifact), validated against full Geant4 at printable cells (`transport_vs_g4.md`). Break = cell where Var(t_actual)/Var(t_straight)=Δκ4_TA/Δκ4_LI departs from 1 by >10%. L=10.0 mm, f=0.4, 20000 protons/point, 20 steps/cell, 2 seeds.

## Break table — cell_break(momentum) [µm]

| p [MeV] | rectilinear | gyroid | y_rms(p) solid [µm] | vs FDM 0.5 mm |
|--:|--:|--:|--:|--:|
| 100 | 559 | 376 | 76 | 1.12× |
| 200 | 260 | 183 | 36 | 0.52× |
| 500 | 111 | 76 | 15 | 0.22× |
| 1000 | 62 | 44 | 9 | 0.12× |

## Ratio map (Var(t_act)/Var(t_str))

**rectilinear**

| p\cell | 3 | 2 | 1 | 0.7 | 0.5 | 0.3 | 0.2 | 0.15 | 0.1 | 0.07 | 0.05 | 0.03 |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 100 | 0.99 | 0.97 | 0.94 | 0.92 | 0.89 | 0.82 | 0.74 | 0.68 | 0.60 | 0.53 | 0.47 | 0.41 |
| 200 | 0.99 | 0.99 | 0.98 | 0.97 | 0.95 | 0.92 | 0.87 | 0.83 | 0.75 | 0.68 | 0.61 | 0.51 |
| 500 | 1.00 | 0.99 | 0.99 | 0.99 | 0.98 | 0.97 | 0.95 | 0.93 | 0.89 | 0.84 | 0.78 | 0.68 |
| 1000 | 1.00 | 1.00 | 0.99 | 0.99 | 0.99 | 0.98 | 0.97 | 0.96 | 0.94 | 0.91 | 0.88 | 0.80 |

**gyroid**

| p\cell | 3 | 2 | 1 | 0.7 | 0.5 | 0.3 | 0.2 | 0.15 | 0.1 | 0.07 | 0.05 | 0.03 |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 100 | 1.00 | 0.99 | 0.97 | 0.96 | 0.93 | 0.87 | 0.81 | 0.76 | 0.68 | 0.61 | 0.55 | 0.49 |
| 200 | 1.00 | 1.00 | 0.99 | 0.99 | 0.98 | 0.95 | 0.91 | 0.88 | 0.82 | 0.75 | 0.69 | 0.59 |
| 500 | 1.00 | 1.00 | 1.00 | 1.00 | 0.99 | 0.99 | 0.97 | 0.96 | 0.93 | 0.89 | 0.84 | 0.76 |
| 1000 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.99 | 0.98 | 0.97 | 0.95 | 0.92 | 0.85 |

## Topology-contrast test

- 100 MeV: rect break 559 µm vs gyroid 376 µm (differ 33%)
- 200 MeV: rect break 260 µm vs gyroid 183 µm (differ 30%)
- 500 MeV: rect break 111 µm vs gyroid 76 µm (differ 32%)
- 1000 MeV: rect break 62 µm vs gyroid 44 µm (differ 29%)

**Topology contrast: SAME cell within ~33% -> the break is a near-universal length-scale criterion set by the lateral-scattering scale, weakly topology-dependent.**

## Scope (Phase D honesty)

The break map uses the transport-aware tool on analytic geometry. It is validated against full Geant4 at printable cells (≥0.5 mm, see `transport_vs_g4.md`): Var(t_actual) and Δκ4 agree within bootstrap error there. Below Geant4's validated step regime (sub-printable, fine cells) NO Geant4 measurement is claimed — the break there rests on the transport-aware map, justified by the printable-cell agreement where the regimes overlap. The FDM print limit is 0.5 mm; fine SLA ~0.05 mm.

