# Floor-model third-form robustness (FLOOR-3FORM; Fable 4.2)

Regenerated on the pod (2026-07-16, Geant4 11.3.2, g4highland) because the raw solid-control
and periodic campaign `.root` files are not committed. Solid controls t={2,3,4,5,8,16} mm at
1e7 protons and the periodic campaign (rectilinear/Schwarz-P 3e6, gyroid 6e6; f=0.20-0.50) at
200 and 500 MeV, locked physics/windowing. Analysis: `analysis/e2_floor_thirdform.py`.

Three floor forms fit to the SAME six solid controls, then the periodic mean scale-mixture
closure k_eff = Delta_kappa4 / (3 a_eff^2 Var(t_pla)) re-derived under each:

| form | 200 MeV k_eff | 500 MeV k_eff | Delta vs quadratic |
|---|--:|--:|--:|
| quadratic  b t + c t^2  (LOCKED, current) | 0.9902 | 1.0088 | 0.00% / 0.00% |
| cubic      b t + c t^2 + d t^3            | 0.9943 | 1.0091 | +0.41% / +0.03% |
| PCHIP monotone spline (non-parametric)    | 1.0006 | 1.0137 | **+1.04%** / +0.48% |

## Verdict (honest branch)

- The periodic closure is robust to the floor functional form at the **~1% level**: <0.5%
  under a cubic polynomial, at most +1.04% under a non-parametric monotone spline (at 200 MeV).
- The scaffold's rigid >1% gate tripped on the PCHIP 200 MeV point (+1.04%), so this is
  FLAGGED (see FLAGGED_FOR_REVIEW.md) rather than reported as "<1% robust".
- **Interpretation:** the ~1% floor-form sensitivity is comparable to the ±0.01
  per-configuration closure uncertainty already quoted, so the floor model does NOT carry the
  periodic signal (unlike the stochastic Voronoi, whose k_eff swings 2.05->0.96 with the floor).
  The periodic headline (k_eff = 1.00±0.01; exponent -0.93) is unaffected at the level that
  matters. §8.2 now states this honestly (cubic <0.5%, spline ~1%), not a false sub-1% claim.
