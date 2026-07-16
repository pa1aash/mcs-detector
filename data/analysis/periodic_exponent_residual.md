# F2 investigation — the periodic-only pooled exponent's residual shallowness vs −1

**Audit ID:** F2 (fresh-referee round). **Date:** 2026-07-16.
**Question:** §8.1 said the periodic-family exponent is "consistent with −1," but the printed
95% CI $-0.932\,[-0.992,-0.872]$ excludes −1. Is there a real mechanism, or is the CI wrong?

Reproduced fresh from `data/analysis/e2_results_combined.json` via the `collapse_fit` routine
(`analysis/e2_analysis.py`); script `scratchpad/f2_investigate.py` + `f2_refine.py`.

## A2 — the CI is correct, not underestimated, not mis-printed
- Reproduced periodic pooled (200+500 MeV, rect+schwarzp+gyroid, n=24) OLS free slope
  **−0.9322 [−0.993, −0.871]**, matching the printed **−0.932 [−0.992, −0.872]** to rounding.
  No printing/computational error.
- Independent **configuration-resampling bootstrap** (resample the 24 configs with replacement,
  refit): **−0.932 [−0.995, −0.869]** — essentially identical to the Birge MC-error-propagation
  CI. The Birge CI is **not** underestimated.
- −1 sits at the ~1.5–1.7 percentile of the slope distribution (~2.1σ one-sided): genuinely
  just outside 95%. The "excludes −1/2 and −2" claims are safe (far outside).

## A1 — the shallowness is a structured straight-chord-N_eff residual, NOT noise
Per-topology mean C·N_eff (± SEM, n=8 each), monotonically ordered by N_eff:

| topology    | N_eff | ⟨C·N_eff⟩ | note |
|-------------|------:|----------:|------|
| rectilinear | ~2.09 | 0.955 ± 0.013 | lowest N_eff, sits ~5% below unity |
| Schwarz-P   | ~2.25 | 0.998 ± 0.011 | ~unity |
| gyroid      | ~4.38 | 1.023 ± 0.017 | highest N_eff, ~2% above unity |

- rect-to-gyroid separation 0.068, combined SEM 0.021 → **~3.2σ** (not sampling noise).
- The offset aligns with N_eff, so a global log-log slope reads it as a **+0.068 trend of
  C·N_eff with N_eff** (== free_slope+1). Over the periodic N_eff span (×2.7) this is a +6.9%
  rise. The same +0.065 residual recurs **within** the gyroid family over its own N_eff range
  (gyroid-only slope −0.935), so it is a general property of the scalar, not one outlier.
- **Energy dependence:** gyroid ⟨C·N_eff⟩ rises 0.998 (200 MeV) → 1.049 (500 MeV); 200-only
  periodic slope −0.974, 500-only −0.891. The shallowness grows mildly with energy.
- Single-point leave-one-out: max slope influence ±0.03 (the two highest-N_eff gyroid f=0.20
  points, opposite signs); the effect is not a single-point artifact.

## Interpretation (honest branch = MECHANISM FOUND)
The straight-chord $N_\mathrm{eff}=f(1-f)L^2/\mathrm{Var}(t)$ is a geometric proxy that
captures the collapse to ≈5–7% but does **not** perfectly equalise the periodic families:
rectilinear's effective cell count is slightly under-counted relative to gyroid's, leaving a
~3σ topology-prefactor spread that a free log-log slope reads as −0.93 rather than −1. The
parameter-free −1 is recovered **not** in the free slope but in the measured-path closure
($k_\mathrm{eff}=\Delta\kappa_4/[3a_\mathrm{eff}^2\mathrm{Var}(t_\mathrm{pla})]=1.00\pm0.01$)
and with the slope fixed at −1 (periodic prefactor 0.99; geometric-mean C·N_eff 0.99).

**Consequence for the text:** the "consistent with −1" language for the periodic *free* slope
is replaced by (i) the quantified topology-prefactor mechanism (§8.1) and (ii) the honest
statement that the free periodic slope excludes −1 at ~2σ while −1 survives as the closure
statement. Applied in §8.1, §5.3, Introduction, Conclusions, and abstract.
