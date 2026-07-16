# `analysis/` — shared analysis library (built S3+)

`lib/` holds the shared cumulant / MLE utilities every experiment imports. No
analysis logic is duplicated in figure scripts or run drivers — it lives here.

> **Status: empty (S0).** `lib/` populated from S3 onward.

## Intended `lib/` contents

- **Cumulant estimators** — `κ₄` with convolution / **deconvolution** algebra:
  subtract telescope ⊕ air at the `κ₄` level (fourth cumulants add under
  convolution; fourth excess-kurtoses do not). This is the headline-statistic
  machinery (`ROADMAP_SOURCE.md` §3).
- **Bootstrap CIs** for the cumulant and shape statistics.
- **Robust quantile-based shape statistic** (the triad's second leg).
- **Two-component Gaussian-mixture MLE** with fixed resolution kernel — mixture
  weight estimates `f` → closure test against mass-measured `f`.
- **Anderson–Darling GoF** — mixture model vs single Gaussian.
- **Highland / Lynch–Dahl width references** for SC1.

Headline statistic is the deconvolved fourth **cumulant + mixture MLE**, *not*
raw sample kurtosis (`ROADMAP_SOURCE.md` §3, §13).
