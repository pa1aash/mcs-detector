# e5a_rawk4_check.md — S7 / E5(a) raw-κ4 MSC diagnostic (the voronoi-steepening branch)

Decides whether the alt-MSC (Urban/option3) **voronoi** collapse-steepening (all-inclusive
slope −1.18 vs periodic-only −1.00; voronoi closure k_eff → ~0.5) is a floor/deconvolution
fragility (→ option 1) or a real MSC sensitivity in the stochastic family (→ option 3).
500 MeV, existing Urban + WentzelVI campaign outputs, **no re-sim**.

**Premise of the test:** a floor/deconvolution artifact lives only in the *subtracted* floor,
so the RAW (pre-subtraction) structured κ4 should rescale uniformly across families; a real MSC
effect should move the voronoi raw κ4 differently from the periodic anchor.

## Urban/WentzelVI ratios per config (500 MeV)

| family | config | geom frac of raw κ4 (W) | raw κ4_struct U/W (CI) | floor U/W | struct−floor gap | **geom Δκ4 U/W (CI)** |
|--|--|--:|--:|--:|--:|--:|
| voronoi | voronoi_f40 | 6.8% | 0.352 [0.349, 0.355] | 0.346 | +0.006 | **0.44 [0.39, 0.49]** |
| voronoi | voronoi_f50 | 6.0% | 0.346 [0.343, 0.349] | 0.347 | −0.001 | **0.33 [0.28, 0.38]** |
| periodic | rectilinear_f40 | 26.5% | 0.399 [0.389, 0.409] | 0.350 | +0.049 | 0.53 [0.49, 0.58] |
| periodic | schwarzp_f40 | 27.1% | 0.422 [0.412, 0.431] | 0.350 | +0.072 | 0.61 [0.57, 0.66] |
| periodic | gyroid_f40 | 16.0% | 0.390 [0.383, 0.397] | 0.347 | +0.043 | 0.62 [0.56, 0.67] |
| periodic | gyroid_f20 | 14.0% | 0.376 [0.367, 0.386] | 0.343 | +0.033 | 0.58 [0.50, 0.67] |

- **Floor rescales UNIFORMLY** under Urban: U/W = 0.347 ± 0.003 for every geometry (the
  intrinsic κ_M is ~0.35× the WentzelVI floor, geometry-independent).
- **Geometry Δκ4 rescales uniformly for the periodic family**: U/W = 0.587 ± 0.038 (a flat,
  N_eff-independent factor — this is exactly why the **periodic-only exponent is −1.00 under
  both physics lists**). Periodic geometry is 14–27% of the raw κ4, so Δκ4 is well-measured.
- **Voronoi geometry is only 6–7% of the raw κ4** → the raw κ4 is floor-DOMINATED, so the raw
  κ4 ratio (0.349) ≈ the floor ratio (0.347) *by construction*, NOT because of geometry
  physics. **The clean raw-κ4 test is therefore CONFOUNDED for voronoi** — the raw quantity
  cannot isolate the geometry there.

## Branch — CONFOUNDED as a clean A/B proof; the prescription is unambiguous

The voronoi geometry Δκ4 is the residual of a **struct−floor ratio gap of ±0.006** (f40 +0.006,
f50 −0.001) — i.e. the structured κ4 and the floor rescale to within ~0.5% of *each other*, and
the small voronoi geometry signal is that difference. The gap sits at the floor-model-accuracy
level (the κ_M = b·t + c·t² fit, evaluated over the voronoi's skewed t_pla — the documented
small-signal-on-large-floor fragility). A sub-percent floor-model error fully accounts for the
voronoi geometry Δκ4 dropping to ×0.33–0.44 (vs periodic ×0.59). The diagnostic therefore
**cannot cleanly separate a small real MSC differential from a floor-model error** for voronoi —
both produce the same sub-% residual. (Note the simple single-control subtraction shows the same
voronoi extra-drop, ×0.37–0.45 vs periodic ×0.60–0.73, so it is not specific to the all-order
floor — it is intrinsic to deconvolving a 6%-of-κ4 signal under an MSC change.)

**But the prescription is the same under either interpretation:** the voronoi quantitative
Δκ4 / k_eff is **not MSC-robust** — its closure swings k_eff ≈ 1.0 (WentzelVI) → ≈ 0.5 (Urban),
resting on a ±0.006 gap the floor model cannot pin to that precision. It cannot carry
quantitative load. This is the **2nd documented voronoi floor fragility** (the 1st: the spurious
~2× excess under the simple subtraction, S6 Stage-1).

**Resolution (option-1 prescription — surfaced for the user's ruling):**
1. **Headline exponent = PERIODIC-only −1.00** (rect/schwarzp/gyroid), MSC-robust under both
   physics lists — geometry is 14–27% of the raw κ4, rescales by a flat ×0.587, no N_eff trend.
2. **Voronoi DEMOTED to qualitative-only** — the stochastic exemplar for the cell-independence
   N_eff dichotomy (N_eff ∝ 1/cell), with NO load-bearing quantitative k_eff.
3. The all-inclusive **−1.18** is reported as inflated by the floor-dominated, hypersensitive
   voronoi points; the **prefactor MSC band 0.80–0.98** (periodic) confirms the slope/prefactor
   split (exponent robust; the a_eff-calibrated prefactor shifts uniformly) rather than
   threatening the law.
4. Exponent MSC band quoted as **−1.00 (periodic, load-bearing) → −1.18 (all-inclusive)** WITH
   the voronoi-floor-fragility attribution.

**Honest caveat (not "PROVEN"):** the raw-κ4 test was confounded by voronoi's floor-dominance,
so this is NOT a clean proof of "floor artifact over real MSC." It IS a proof that the voronoi
quantitative signal is not MSC-stable and must not carry quantitative load — which is what the
prescription requires. Surfaced for the user's ruling before the budget / S8.
