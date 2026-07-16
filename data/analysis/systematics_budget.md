# systematics_budget.md — S7 robustness + systematics budget

The three S7 systematics (E5 MSC model, E3 affordability, E4 print-realism), each with its
number and scope. The organising principle is the **slope/prefactor split**: the load-bearing
result is the *parameter-free exponent* of the structure-induced kurtosis law
Δκ4 = 3a²Var(t) ∝ N_eff⁻¹ (SC2); the *prefactor/absolute amplitude* is a calibrated,
weaker corroboration that legitimately carries model- and print-systematic bands.

---

## E5 — MSC-model systematic (the gate; run first)

Alternative MSC = `G4EmStandardPhysics_option3` (UrbanMsc-based) vs the locked `option4`
(WentzelVI multiple scattering + single Coulomb). A deliberately LARGE, conservative stress
test: under Urban the proton core narrows (Highland residual −9.2% vs −3.1%), the effective
scattering power **a_eff drops 17%** (ratio 0.834), and the intrinsic κ_M floor is ~3× smaller.
Geometry (ray-traced N_eff) is byte-identical to the locked campaign, isolating the MSC effect.

### E5(a) — THE COLLAPSE (clean regime)  [GATE: PASS]

- **Exponent (load-bearing, parameter-free): the periodic-lattice OLS slope is −1.00 under
  BOTH physics lists** (alt −1.002, locked −0.987; rect/schwarzp/gyroid, N_eff 2–5.5). The
  SC2 N_eff⁻¹ law is **MSC-robust** — the geometry sets N_eff, and the MSC model rescales
  a_eff uniformly, which cancels in the collapse coordinate C = Δκ4/(3 a_eff² f(1−f) L²).
- **Exponent MSC band: −1.00 (periodic, load-bearing) → −1.18 (all-inclusive).** The
  all-inclusive slope steepens to −1.182 [−1.416, −0.946] **solely** because the two
  stochastic-voronoi points' closure k_eff swings 1.0 (WentzelVI) → ~0.5 (Urban). The
  raw-κ4 diagnostic (`e5a_rawk4_check.md`) shows the voronoi geometry signal is only ~6% of
  the raw κ4 (floor-dominated) and rests on a struct−floor ratio gap of ±0.006, at the
  floor-model-accuracy level — i.e. the voronoi quantitative k_eff is **not MSC-stable**
  (its 2nd documented floor fragility). **Ruling (user): demote voronoi to qualitative-only**
  (the stochastic exemplar for the N_eff ∝ 1/cell dichotomy; no load-bearing quantitative
  k_eff). Honest scope: the diagnostic proves voronoi is unreliable, not "floor-artifact over
  small-real-MSC" (the raw test is confounded by voronoi's floor-dominance).
- **Prefactor band: k = 0.98 (locked) → 0.80 (alt), uniform across N_eff** — this confirms the
  slope/prefactor split (the a_eff-calibrated prefactor shifts uniformly without touching the
  exponent), it does not threaten the law.
- **Absolute Δκ4 band: −40 to −70%** under Urban (the expected a_eff⁻rescaling: a² → −30%,
  plus the smaller floor). This is the model-systematic on the absolute amplitude, to be
  REPORTED — distinct from the exponent, which is robust.

### E5(b) — THE BOUNDARY

Re-derived under Urban: a_eff −16% (ratio 0.84 at both 200 & 500 MeV), Highland residual −9%
(vs −3% locked), intrinsic κ_M floor ~0.35×. Two SEPARATE results (do not conflate):

1. **The QUALITATIVE foam-shape failure is MSC-robust (Geant4-anchored).** At the 0.2 mm
   rect foam cell (200 MeV), the Geant4-measured excess kurtosis is γ2 = 1.97 (WentzelVI)
   and **γ2 = 0.90 (Urban)** — both far from the homogeneous γ2 = 0 (geometry is 17–24% of the
   raw κ4, a reliable fraction — not the voronoi small-difference fragility). **The broad
   homogeneous-shape failure does NOT vanish under the MSC change**, and it does not depend on
   the transport tool. This is the load-bearing qualitative claim, and it survives.
2. **The PRECISE magnitude / cell_homog carries an MSC systematic the fine-cell regime forbids
   bounding cleanly (stated limitation).** The foam γ2 magnitude has a **~2× MSC band (0.9–2.0)**.
   The transport tool predicts γ2@0.2mm is a-independent (2.13 under both physics — it carries
   `a` only through the wandering), but the Geant4 measurement shows the magnitude IS
   MSC-sensitive: the scale-mixture closure the transport tool rests on degrades toward fine
   cells under Urban (k_eff ≈ 0.8 at the 2.5 mm printable cell → ≈ 0.4 at the 0.2 mm foam cell).
   So the precise cell_homog (~µm) **cannot be cleanly MSC-bounded** — it is reported as ~µm
   (sub-foam under both physics, so the conclusion is qualitatively unchanged) WITH this
   explicit transport-tool MSC limitation, not as a clean ±X number. cell_homog itself shifts
   only ~30% in the tool (e.g. rect 0.06→0.08 µm, gyroid 0.13→0.16 µm at 200 MeV), but that
   tool number inherits the scale-mixture MSC band above.

---

## E3 — Affordability of the shape channel (calibration-anchored)

From the Phase-0c calibration + campaign (measured deconvolved-Δκ4 CI scaled by the
1/√N moment law; the 470M-track figure is BANNED). Flux anchor (verified): the LLUMC/UCSC
phase-II proton-CT scanner, **>1 MHz** individually-measured protons at ~200 MeV through a
Si-strip tracker, full scan ≤6 min [Johnson et al., *Phys. Procedia* **90** (2017) 209].

- **Protons / beam-time to resolve the geometry Δκ4** (representative configs, 1 MHz):
  - rectilinear f40 (N_eff 2): 30% → 3.3e5 (0.3 s), 10% → 2.9e6 (2.9 s).
  - gyroid f40 (N_eff 4): 30% → 1.0e6 (1.0 s), 10% → 9.3e6 (9.3 s).
  - voronoi f30 (N_eff 13.5): 30% → 1.9e7 (19 s), 10% → 1.7e8 (2.8 min).
  - **Printable/foam-relevant configs are resolvable to 10% in <10 s of 1-MHz beam; the
    hardest high-N_eff case in ~minutes. Feasible at a proton-CT-class facility.**
- **Cost RATIO vs the standard WIDTH (Highland) channel** (equal fractional precision, both
  ∝1/√N): **~10³–10⁴× for the printable configs (rect/gyroid), rising to ~10⁶× for the
  hardest high-N_eff voronoi.** Decomposed = moment-order penalty (4th vs 2nd moment, ~30–70×)
  × small-signal/deconvolution penalty (Δκ4 is a small excess on a large floor, ~40–2700×).
- **Honest:** the shape channel is MORE expensive than width — the multiplier is the real
  price — but affordable in ABSOLUTE terms (seconds-to-minutes at 1 MHz) because width is so
  cheap. See `e3_affordability.md`.

---

## E4 — Print-realism (FDM artifacts on the printable collapse)

Realistic FDM artifacts perturb the campaign voxel geometry on a subset (rect/gyroid/voronoi
f40), quantified through the Geant4-validated ray-tracer (at the 2.5 mm cell the straight-chord
Var(t) matches Geant4 to <1%, S5); the dominant artifact is Geant4-spot-checked.

- **Dominant systematic = dimensional tolerance (strut/wall over/under-extrusion): Δκ4 shifts
  up to ~18%** (it moves f by +12 to −32% and Var(t) with it). Surface roughness / layer lines
  ~−8 to −9%; microporosity ~−4 to −5%; local infill variation ~±2–4%.
- **Two honest statements:** (1) the ABSOLUTE Δκ4 is robust to ~5% under the uncontrollable
  artifacts but ~18% under dimensional tolerance, which must be controlled (or the wall
  thickness independently measured). (2) the COLLAPSE point is more robust: N_eff = f(1−f)L²/
  Var(t) and Δκ4 co-move under a tolerance offset, so the **collapse coordinate N_eff shifts
  ≤10%** across all artifacts — the points move along the N_eff⁻¹ line, and the parameter-free
  exponent is not threatened. See `e4_print_realism.md` (+ `e4_validate.json` Geant4 check).

---

## The qualitative-vs-precise foam-failure distinction (carried explicitly)

The two must not be conflated:

- **MSC-ROBUST (Geant4-anchored, load-bearing):** the *qualitative* broad homogeneous-shape
  failure — γ2 = O(1) ≫ 0 at the ~0.2 mm foam cell — is directly Geant4-confirmed under BOTH
  physics lists (γ2 = 1.97 WentzelVI, 0.90 Urban), independent of the transport tool. Foams
  (~100–500 µm) fail the homogeneous *kurtosis* across all practical momenta; this conclusion
  does not move under the MSC change.
- **MSC-SYSTEMATIC-LIMITED (transport-derived, reported with a band/limitation):** the *precise*
  magnitude (foam γ2 ~ 0.9–2.0, a ~2× band) and the precise cell_homog (~µm) rest on the
  transport tool, whose scale-mixture basis itself carries an MSC band that grows toward fine
  cells (k_eff 0.8→0.4). The precise µm boundary cannot be cleanly MSC-bounded — stated as a
  limitation, not a number.

## Overall verdict

| systematic | number | gate |
|--|--|--|
| **E5 MSC — exponent** | periodic **−1.00** at 200 & 500 MeV (MSC-robust); all-inclusive 500 −1.18 (voronoi-floor, demoted) | **PASS** (collapse not overturned) |
| **E5 MSC — prefactor / abs Δκ4** | k 0.98→0.80; absolute Δκ4 −40 to −70% (a_eff −16%) | reported band |
| **E5 MSC — foam failure** | qualitative γ2≫0 robust (1.97→0.90, Geant4 both lists); magnitude band ~0.9–2.0; precise cell_homog ~µm MSC-limited | **PASS** (qualitative failure not overturned) |
| **E3 affordability** | printable configs to 10% in <10 s @1 MHz (pCT); shape/width cost ~10³–10⁴× (rect/gyroid) | traced to calibration |
| **E4 print-realism** | dominant = dimensional tolerance, Δκ4 ±~10–18% (Geant4-confirmed); N_eff collapse-point ≤10% | dominant quantified |

**Bottom line.** The load-bearing results survive the systematics: the SC2 N_eff⁻¹ **exponent
is MSC-robust (−1.00, periodic, two energies and two physics lists)** and the **qualitative
homogeneous foam-failure is Geant4-anchored under both physics lists**. The systematics live in
the *amplitudes/prefactors*, exactly where the slope/prefactor split predicts they should: the
absolute Δκ4 carries a −40…−70% MSC band and a ±~10–18% print-tolerance band; the foam γ2
magnitude carries a ~2× MSC band; the precise µm cell_homog is transport-derived and not cleanly
MSC-boundable. The stochastic-voronoi quantitative k_eff is demoted to qualitative (2nd
documented floor fragility). The shape channel is affordable (seconds-to-minutes of proton-CT-
class beam) but ~10³–10⁴× costlier than the width channel. No systematic overturns the headline.
