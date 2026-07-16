#!/usr/bin/env python3
"""e5a_rawk4_check.py -- S7 / E5(a) DECISIVE diagnostic for the voronoi steepening.

The alt-MSC (Urban/option3) all-inclusive collapse slope is -1.18 (outside the locked
CI), steepened SOLELY by the stochastic-voronoi family whose closure k_eff drops to ~0.5.
Two hypotheses:
  (A) floor/deconvolution artifact -- the divergence lives in the SUBTRACTED floor (the
      Urban kappa_M is ~3x smaller; the all-order subtraction over-corrects the small-signal
      voronoi). Then the RAW (pre-subtraction) structured kappa4 rescales UNIFORMLY across
      families: ratio_raw(voronoi) ~= ratio_raw(periodic).  -> option 1.
  (B) real MSC sensitivity in the stochastic family -- the RAW structured kappa4 itself moves
      differently for voronoi than for the periodic anchor: ratio_raw(voronoi) departs from
      ratio_raw(periodic) beyond stats.  -> option 3.

Decisive test (existing data only; NO re-sim): for each config compute, under Urban and under
WentzelVI (locked), with bootstrap CIs:
  - k4_struct   = windowed structured kappa4, NO subtraction (the raw measured quantity)
  - floor       = <kappa_M(tpla)> (all-order), per physics list
  - dk4_simple  = k4_struct - k4_solid@f*L (simple single-control subtraction)
  - dk4_allord  = k4_struct - floor (the all-order deconvolution used in the collapse)
Then the Urban/WentzelVI RATIO of each, per config, and the family comparison.

Writes data/analysis/e5a_rawk4_check.{md,json}.  Run inside g4highland.
"""
from __future__ import annotations
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HERE, "lib"))
sys.path.insert(0, os.path.join(ROOT, "geom"))
import kink_stats as ks   # noqa
import e2_analysis as e2  # noqa

AOUT = os.path.join(ROOT, "data", "analysis")
LOCKED = os.path.join(ROOT, "data", "runs")
URBAN = os.path.join(ROOT, "data", "runs_e5", "opt3")
E = 500
W = e2.W[E]
# the two Urban voronoi infills that were run + the periodic anchors run under Urban
VORONOI = [("voronoi", 0.40), ("voronoi", 0.50)]
PERIODIC = [("rectilinear", 0.40), ("schwarzp", 0.40), ("gyroid", 0.40), ("gyroid", 0.20)]


def measure(runs_dir, topo, infill):
    """Raw struct kappa4 (+SE), matched solid kappa4, all-order floor, both subtractions."""
    e2.RUNS = runs_dir                      # redirect e2's build_floor / load_tpla / solids
    kM = e2.build_floor(E)
    tag = f"camp_{topo}_f{int(round(infill*100)):02d}_E{E}"
    ang = ks.load_run(os.path.join(runs_dir, tag + ".root")).angles
    tpla = e2.load_tpla(tag)
    k2, _ = ks.cumulants_in_window(ang, W)
    k4_struct, se = e2._k4_pt_se(ang, W, np.random.default_rng(11))
    tmm = round(infill * e2.L)
    sp = os.path.join(runs_dir, f"solid_E{E}_t{tmm}.root")
    solid = ks.load_run(sp).angles
    _, k4_solid = ks.cumulants_in_window(solid - np.mean(solid), W)
    floor = float(np.mean(kM(tpla)))
    return dict(k4_struct=float(k4_struct), se=float(se), k4_solid=float(k4_solid),
                floor=float(floor), dk4_simple=float(k4_struct - k4_solid),
                dk4_allord=float(k4_struct - floor), k2=float(k2),
                tpla_mean=float(np.mean(tpla)))


def ratio_ci(uU, seU, uW, seW):
    """ratio U/W and its ~95% CI via independent-error propagation."""
    r = uU / uW
    rel = np.sqrt((seU / uU) ** 2 + (seW / uW) ** 2) if uU and uW else float("nan")
    return r, r * (1 - 1.96 * rel), r * (1 + 1.96 * rel)


def main():
    rows = []
    for fam, lst in (("voronoi", VORONOI), ("periodic", PERIODIC)):
        for topo, infill in lst:
            U = measure(URBAN, topo, infill)
            Wm = measure(LOCKED, topo, infill)
            r_raw, lo_raw, hi_raw = ratio_ci(U["k4_struct"], U["se"], Wm["k4_struct"], Wm["se"])
            # floor SE is sub-dominant (solids 1e7); treat floor ratio as a point value
            r_floor = U["floor"] / Wm["floor"]
            r_simple = U["dk4_simple"] / Wm["dk4_simple"] if Wm["dk4_simple"] else float("nan")
            r_allord = U["dk4_allord"] / Wm["dk4_allord"] if Wm["dk4_allord"] else float("nan")
            rows.append(dict(family=fam, config=f"{topo}_f{int(infill*100)}",
                             ratio_raw_struct=r_raw, ratio_raw_ci=[lo_raw, hi_raw],
                             ratio_floor=r_floor, ratio_dk4_simple=r_simple,
                             ratio_dk4_allord=r_allord,
                             U=U, W=Wm))
            print(f"{fam:8s} {topo}_f{int(infill*100):<3} "
                  f"raw_k4 U/W={r_raw:.3f}[{lo_raw:.3f},{hi_raw:.3f}]  "
                  f"floor U/W={r_floor:.3f}  dk4(simple) U/W={r_simple:.2f}  "
                  f"dk4(allord) U/W={r_allord:.2f}")

    per_raw = np.array([r["ratio_raw_struct"] for r in rows if r["family"] == "periodic"])
    vor_raw = np.array([r["ratio_raw_struct"] for r in rows if r["family"] == "voronoi"])
    per_mean, per_sd = float(per_raw.mean()), float(per_raw.std(ddof=1))
    vor_mean = float(vor_raw.mean())
    # voronoi raw ratios within the periodic raw-ratio band (+ their own CIs)?
    vor_ci = [r["ratio_raw_ci"] for r in rows if r["family"] == "voronoi"]
    agree = all(ci[0] <= per_mean + 3 * per_sd and ci[1] >= per_mean - 3 * per_sd
                for ci in vor_ci)
    branch = ("A: floor/deconvolution artifact (raw kappa4 rescales uniformly; divergence is "
              "in the subtraction) -> OPTION 1") if agree else \
             ("B: real MSC sensitivity in the stochastic family (raw kappa4 itself moves) "
              "-> OPTION 3 (STOP)")
    res = dict(E=E, rows=rows, periodic_raw_ratio_mean=per_mean, periodic_raw_ratio_sd=per_sd,
               voronoi_raw_ratio_mean=vor_mean, raw_agrees=bool(agree), branch=branch)
    json.dump(res, open(os.path.join(AOUT, "e5a_rawk4_check.json"), "w"), indent=1)
    write_md(res)
    print(f"\nperiodic raw_k4 U/W = {per_mean:.3f} +/- {per_sd:.3f};  "
          f"voronoi raw_k4 U/W = {vor_mean:.3f}")
    print(f"BRANCH -> {branch}")
    return 0


def write_md(res):
    L = ["# e5a_rawk4_check.md -- S7 / E5(a) raw-kappa4 MSC diagnostic\n",
         "Decides whether the alt-MSC (Urban) voronoi collapse-steepening (all-inclusive slope "
         "-1.18 vs periodic-only -1.00; voronoi k_eff -> ~0.5) is a floor/deconvolution artifact "
         "(option 1) or a real MSC sensitivity in the stochastic family (option 3). "
         "**A floor artifact lives ONLY in the subtracted floor**, so the RAW (pre-subtraction) "
         "structured kappa4 must rescale UNIFORMLY across families; a real MSC effect moves the "
         "voronoi raw kappa4 differently from the periodic anchor. 500 MeV, existing data, no re-sim.\n",
         "## Urban/WentzelVI ratios per config\n",
         "| family | config | raw kappa4_struct U/W (95% CI) | floor U/W | dk4_simple U/W | "
         "dk4_allord U/W |", "|--|--|--|--:|--:|--:|"]
    for r in res["rows"]:
        ci = r["ratio_raw_ci"]
        L.append(f"| {r['family']} | {r['config']} | {r['ratio_raw_struct']:.3f} "
                 f"[{ci[0]:.3f}, {ci[1]:.3f}] | {r['ratio_floor']:.3f} | "
                 f"{r['ratio_dk4_simple']:.2f} | {r['ratio_dk4_allord']:.2f} |")
    L += [f"\n**Periodic raw kappa4 U/W = {res['periodic_raw_ratio_mean']:.3f} "
          f"+/- {res['periodic_raw_ratio_sd']:.3f}; voronoi raw kappa4 U/W = "
          f"{res['voronoi_raw_ratio_mean']:.3f}.**\n",
          "## Branch\n",
          f"Raw kappa4 rescales uniformly across families (voronoi raw ratio within the periodic "
          f"band): **{res['raw_agrees']}**.\n\n**{res['branch']}**\n",
          "### Reading\n",
          "- If the RAW structured kappa4 ratio (Urban/WentzelVI) is the SAME for voronoi and the "
          "periodic anchor, then Urban rescales the directly-measured 4th cumulant uniformly "
          "regardless of geometry -- the family-specific deconvolved-Delta_kappa4 deviation is "
          "injected by the SUBTRACTION (the floor), confirming the floor/deconvolution artifact. "
          "The voronoi is then DEMOTED to qualitative-only (its 2nd documented floor fragility); "
          "the headline exponent is the PERIODIC-only -1.00 (MSC-robust), and the prefactor MSC "
          "band 0.80-0.98 confirms the slope/prefactor split rather than threatening it.\n",
          "- If the raw voronoi ratio DEPARTS from the periodic anchor beyond stats, the stochastic "
          "structured kappa4 itself is MSC-sensitive -- a real SC2 tension; STOP and surface before "
          "the budget / S8.\n"]
    open(os.path.join(AOUT, "e5a_rawk4_check.md"), "w").write("\n".join(L) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
