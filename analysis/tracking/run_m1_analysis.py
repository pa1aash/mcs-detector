#!/usr/bin/env python3
"""run_m1_analysis.py -- assemble the M1 reconstruction-impact results.

Reads the telescope runs in data/runs_m1/, computes the struct-vs-matched-slab observables
(tail-fraction ratios with bootstrap CIs, track loss, pulls, impact-parameter proxy),
demonstrates the target-out subtraction on the Si-plane variant, and writes
results/M1/{M1_results.json, T3_impact.md, f6_data.npz}. Frozen procedure: ANALYSIS_PLAN M1.
"""
from __future__ import annotations
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
import track_impact as ti  # noqa: E402

RUNS = os.path.join(ROOT, "data", "runs_m1")
OUT = os.path.join(ROOT, "results", "M1")
T_FL = 3.0  # matched mean material (mm)


def rp(tag):
    return os.path.join(RUNS, tag + ".root")


def one(tag, E):
    p = rp(tag)
    return ti.analyze(p, E=E, t_fL=T_FL) if os.path.exists(p) else None


def main():
    os.makedirs(OUT, exist_ok=True)
    results = {}
    for E in (200, 1000):
        slab = one(f"m1_slab_massless_pencil_E{E}", E)
        if slab is None:
            print(f"  slab E{E} missing -- skip"); continue
        tk_s, pas_s, sig_s = ti.pooled_raw(slab)
        for struct in ("rectilinear", "gyroid"):
            o = one(f"m1_{struct}_massless_pencil_E{E}", E)
            if o is None:
                continue
            tk, pas, sig = ti.pooled_raw(o)
            row = dict(E=E, struct=struct, n=int(tk.size),
                       sigma_var_ratio=float(sig / sig_s),      # Result 1 (~1)
                       track_loss=float(1 - pas.mean()),
                       track_loss_slab=float(1 - pas_s.mean()),
                       track_loss_ratio=float((1 - pas.mean()) / (1 - pas_s.mean())),
                       pull_core=float(o["x"]["pull_core"]),
                       pull_kurt_win=float(o["x"]["pull_kurt_win"]),
                       pull_tail3=float(o["x"]["pull_tail3"]),
                       pull_tail3_slab=float(slab["x"]["pull_tail3"]),
                       pull_tail3_ratio=float(o["x"]["pull_tail3"] /
                                              max(slab["x"]["pull_tail3"], 1e-9)),
                       ip_tail_ratio=float(o["x"]["ip_tail_3sig"] /
                                           max(slab["x"]["ip_tail_3sig"], 1e-9)))
            for k in (3, 4, 5):
                fs, lo_s, hi_s = ti.boot_tailfrac(tk_s, sig_s, k)
                f_, lo, hi = ti.boot_tailfrac(tk, sig, k)
                row[f"tail{k}_struct"] = f_
                row[f"tail{k}_slab"] = fs
                row[f"tail{k}_ratio"] = float(f_ / max(fs, 1e-9))
                row[f"tail{k}_ratio_lo"] = float(lo / max(hi_s, 1e-9))
                row[f"tail{k}_ratio_hi"] = float(hi / max(lo_s, 1e-9))
                row[f"tail{k}_pred"] = o["x"]["tail_pred"][k]
            results[f"{struct}_E{E}"] = row
            print(f"  {struct} E{E}: sigvar_ratio={row['sigma_var_ratio']:.3f} "
                  f"tail3x={row['tail3_ratio']:.2f} loss_ratio={row['track_loss_ratio']:.2f}")

    # Si-plane target-out subtraction demo (E200): kappa4(struct+Si) - kappa4(empty+Si)
    si = {}
    for tag in ("rectilinear", "slab", "empty"):
        o = one(f"m1_{tag}_si_pencil_E200", 200)
        if o is not None:
            tk, _, _ = ti.pooled_raw(o)
            # windowed kappa4 of the reconstructed (triplet) kink
            sys.path.insert(0, os.path.join(ROOT, "analysis", "lib")); import kink_stats as ks
            _, k4 = ks.cumulants_in_window(tk, 37.84e-3)
            si[tag] = float(k4)
    if si:
        results["si_subtraction_E200"] = dict(
            k4_struct_si=si.get("rectilinear"), k4_slab_si=si.get("slab"),
            k4_empty_si=si.get("empty"),
            k4_struct_sub=si.get("rectilinear", 0) - si.get("empty", 0),
            note="target-out subtraction removes the Si-telescope kappa4")

    # F6 data: rect vs slab triplet-kink at 200 MeV (money plot)
    r2 = one("m1_rectilinear_massless_pencil_E200", 200)
    s2 = slab if (slab := one("m1_slab_massless_pencil_E200", 200)) else None
    if r2 and s2:
        tkr, _, sr = ti.pooled_raw(r2); tks, _, ss = ti.pooled_raw(s2)
        np.savez(os.path.join(OUT, "f6_data.npz"),
                 tk_rect=tkr, tk_slab=tks, sig_rect=sr, sig_slab=ss,
                 tpla_rect=r2["tpla"] if "tpla" in r2 else np.array([]))

    json.dump(results, open(os.path.join(OUT, "M1_results.json"), "w"), indent=1)

    # T3 table
    L = ["# T3 -- M1 reconstruction impact (struct vs matched slab)\n",
         "Massless telescope, pencil beam, 5e6 protons/config. sigma-var ratio ~1 confirms "
         "Result 1 on the reconstructed kinks; tail ratios (bootstrap 68% CI) and track-loss "
         "ratio quantify the geometry-induced cost. Numbers written by run_m1_analysis.py.\n",
         "| target | E [MeV] | sig-ratio | tail>3s | tail>4s | track-loss ratio | pull>3 ratio |",
         "|--|--:|--:|--:|--:|--:|--:|"]
    for key, r in results.items():
        if "_E" not in key or "si_" in key:
            continue
        L.append(f"| {r['struct']} | {r['E']} | {r['sigma_var_ratio']:.3f} | "
                 f"{r['tail3_ratio']:.2f}x [{r['tail3_ratio_lo']:.2f},{r['tail3_ratio_hi']:.2f}] | "
                 f"{r['tail4_ratio']:.2f}x | {r['track_loss_ratio']:.2f}x "
                 f"({r['track_loss']:.3f} vs {r['track_loss_slab']:.3f}) | "
                 f"{r['pull_tail3_ratio']:.2f}x |")
    open(os.path.join(OUT, "T3_impact.md"), "w").write("\n".join(L) + "\n")
    print(f"wrote {OUT}/M1_results.json + T3_impact.md + f6_data.npz")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
