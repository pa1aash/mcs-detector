#!/usr/bin/env python3
"""make_impact_hist.py -- derive the small, git-trackable figure input for fig_impact
from the (large, untracked) results/M1/f6_data.npz: normalised-kink histogram counts
in the fixed 0.1-sigma binning, the fitted sigmas, and the integer-k tail fractions.
Also computes p_void = P(t_pla == 0) for the M1 rectilinear target directly from the
pencil run's tpla branch, so the figure's void-channel annotation is artifact-backed
rather than a hard-coded constant.
Run after analysis/tracking/run_m1_analysis.py. Output: results/M1/f6_hist.json.
"""
from __future__ import annotations
import json, os
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NPZ = os.path.join(ROOT, "results", "M1", "f6_data.npz")
PENCIL = os.path.join(ROOT, "data", "runs_m1", "m1_rectilinear_massless_pencil_E200.root")
OUT = os.path.join(ROOT, "results", "M1", "f6_hist.json")

d = np.load(NPZ)
tkr, tks = d["tk_rect"], d["tk_slab"]
sr, ss = float(d["sig_rect"]), float(d["sig_slab"])
bins = np.linspace(-6, 6, 121)          # 0.1-sigma bins (frozen)
out = {"bins": bins.tolist(), "sig_rect": sr, "sig_slab": ss,
       "n_rect": int(tkr.size), "n_slab": int(tks.size)}
for key, a, s in (("rect", tkr, sr), ("slab", tks, ss)):
    x = a / s
    h, _ = np.histogram(x, bins=bins)
    out[f"counts_{key}"] = h.tolist()
    out[f"tail_{key}"] = {str(k): float(np.mean(np.abs(x) > k)) for k in (2, 3, 4, 5)}

if os.path.exists(PENCIL):
    import uproot
    tpla = uproot.open(PENCIL)["kinks"]["tpla"].array(library="np")
    out["p_void"] = float(np.mean(tpla == 0.0))
    out["p_void_source"] = os.path.relpath(PENCIL, ROOT)
else:
    print("WARNING: pencil run absent; p_void not written:", PENCIL)

json.dump(out, open(OUT, "w"))
print("wrote", OUT, "p_void=%.4f" % out.get("p_void", float("nan")))
