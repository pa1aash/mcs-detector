#!/usr/bin/env python3
"""w_sensitivity.py -- acceptance-window systematic (referee request): every
fourth-cumulant quantity is defined inside the fixed window W(E) = 5*sigma_core(16mm).
This check re-runs the ENTIRE chain (a_eff calibration, kappa_M floor fit, all-order
deconvolution) at 4, 5, and 6 sigma_core on the anchored 200 MeV configurations and
reports the shift of the deconvolved Delta_kappa4 and of gamma_2.
Writes data/analysis/w_sensitivity.json.
"""
from __future__ import annotations
import json, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "analysis"))
sys.path.insert(0, os.path.join(ROOT, "analysis", "lib"))
import kink_stats as ks   # noqa
import e2_analysis as e2  # noqa

E = 200
RUNS = os.path.join(ROOT, "data", "runs")
CONFIGS = ("camp_rectilinear_f40_E200", "camp_gyroid_f40_E200")
SOLID_T = (2, 3, 4, 5, 8, 16)


def chain_at(Wp, cache):
    """a_eff, floor, and per-config (dk4, gamma2) all evaluated inside window Wp."""
    a_vals, kt, kv = [], [0.0], [0.0]
    for t in SOLID_T:
        ang = cache[f"solid_{t}"]
        k2, k4 = ks.cumulants_in_window(ang, Wp)
        if t <= 5:
            a_vals.append(k2 / t)
        kt.append(float(t)); kv.append(k4)
    a = float(np.mean(a_vals))
    kt, kv = np.array(kt), np.array(kv)
    b, c = np.linalg.lstsq(np.vstack([kt, kt ** 2]).T, kv, rcond=None)[0]
    kM = lambda tt: b * np.clip(tt, 0, 16.0) + c * np.clip(tt, 0, 16.0) ** 2
    out = {}
    for tag in CONFIGS:
        ang, tpla = cache[tag]
        k2, k4 = ks.cumulants_in_window(ang, Wp)
        dk4 = k4 - float(np.mean(kM(tpla)))
        out[tag] = dict(dk4=float(dk4), gamma2=float(dk4 / k2 ** 2))
    return a, out


def main():
    import uproot
    cache = {}
    for t in SOLID_T:
        cache[f"solid_{t}"] = ks.load_run(
            os.path.join(RUNS, f"solid_E{E}_t{t}.root")).angles
    for tag in CONFIGS:
        with uproot.open(os.path.join(RUNS, tag + ".root")) as fr:
            tr = fr["kinks"]
            cache[tag] = (np.concatenate([tr["thetax"].array(library="np"),
                                          tr["thetay"].array(library="np")]),
                          tr["tpla"].array(library="np"))

    sc16 = ks.core_sigma(cache["solid_16"])
    res = {"E": E, "sigma_core_16": float(sc16), "windows": {}}
    base = None
    for k in (4.0, 5.0, 6.0):
        a, cfg = chain_at(k * sc16, cache)
        res["windows"][f"{k:g}sigma"] = dict(W=float(k * sc16), a_eff=a, configs=cfg)
        if k == 5.0:
            base = cfg
    for k in ("4sigma", "6sigma"):
        for tag in CONFIGS:
            d = res["windows"][k]["configs"][tag]
            d["dk4_shift_pct"] = round(
                100 * (d["dk4"] / base[tag]["dk4"] - 1), 2)
            d["gamma2_shift_pct"] = round(
                100 * (d["gamma2"] / base[tag]["gamma2"] - 1), 2)
            print(f"{k:7s} {tag:28s} dk4 shift {d['dk4_shift_pct']:+6.2f}%   "
                  f"gamma2 shift {d['gamma2_shift_pct']:+6.2f}%")
    json.dump(res, open(os.path.join(ROOT, "data", "analysis",
                                     "w_sensitivity.json"), "w"), indent=1)
    print("wrote data/analysis/w_sensitivity.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
