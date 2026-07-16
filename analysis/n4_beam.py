#!/usr/bin/env python3
"""n4_beam.py -- N4 beam systematics: how beam optics + momentum spread move the
geometry-induced Delta_kappa4, and why the struct-minus-control subtraction absorbs them.

Two beam imperfections, tested on anchored Geant4 samples (rect f0.30 @200 MeV struct +
its matched t=3 mm PLA solid control), analysis-only:

  1. Angular spread (spot size + divergence + hit resolution) -> convolve every measured
     angle with an independent Gaussian N(0, sigma_div^2). A Gaussian has kappa4 = 0 and
     cumulants add, so the geometry kappa4 is invariant under Gaussian beam optics (up to a
     small fixed-window leakage). Delta_kappa4 = kappa4(struct) - kappa4(control) is
     unchanged: BOTH samples get the same convolution.

  2. Momentum spread delta_p/p -> the Highland angle scales as 1/(beta c p), so a proton of
     momentum p0(1+delta) has all its angles scaled by ~ (1 - delta). A per-event global
     scale s_i ~ N(1, sigma_p^2) is a SCALE MIXTURE that ADDS kappa4 -- but it is common-mode
     (the same beam hits struct and control), so it cancels in Delta_kappa4. This is exactly
     the momentum-mixing channel of arXiv:2509.12800, here a nuisance the control removes.

Both are demonstrated numerically: the geometry signal Delta_kappa4 = kappa4(struct) -
kappa4(control) shifts by << its bootstrap CI under realistic sigma_div, sigma_p. Writes
results/N4/beam_sys.json. Run inside g4highland: python analysis/n4_beam.py
"""
from __future__ import annotations
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "analysis", "lib"))
import kink_stats as ks  # noqa

OUT = os.path.join(ROOT, "results", "N4")
M2 = os.path.join(ROOT, "data", "runs_m2")
STRUCT = os.path.join(M2, "m2_rectilinear_tilt0.root")   # rect f0.30 @200 (ANCHORED)
CONTROL = os.path.join(M2, "m2_solid_t3.root")           # matched t=fL=3 mm PLA solid
W = 37.84e-3                                              # locked W(200)


def dk4(theta_s, theta_c):
    """Geometry signal proxy = kappa4(struct) - kappa4(control) in the locked window."""
    _, k4s = ks.cumulants_in_window(theta_s, W)
    _, k4c = ks.cumulants_in_window(theta_c, W)
    return k4s - k4c, k4s, k4c


def boot_ci(theta_s, theta_c, n_boot=400, n_work=300000, seed=5):
    rng = np.random.default_rng(seed)
    ws = theta_s if theta_s.size <= n_work else rng.choice(theta_s, n_work, replace=False)
    wc = theta_c if theta_c.size <= n_work else rng.choice(theta_c, n_work, replace=False)
    out = np.empty(n_boot)
    for i in range(n_boot):
        out[i] = dk4(ws[rng.integers(0, ws.size, ws.size)],
                     wc[rng.integers(0, wc.size, wc.size)])[0]
    scale = np.sqrt(min(theta_s.size, n_work) / theta_s.size)  # rescale SE to full n
    return float(np.std(out, ddof=1) * scale)


def main():
    os.makedirs(OUT, exist_ok=True)
    ts = ks.load_run(STRUCT).angles
    tc = ks.load_run(CONTROL).angles
    rng = np.random.default_rng(11)
    base, k4s0, k4c0 = dk4(ts, tc)
    ci = boot_ci(ts, tc)
    print(f"baseline: dk4 = {base:.3e}  (kappa4 struct {k4s0:.3e}, control {k4c0:.3e}); "
          f"bootstrap SE = {ci:.3e} ({100*ci/base:.1f}%)")
    res = {"baseline_dk4": base, "dk4_se": ci, "window_mrad": W * 1e3,
           "struct": os.path.relpath(STRUCT, ROOT), "control": os.path.relpath(CONTROL, ROOT),
           "divergence": [], "momentum": []}

    # 1. Gaussian angular spread (spot/divergence/hit resolution), applied to BOTH samples
    for sig in (1e-3, 2e-3, 3e-3):        # 1, 2, 3 mrad
        d = dk4(ts + rng.normal(0, sig, ts.size), tc + rng.normal(0, sig, tc.size))[0]
        res["divergence"].append(dict(sigma_mrad=sig * 1e3, dk4=d,
                                      frac_shift=(d - base) / base, n_sigma=(d - base) / ci))
        print(f"  divergence sigma={sig*1e3:.0f} mrad: dk4={d:.3e}  "
              f"shift={100*(d-base)/base:+.1f}% ({(d-base)/ci:+.2f} sigma)")

    # 2. Momentum spread: per-event global scale s ~ N(1, sigma_p^2), common-mode to both
    for sp in (0.02, 0.05):               # 2%, 5% delta_p/p
        ss = 1.0 + rng.normal(0, sp, ts.size)
        sc = 1.0 + rng.normal(0, sp, tc.size)
        d = dk4(ts * ss, tc * sc)[0]
        # also show what momentum spread does to the RAW kappa4 (uncancelled), struct only
        _, k4s_p = ks.cumulants_in_window(ts * ss, W)
        res["momentum"].append(dict(sigma_p=sp, dk4=d, frac_shift=(d - base) / base,
                                    n_sigma=(d - base) / ci,
                                    k4_struct_raw_shift=(k4s_p - k4s0) / k4s0))
        print(f"  momentum dp/p={sp*100:.0f}%: dk4={d:.3e}  shift={100*(d-base)/base:+.1f}% "
              f"({(d-base)/ci:+.2f} sigma); raw kappa4(struct) alone moves "
              f"{100*(k4s_p-k4s0)/k4s0:+.1f}%")
    json.dump(res, open(os.path.join(OUT, "beam_sys.json"), "w"), indent=1)
    print(f"wrote {OUT}/beam_sys.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
