#!/usr/bin/env python3
"""voxel_vs_tess.py -- S5.5 Task 5: kink-angle kappa2/kappa4 agreement between the
voxelised and tessellated representations of the SAME Voronoi block (200 MeV
protons). A common fixed angular window (5*sigma_core of the tessellated run) is
applied to both; agreement within bootstrap error => the voxelisation is
artifact-free (the S4 decision to voxelise the 461k-facet Voronoi is validated).
"""
from __future__ import annotations
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "lib"))
import kink_stats as ks  # noqa: E402


def main():
    vox = ks.load_run("/tmp/vox_p.root").angles
    tess = ks.load_run("/tmp/tess_p.root").angles
    W = ks.ACCEPT_K * ks.core_sigma(tess)        # common fixed window
    out = {}
    for name, a in (("tessellated", tess), ("voxel", vox)):
        k2, k4 = ks.cumulants_in_window(a, W)
        k4m, lo, hi = ks.bootstrap_kappa4(a, W)
        sc = ks.core_sigma(a)
        out[name] = dict(n=a.size, sigma_core_mrad=sc * 1e3, k2=k2, k4=k4m,
                         k4_lo=lo, k4_hi=hi)
    # agreement
    t, v = out["tessellated"], out["voxel"]
    dk2 = v["k2"] / t["k2"] - 1.0
    dk4 = v["k4"] / t["k4"] - 1.0
    # do the bootstrap CIs overlap?
    k4_overlap = not (v["k4_hi"] < t["k4_lo"] or t["k4_hi"] < v["k4_lo"])
    print(f"window W = {W*1e3:.2f} mrad (5*sigma_core tessellated)")
    for name in ("tessellated", "voxel"):
        o = out[name]
        print(f"  {name:12s}: N={o['n']} sigma_core={o['sigma_core_mrad']:.3f} mrad "
              f"k2={o['k2']:.3e} k4={o['k4']:.3e} [{o['k4_lo']:.2e},{o['k4_hi']:.2e}]")
    print(f"  d(k2)={dk2*100:+.2f}%  d(k4)={dk4*100:+.2f}%  "
          f"k4 bootstrap CIs overlap: {k4_overlap}")
    verdict = "PASS" if (abs(dk2) < 0.03 and k4_overlap) else "CHECK"
    print(f"  VOXEL-VS-TESSELLATED: {verdict}")
    return out, dict(dk2=dk2, dk4=dk4, k4_overlap=bool(k4_overlap),
                     window_mrad=W * 1e3, verdict=verdict)


if __name__ == "__main__":
    main()
