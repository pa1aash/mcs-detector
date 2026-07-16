#!/usr/bin/env python3
"""m4_carbon.py -- M4 material transfer: amorphous carbon (rho 1.7).

Confirms (i) width invariance kappa2(carbon lattice)/kappa2(carbon solid@fL) ~ 1 (Result 1
is material-independent), and (ii) the geometry Delta_kappa4 collapses onto C = 1/N_eff with
ONLY the scattering power a rescaled to carbon (geometry f, N_eff unchanged). a_eff and the
floor come from the M4 carbon solids. N_eff from the as-built geometry (reuses e2_analysis).
Writes results/M4/carbon.json.
"""
from __future__ import annotations
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "analysis", "lib"))
sys.path.insert(0, HERE)
import kink_stats as ks   # noqa
import e2_analysis as e2  # noqa  (geom_asbuilt, W, X0)

M4 = os.path.join(ROOT, "data", "runs_m4")
OUT = os.path.join(ROOT, "results", "M4")
W, L, E = 37.84e-3, 10.0, 200
CONFIGS = [("gyroid", 0.05), ("gyroid", 0.30), ("voronoi", 0.05), ("voronoi", 0.30)]


def k2k4(tag):
    a = ks.load_run(os.path.join(M4, tag + ".root")).angles
    return (*ks.cumulants_in_window(a, W), a)


def load_tpla(tag):
    import uproot
    with uproot.open(os.path.join(M4, tag + ".root")) as f:
        return f["kinks"]["tpla"].array(library="np")


def main():
    os.makedirs(OUT, exist_ok=True)
    # carbon floor + a_eff from the M4 carbon solids
    kt, kv = [0.0], [0.0]
    for t in (1, 2, 3, 4, 8, 16):    # D8: full-range floor (t up to 16 mm)
        k2, k4, _ = k2k4(f"m4_csolid_t{t}")
        kt.append(float(t)); kv.append(k4)
    kt, kv = np.array(kt), np.array(kv)
    b, c = np.linalg.lstsq(np.vstack([kt, kt ** 2]).T, kv, rcond=None)[0]
    tmax = kt.max()
    kM = lambda tt: b * np.clip(tt, 0, tmax) + c * np.clip(tt, 0, tmax) ** 2
    aeff = float(np.mean([k2k4(f"m4_csolid_t{t}")[0] / t for t in (1, 2, 3, 4)]))
    # solid kappa2 at the matched fL, for the width-invariance ratio
    def k2_solid_at(t):  # interpolate kappa2(t) ~ aeff*t (linear through origin)
        return aeff * t
    print(f"a_eff(carbon,200) = {aeff:.4e}  (PLA was ~3.88e-6; carbon scatters more)")

    rows = []
    for topo, f in CONFIGS:
        tag = f"m4_{topo}_f{int(round(f*100)):02d}"
        if not os.path.exists(os.path.join(M4, tag + ".root")):
            continue
        k2, k4, _ = k2k4(tag)
        tpla = load_tpla(tag)
        dk4 = float(k4 - np.mean(kM(tpla)))
        fb, var_t, neff, _ = e2.geom_asbuilt(topo, f)     # as-built f, Var(t), N_eff
        width_ratio = k2 / k2_solid_at(f * L)             # Result 1 (carbon)
        C = dk4 / (3 * aeff ** 2 * fb * (1 - fb) * L ** 2)
        rows.append(dict(topo=topo, f=f, f_built=float(fb), N_eff=float(neff),
                         k2=k2, dk4=dk4, C=float(C), C_theory=float(1 / neff),
                         width_ratio=float(width_ratio), k_prefactor=float(C * neff)))
        print(f"  {topo:8s} f={f:.2f}  width_ratio={width_ratio:.3f}  N_eff={neff:6.2f}  "
              f"dk4={dk4:.3e}  C*N_eff={C*neff:.2f}")
    json.dump({"a_eff_carbon": aeff, "rows": rows}, open(os.path.join(OUT, "carbon.json"), "w"),
              indent=1)
    wr = [r["width_ratio"] for r in rows]
    kpf = [r["k_prefactor"] for r in rows if r["N_eff"] < 50 and r["dk4"] > 0]
    print(f"\nwidth invariance (carbon): mean={np.mean(wr):.3f} range [{min(wr):.3f},{max(wr):.3f}]")
    print(f"collapse prefactor C*N_eff (carbon): {[round(k,2) for k in kpf]} "
          f"(~const + ~PLA value => geometry law material-transfers with rescaled a)")
    print(f"wrote {OUT}/carbon.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
