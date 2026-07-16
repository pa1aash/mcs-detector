#!/usr/bin/env python3
"""n3_bayes.py -- N3 metrology: recovery of N_eff and f vs proton count (statistics,
NOT ML). The "feasibility teeth" of the fourth-cumulant channel.

The Gaussian scale mixture has EXACT sufficient statistics: kappa2 = a<t> and
kappa4 = 3 a^2 Var(t). So (f, N_eff) are recovered model-free from the two cumulants --
kappa2 fixes <t> = fL, kappa4 fixes Var(t) = f(1-f)L^2/N_eff -- with NO mixing-density
model to bias the estimate. The posterior is the bootstrap distribution of (f, N_eff).

REWORK (D9, refinement round): the previous version drew ONE capped 2e5-event bootstrap
and rescaled its SD by sqrt(m/N) -- the plotted points then followed the 1/sqrt(N) law
by construction. This version draws GENUINE bootstrap resamples at every plotted size
(N = 1e5, 1e6, 3e6 protons; 2N projected angles per resample from the run's 6e6-angle
pool), so the scaling is measured, not imposed; N = 1e7 is an explicit extrapolation
from the measured 3e6 point and is flagged as such in the output. The anchor is the
fresh seeded P7 run (the original m2 file was regenerable and absent), and the floor +
scattering power come from the SAME all-order machinery as the headline analysis
(e2_analysis.build_floor / a_eff_of_E; the old local floor was clipped at t=5 mm, the
D8 bug pattern, and carried a stale 2nd-order D4 renorm).

DEGENERACY / IDENTIFIABILITY: kappa2 and kappa4 are the only two handles, so (f, N_eff)
are identifiable but the TOPOLOGY CLASS is not (the paper's own universality law).
Writes results/N3/posterior.json.
"""
from __future__ import annotations
import csv, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "analysis"))
sys.path.insert(0, os.path.join(ROOT, "analysis", "lib"))
import kink_stats as ks   # noqa
import e2_analysis as e2  # noqa

OUT = os.path.join(ROOT, "results", "N3")
TAG = "camp_rectilinear_f30_E200"       # rect f0.30 @200 (ANCHORED, seeded P7 run)
DATA = os.path.join(ROOT, "data", "runs", TAG + ".root")
E, L = 200, 10.0
N_BOOT = 400
SIZES_MEASURED = (100_000, 1_000_000, 3_000_000)   # protons; 3e6 = the full run
N_EXTRAP = 10_000_000


def run_seed():
    with open(os.path.join(ROOT, "RUNLOG.csv")) as fh:
        for row in csv.reader(fh):
            if len(row) > 4 and row[0] == "P7_repro" and row[1] == TAG:
                return int(row[3])
    return None


def recover(angles, fmean, W, A):
    """(f, N_eff) from the two windowed cumulants; all-order floor subtraction."""
    k2, k4 = ks.cumulants_in_window(angles, W)
    dk4 = k4 - fmean
    f = np.clip(k2 / A / L, 1e-3, 0.999)
    var_t = dk4 / (3.0 * A ** 2)
    neff = f * (1 - f) * L ** 2 / var_t if var_t > 0 else np.inf
    return float(f), float(neff)


def boot_at_size(pool, N_protons, fmean, W, A, seed):
    """Genuine bootstrap at the plotted size: N_BOOT resamples of 2*N_protons angles."""
    rng = np.random.default_rng(seed)
    n = 2 * N_protons
    fs, nes = np.empty(N_BOOT), np.empty(N_BOOT)
    for i in range(N_BOOT):
        s = pool[rng.integers(0, pool.size, n)]
        fs[i], nes[i] = recover(s, fmean, W, A)
    ok = np.isfinite(nes)
    return (float(np.mean(fs)), float(np.std(fs, ddof=1)),
            float(np.mean(nes[ok])), float(np.std(nes[ok], ddof=1)))


def main():
    os.makedirs(OUT, exist_ok=True)
    W, A = e2.W[E], e2.a_eff_of_E(E)
    kM = e2.build_floor(E)

    import uproot
    with uproot.open(DATA) as fr:
        tree = fr["kinks"]
        pool = np.concatenate([tree["thetax"].array(library="np"),
                               tree["thetay"].array(library="np")])
        tpla = tree["tpla"].array(library="np")
    fmean = float(np.mean(kM(tpla)))
    f0, ne0 = recover(pool, fmean, W, A)
    print(f"data: {tpla.size:.0f} protons ({pool.size} angles), {TAG} "
          f"(design f=0.30).  a_eff={A:.4e}, <kappa_M>={fmean:.3e}.  "
          f"point estimate: f={f0:.3f}, N_eff={ne0:.2f}")

    res = []
    for k, N in enumerate(SIZES_MEASURED):
        fm, fsd, nm, nsd = boot_at_size(pool, N, fmean, W, A, seed=11 + k)
        res.append(dict(N=int(N), measured=True,
                        f=(round(fm, 5), round(fsd, 6)),
                        N_eff=(round(nm, 4), round(nsd, 5))))
        print(f"  N={N:.0e} (measured): N_eff = {nm:.2f} +/- {nsd:.3f} "
              f"({100*nsd/nm:.1f}%);  f = {fm:.4f} +/- {fsd:.5f} ({100*fsd/fm:.2f}%)")

    # N = 1e7: extrapolated from the measured full-run point by the sampling law,
    # and drawn as such (open marker) in the figure.
    base = res[-1]
    s = np.sqrt(SIZES_MEASURED[-1] / N_EXTRAP)
    res.append(dict(N=int(N_EXTRAP), measured=False,
                    f=(base["f"][0], round(base["f"][1] * s, 6)),
                    N_eff=(base["N_eff"][0], round(base["N_eff"][1] * s, 5))))
    print(f"  N={N_EXTRAP:.0e} (EXTRAPOLATED from 3e6 by sqrt-N): "
          f"N_eff +/- {res[-1]['N_eff'][1]:.3f} "
          f"({100*res[-1]['N_eff'][1]/base['N_eff'][0]:.1f}%);  "
          f"f +/- {res[-1]['f'][1]:.5f}")

    json.dump({"data": os.path.relpath(DATA, ROOT), "seed": run_seed(),
               "a_eff": A, "floor_mean": fmean, "n_boot": N_BOOT,
               "convention": "1 sigma bootstrap width; N counts protons "
                             "(2N projected angles); all-order e2 floor",
               "point": dict(f=f0, N_eff=ne0), "posteriors": res},
              open(os.path.join(OUT, "posterior.json"), "w"), indent=1)
    print(f"wrote {OUT}/posterior.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
