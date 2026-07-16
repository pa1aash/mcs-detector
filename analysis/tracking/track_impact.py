#!/usr/bin/env python3
"""track_impact.py -- M1 downstream reconstruction-impact analysis.

Reads the Geant4 "telescope" ntuple (6 plane hits + truth kink + tpla per proton),
applies a Gaussian hit resolution AT ANALYSIS TIME, and reconstructs each track with
a broken-line / Kalman-equivalent fit whose multiple-scattering process noise is taken
from the HOMOGENEOUS Highland model (the deliberate model error). Quality cuts and the
M1 observables (tail fractions, track loss, pulls, impact-parameter proxy) follow.

Everything is per projection (x, y decouple in the small-angle regime) and vectorised
over events. The fit design matrix is identical for every event (planes fixed), so its
normal-equation inverse is precomputed once.

Frozen procedure: see ANALYSIS_PLAN.md M1.
"""
from __future__ import annotations
import os, sys
import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import kink_stats as ks  # noqa: E402  (core_sigma, highland_theta0, proton kinematics)

PLANE_Z = np.array([-300., -200., -100., 100., 200., 300.])  # mm
Z_VERTEX = -400.0        # impact-parameter proxy plane (mm), upstream of plane 0
SIGMA_HIT = 0.005        # 5 um Gaussian hit resolution (mm)
CHI2_PVAL_CUT = 0.01     # track-quality cut
OUTLIER_NSIG = 3.0       # iterative hit-outlier rejection threshold


# --------------------------------------------------------------------------
def load_telescope(root_path):
    """Load the telescope ntuple -> dict of (N,6) hit arrays [mm] + kink/tpla."""
    import uproot
    with uproot.open(root_path) as f:
        a = f["telescope"].arrays(library="np")
    n = a["x0"].size
    X = np.stack([a[f"x{i}"] for i in range(6)], axis=1)   # (N,6)
    Y = np.stack([a[f"y{i}"] for i in range(6)], axis=1)
    return dict(X=X, Y=Y, tpla=a["tpla"], thetax=a["thetax"], thetay=a["thetay"], n=n)


def smear(hits, sigma_hit=SIGMA_HIT, rng=None):
    """Add a Gaussian hit resolution (analysis-time), keeping truth separate."""
    rng = rng or np.random.default_rng(0)
    return hits + rng.normal(0.0, sigma_hit, size=hits.shape)


# --------------------------------------------------------------------------
# Broken-line / Kalman-equivalent fit (one projection), vectorised over events.
#   model:  m(z) = x0 + s*z + k*max(z,0)      (kink k at the z=0 scattering node)
#   prior:  k ~ N(0, theta0^2)                (HOMOGENEOUS Highland process noise)
# --------------------------------------------------------------------------
def _design(z):
    dn = np.where(z > 0, z, 0.0)
    return np.stack([np.ones_like(z), z, dn], axis=1)      # (6,3)


def brokenline_fit(M, z, theta0, sigma_hit=SIGMA_HIT):
    """M: (N,6) measured positions. Returns dict with params, kink, chi2, pval, pulls.

    dof = 6 hits + 1 prior pseudo-measurement - 3 params = 4."""
    nz = z.shape[0]                                         # 6 normally; fewer after rejection
    A = _design(z)                                          # (nz,3)
    Aa = np.vstack([A, [0., 0., 1.]])                       # (nz+1,3): +prior row
    w = np.concatenate([np.full(nz, 1.0 / sigma_hit ** 2), [1.0 / theta0 ** 2]])
    AtWA_inv = np.linalg.inv(Aa.T @ (w[:, None] * Aa))      # (3,3), same all events
    N = M.shape[0]
    B = np.concatenate([M, np.zeros((N, 1))], axis=1)       # (N,nz+1), prior target 0
    rhs = (w[None, :] * B) @ Aa                             # (N,3)
    par = rhs @ AtWA_inv.T                                  # (N,3): x0, s, k
    pred = par @ A.T                                        # (N,nz) at the hits
    resid = M - pred
    pulls = resid / sigma_hit
    dof = max(1, nz + 1 - 3)                                # nz hits + 1 prior - 3 params
    chi2 = np.sum(pulls ** 2, axis=1) + (par[:, 2] / theta0) ** 2
    pval = stats.chi2.sf(chi2, df=dof)
    return dict(x0=par[:, 0], s=par[:, 1], kink=par[:, 2],
                chi2=chi2, pval=pval, pulls=pulls)


def triplet_kink(M, z):
    """Model-independent kink: downstream slope - upstream slope (straight lines)."""
    up, dn = z < 0, z > 0
    def slope(mask):
        zz = z[mask]
        A = np.stack([np.ones_like(zz), zz], axis=1)
        binv = np.linalg.inv(A.T @ A) @ A.T                 # (2, n_mask)
        return (M[:, mask] @ binv.T)[:, 1]                  # slope column
    return slope(dn) - slope(up)


# --------------------------------------------------------------------------
# Iterative 3-sigma hit-outlier rejection (frozen). Refit only flagged tracks.
# --------------------------------------------------------------------------
def outlier_reject(M, z, theta0, sigma_hit=SIGMA_HIT, nsig=OUTLIER_NSIG, max_iter=3):
    """Iterative 3-sigma hit-outlier rejection, VECTORISED by keep-pattern. Each iteration:
    tracks with a kept-hit pull > nsig (and >4 hits) drop their worst hit; tracks sharing a
    keep-pattern are refit together (one brokenline_fit per distinct pattern, not per track).
    Returns (fit_dict, n_rejected, keep)."""
    N, nh = M.shape
    keep = np.ones((N, nh), bool)
    fit = brokenline_fit(M, z, theta0, sigma_hit)          # full-hit fit (vectorised)
    x0, s, kink = fit["x0"].copy(), fit["s"].copy(), fit["kink"].copy()
    chi2, pval, pulls = fit["chi2"].copy(), fit["pval"].copy(), fit["pulls"].copy()
    nrej = np.zeros(N, int)
    bit = (1 << np.arange(nh))
    for _ in range(max_iter):
        pk = np.where(keep, np.abs(pulls), -np.inf)
        worst = pk.max(axis=1)
        flag = (worst > nsig) & (keep.sum(axis=1) > 4)
        if not flag.any():
            break
        idx = np.where(flag)[0]
        wcol = pk[idx].argmax(axis=1)
        keep[idx, wcol] = False
        nrej[idx] += 1
        pat = (keep[idx] * bit).sum(axis=1)                # bitmask of kept hits
        for p in np.unique(pat):
            sel = idx[pat == p]
            cols = np.where((p >> np.arange(nh)) & 1)[0]
            r = brokenline_fit(M[np.ix_(sel, cols)], z[cols], theta0, sigma_hit)
            x0[sel], s[sel], kink[sel] = r["x0"], r["s"], r["kink"]
            chi2[sel], pval[sel] = r["chi2"], r["pval"]
            pulls[np.ix_(sel, cols)] = r["pulls"]           # kept-column residual pulls
    return dict(x0=x0, s=s, kink=kink, chi2=chi2, pval=pval, pulls=pulls), nrej, keep


# --------------------------------------------------------------------------
# Scale-mixture tail prediction from the measured p(t): P(|theta|>x) analytic.
# --------------------------------------------------------------------------
def mixture_tail(tpla, a_scatter, x_over_sigma, sigma):
    """Exact Gaussian scale-mixture tail P(|theta| > k*sigma) from the measured tpla.
    theta|t ~ N(0, a*t); marginalise over the sample tpla."""
    t = np.asarray(tpla)
    t = t[t > 0]
    sig_t = np.sqrt(a_scatter * t)                          # per-proton width
    x = x_over_sigma * sigma
    return float(np.mean(2.0 * stats.norm.sf(x / sig_t)))


# --------------------------------------------------------------------------
def analyze(root_path, E, t_fL, seed=0):
    """Full M1 observable set for one telescope run. t_fL = mean material (mm) that
    sets the HOMOGENEOUS Highland process noise theta0 (the model error)."""
    d = load_telescope(root_path)
    rng = np.random.default_rng(seed)
    X0_PLA = 315.423                                        # mm (GetRadlen)
    theta0 = ks.highland_theta0(E, t_fL, X0_PLA)            # homogeneous process noise
    Xs, Ys = smear(d["X"], rng=rng), smear(d["Y"], rng=rng)

    W = {100: 79.85e-3, 200: 37.84e-3, 500: 16.22e-3, 1000: 8.95e-3}[E]
    a_scatter = theta0 ** 2 / t_fL                          # solid scattering power
    out = {}
    for proj, M in (("x", Xs), ("y", Ys)):
        # (a) tail fractions from the MODEL-INDEPENDENT triplet kink (no fit shrinkage),
        #     normalised by the VARIANCE-based sigma = sqrt(kappa2 in window) -- equal for
        #     struct and matched slab by Result 1, so the excess is honest.
        tk = triplet_kink(M, PLANE_Z)
        k2, _ = ks.cumulants_in_window(tk, W)
        sig = float(np.sqrt(k2))                            # RMS-equivalent (Result-1 sigma)
        tails = {k: float(np.mean(np.abs(tk) > k * sig)) for k in (3, 4, 5)}
        pred = {k: mixture_tail(d["tpla"], a_scatter, k, sig) for k in (3, 4, 5)}
        # (b) reconstruction cost from the broken-line fit with HOMOGENEOUS process noise.
        fit, nrej, keep = outlier_reject(M, PLANE_Z, theta0)
        passcut = fit["pval"] > CHI2_PVAL_CUT
        pulls_flat = fit["pulls"][keep]
        # Pull non-Gaussianity is a 4th-moment quantity -> the RAW kurtosis is dominated by
        # rare single-scattering pulls (esp. at high p, where theta0 is small and the prior
        # fights real kinks). Report ROBUST measures: a tail-insensitive pull width, the
        # excess kurtosis WITHIN a pull window, and the >3 pull-tail fraction (the clean
        # non-Gaussianity signature).
        pw = pulls_flat[np.abs(pulls_flat) < 8.0]
        out[proj] = dict(
            sigma_var=sig, sigma_core=float(ks.core_sigma(tk)), tails=tails, tail_pred=pred,
            track_loss=float(1.0 - passcut.mean()),
            rej_frac=float((nrej > 0).mean()),
            slope_bias=float(np.mean(fit["s"])),            # true upstream slope = 0
            intercept_bias=float(np.mean(fit["x0"])),
            pull_core=float(ks.core_sigma(pulls_flat)),     # tail-robust pull width (~1)
            pull_kurt_win=float(stats.kurtosis(pw, fisher=True)),  # windowed excess kurtosis
            pull_tail3=float(np.mean(np.abs(pulls_flat) > 3.0)),   # P(|pull|>3)
            ks_pval_uniform=float(stats.kstest(fit["pval"][passcut], "uniform").statistic),
            n=int(M.shape[0]))
        # impact-parameter proxy: downstream-only line extrapolated to Z_VERTEX
        dn = PLANE_Z > 0
        A = np.stack([np.ones(dn.sum()), PLANE_Z[dn]], axis=1)
        binv = np.linalg.inv(A.T @ A) @ A.T
        p = M[:, dn] @ binv.T                               # (N,2): x0_dn, slope_dn
        ip = p[:, 0] + p[:, 1] * Z_VERTEX                   # extrapolated position at vertex
        ip_sig = ks.core_sigma(ip)
        out[proj]["ip_tail_3sig"] = float(np.mean(np.abs(ip) > 3 * ip_sig))
        out[proj]["theta0_homog"] = theta0
        # raw per-event arrays for bootstrap (driver level)
        out[proj]["_tk"] = tk
        out[proj]["_sig"] = sig
        out[proj]["_pass"] = passcut
    return out


def pooled_raw(o):
    """Pool x+y raw arrays from an analyze() result for bootstrapping."""
    tk = np.concatenate([o["x"]["_tk"], o["y"]["_tk"]])
    pas = np.concatenate([o["x"]["_pass"], o["y"]["_pass"]])
    sig = 0.5 * (o["x"]["_sig"] + o["y"]["_sig"])
    return tk, pas, sig


def boot_tailfrac(tk, sig, k, n_boot=2000, seed=3, n_work=200000):
    """68% interval of the tail fraction P(|tk| > k*sig). Point estimate on the FULL
    sample; the bootstrap SE is taken on a capped working subsample (n_work, without
    replacement) and rescaled SE_full = SE_work * sqrt(n_work/n_full) -- a proportion's
    SE scales as 1/sqrt(n), so this is exact in expectation and keeps >=2000 resamples
    fast at 1e7 statistics (same device as e2_analysis._k4_pt_se)."""
    rng = np.random.default_rng(seed)
    n = tk.size
    a = np.abs(tk)
    p_full = float(np.mean(a > k * sig))
    work = a if n <= n_work else rng.choice(a, n_work, replace=False)
    m = work.size
    vals = np.empty(n_boot)
    thr = k * sig
    for i in range(n_boot):
        vals[i] = np.mean(work[rng.integers(0, m, m)] > thr)
    se = float(np.std(vals, ddof=1)) * np.sqrt(m / n)
    return p_full, p_full - se, p_full + se
