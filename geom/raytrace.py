"""raytrace.py -- S4 beam-direction ray statistics for the line-integral picture.

Given the binary material indicator chi(x,y,z) sampled along +z for many
transverse ray positions (a 2-D array, one row per ray, columns = z samples),
compute everything the N_eff theory needs:

  * p(t): distribution of the material line-integral t = integral chi dz, and its
    first four cumulants.
  * C(u): along-z autocovariance of chi, averaged over rays.
  * l_int: directional integral correlation length = (1/C(0)) * int_0^inf C(u) du.
  * Var(t), N_eff_exact = f(1-f)L^2/Var(t), N_eff_asymp = L/(2 l_int).

The same estimator is validated against the analytic telegraph (theory.py B10)
before it is trusted on any printed geometry (see validate_estimator.py).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class RayStats:
    f: float                 # solid fraction <chi>
    L: float                 # target depth (mm)
    t_mean: float            # <t>  (= f L)
    var_t: float             # Var(t)
    kappa: tuple             # (k1,k2,k3,k4) cumulants of p(t)
    l_int: float             # integral correlation length (mm)
    N_eff_exact: float       # f(1-f)L^2 / Var(t)
    N_eff_asymp: float       # L / (2 l_int)
    C: np.ndarray = field(repr=False)   # autocovariance C(u), u = 0..len-1 in dz
    u: np.ndarray = field(repr=False)   # lag axis (mm)


def cumulants(x: np.ndarray) -> tuple:
    """First four cumulants of samples x (k1 mean, k2 var, k3, k4=mu4-3 mu2^2)."""
    x = np.asarray(x, float)
    m = x.mean()
    d = x - m
    mu2 = np.mean(d ** 2)
    mu3 = np.mean(d ** 3)
    mu4 = np.mean(d ** 4)
    k1, k2, k3 = m, mu2, mu3
    k4 = mu4 - 3.0 * mu2 ** 2
    return (float(k1), float(k2), float(k3), float(k4))


def autocovariance(chi: np.ndarray, f: float, max_lag: int) -> np.ndarray:
    """Unbiased along-z autocovariance C[k], k=0..max_lag, averaged over rays.

    C[k] = < (chi[z]-f)(chi[z+k]-f) >  over all rays and valid z (FFT per row).
    Unbiased: lag k divided by the (Nz-k) overlapping z-pairs (x n_rays), NOT by
    Nz. The biased ÷Nz form imposes a Bartlett triangular window C_true*(1-k/Nz)
    that biases the integral correlation length low by ~l_int/L -- unacceptable
    for the l_int = int C(u) du estimate. The (Nz-k) divisor removes that window;
    large-lag noise is bounded because n_rays pairs remain at every lag.
    """
    a = chi.astype(float) - f
    n_rays, nz = a.shape
    nfft = 1 << int(np.ceil(np.log2(2 * nz)))
    F = np.fft.rfft(a, n=nfft, axis=1)
    ac = np.fft.irfft(F * np.conj(F), n=nfft, axis=1)[:, :nz]
    ac = ac.sum(axis=0)                            # sum over rays and z-pairs
    counts = (nz - np.arange(nz)) * n_rays         # unbiased # of pairs per lag
    ac = ac / counts
    return ac[:max_lag + 1]


def integral_corr_length(C: np.ndarray, dz: float) -> float:
    """l_int = (1/C(0)) int_0^inf C(u) du, left-Riemann sum (-> lambda_c for a
    telegraph field as dz->0). Integrates the full supplied C (caller truncates)."""
    if C[0] <= 0:
        return 0.0
    return float(C.sum() * dz / C[0])


def stats_from_chi(chi: np.ndarray, dz: float, L: float,
                   corr_frac: float = 0.5) -> RayStats:
    """Full ray statistics from a (n_rays, nz) indicator array.

    corr_frac: integrate C(u) over the first `corr_frac` of the z-range (the tail
    is pair-starved and noisy; 0.5 captures the decay for all geometries here).
    """
    n_rays, nz = chi.shape
    f = float(chi.mean())
    t = chi.sum(axis=1) * dz                       # line-integral per ray
    k = cumulants(t)
    var_t = k[1]
    max_lag = max(1, int(corr_frac * nz))
    C = autocovariance(chi, f, max_lag)
    u = np.arange(C.size) * dz
    l_int = integral_corr_length(C, dz)
    fL2 = f * (1.0 - f) * L ** 2
    N_eff_exact = fL2 / var_t if var_t > 0 else np.inf
    N_eff_asymp = L / (2.0 * l_int) if l_int > 0 else np.inf
    return RayStats(f=f, L=L, t_mean=k[0], var_t=var_t, kappa=k, l_int=l_int,
                    N_eff_exact=N_eff_exact, N_eff_asymp=N_eff_asymp, C=C, u=u)


# --------------------------------------------------------------------------
# Synthetic telegraph (two-state Markov) field -- estimator validation harness
# --------------------------------------------------------------------------
def make_telegraph(f: float, lam_c: float, L: float, dz: float, n_rays: int,
                   rng: np.random.Generator) -> np.ndarray:
    """Sample n_rays independent stationary telegraph realisations chi(z) on a
    z-grid of spacing dz over [0,L].

    Mean chords: Lambda_s = lam_c/(1-f), Lambda_v = lam_c/f  (so f=Ls/(Ls+Lv),
    lam_c = Ls Lv/(Ls+Lv)). Stationary start: solid w.p. f. Exponential segment
    lengths -> discrete switch probability p_switch = dz/Lambda per step.
    """
    nz = int(round(L / dz))
    Lambda_s = lam_c / (1.0 - f)     # mean solid chord
    Lambda_v = lam_c / f             # mean void chord
    ps = dz / Lambda_s               # P(switch out of solid) per step
    pv = dz / Lambda_v               # P(switch out of void) per step
    chi = np.empty((n_rays, nz), dtype=np.uint8)
    state = (rng.random(n_rays) < f)                 # True = solid
    u = rng.random((n_rays, nz))
    for j in range(nz):
        chi[:, j] = state
        sw = np.where(state, u[:, j] < ps, u[:, j] < pv)
        state = np.where(sw, ~state, state)
    return chi
