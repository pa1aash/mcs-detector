"""theory.py -- numerically-evaluable forms of the S2 MCS scale-mixture theory.

Encodes the verified derivation (verification log: THEORY_CHECK.md; typeset:
paper/sections/theory.tex). All quantities are leading order in t/X0 <~ 0.05.

Symbols
-------
f      : solid (infill) volume fraction, <chi>.
L      : target depth along the beam (length).
lam_c  : telegraph correlation length lambda_c = Ls*Lv/(Ls+Lv) (length).
l_int  : directional integral correlation length ell_int (length).
Var_t  : variance of the beam-direction line-integral t (length^2).
a      : Highland scattering power, angle^2 / length  [Eq. (a)].
N_eff  : effective cell count (dimensionless).
gamma2 : excess kurtosis kappa4 / kappa2^2 (dimensionless).

Core relations
--------------
RESULT 1 :  kappa2 = a * <t> = a * f * L                       (width invariance)
RESULT 2 :  Delta_kappa4 = 3 a^2 Var(t);  gamma2 = 3(1-f)/(f N_eff);
            N_eff = f(1-f)L^2/Var(t) = L/(2 l_int) [many-cell] = 1 [extruded]
RESULT 3 :  kappa4 additive; Delta_kappa4 = kappa4_struct - kappa4_solid = 3 a^2 Var(t)
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "highland_a",
    "kappa2",
    "t_mean",
    "var_t_markov",
    "l_int_markov",
    "var_t_extruded",
    "N_eff_exact",
    "N_eff_asymp",
    "gamma2",
    "gamma2_extruded",
    "Delta_kappa4",
    "delta_kappa4_curve",
    "fit_exponent",
]


# --------------------------------------------------------------------------
# Highland scattering power and the width-invariance lemma (RESULT 1)
# --------------------------------------------------------------------------
def highland_a(beta_c_p_MeV: float, X0: float, z: float = 1.0) -> float:
    """Highland scattering power a = (13.6 MeV / (beta c p))^2 z^2 / X0.

    Units: angle^2 / length (X0 in the same length units as t, L).
    `beta_c_p_MeV` is beta*(p c) in MeV; `z` is the projectile charge number.
    """
    return (13.6 / beta_c_p_MeV) ** 2 * z ** 2 / X0


def t_mean(f: float, L: float) -> float:
    """Mean line-integral <t> = f L  [Eq. (tmean)]."""
    return f * L


def kappa2(a: float, f: float, L: float) -> float:
    """RESULT 1: kappa2(theta) = a <t> = a f L (topology-independent)."""
    return a * t_mean(f, L)


# --------------------------------------------------------------------------
# Var(t): exponential-chord / telegraph (Markov) and extruded forms
# --------------------------------------------------------------------------
def var_t_markov(f: float, L: float, lam_c: float) -> float:
    """Exact finite-L variance for the binary telegraph (Markov) chord model:

        Var(t) = 2 f (1-f) [ L lam_c - lam_c^2 (1 - exp(-L/lam_c)) ]   [Eq. (varmarkov)]

    Limits: L >> lam_c -> 2 f(1-f) lam_c L (N_eff -> L/2lam_c);
            L << lam_c -> f(1-f) L^2       (N_eff -> 1, extruded limit).
    """
    return 2.0 * f * (1.0 - f) * (L * lam_c - lam_c ** 2 * (1.0 - np.exp(-L / lam_c)))


def l_int_markov(lam_c: float) -> float:
    """Integral correlation length for the telegraph model: l_int = lam_c."""
    return lam_c


def var_t_extruded(f: float, L: float) -> float:
    """Extruded / binary cell: t = L * Bernoulli(f) -> Var(t) = f(1-f) L^2."""
    return f * (1.0 - f) * L ** 2


# --------------------------------------------------------------------------
# N_eff: exact definition and many-cell asymptote (RESULT 2)
# --------------------------------------------------------------------------
def N_eff_exact(f: float, L: float, Var_t: float) -> float:
    """Exact effective cell count  N_eff = f(1-f) L^2 / Var(t)  [Eq. (neffexact)].

    Primary definition. Bridges extruded (N_eff=1) and many-cell regimes.
    """
    return f * (1.0 - f) * L ** 2 / Var_t


def N_eff_asymp(L: float, l_int: float) -> float:
    """Many-cell asymptote  N_eff = L / (2 l_int)  [Eq. (neffasymp)], valid L >> l_int."""
    return L / (2.0 * l_int)


# --------------------------------------------------------------------------
# Shape law: gamma2 and Delta_kappa4 (RESULT 2 / RESULT 3)
# --------------------------------------------------------------------------
def gamma2(f: float, N_eff: float) -> float:
    """Excess kurtosis shape law  gamma2 = 3 (1-f) / (f N_eff)  [Eq. (shapelaw)]."""
    return 3.0 * (1.0 - f) / (f * N_eff)


def gamma2_extruded(f: float) -> float:
    """Extruded limit (N_eff=1): gamma2 = 3 (1-f) / f."""
    return 3.0 * (1.0 - f) / f


def Delta_kappa4(a: float, Var_t: float) -> float:
    """Geometry-induced fourth-cumulant signal  Delta_kappa4 = 3 a^2 Var(t)  [Eq. (dk4result)]."""
    return 3.0 * a ** 2 * Var_t


def delta_kappa4_curve(a: float, f: float, L: float, N_eff):
    """Predicted Delta_kappa4 as a function of N_eff at fixed (a, f, L):

        Delta_kappa4 = 3 a^2 f (1-f) L^2 / N_eff      (theoretical log-log slope -1).

    Accepts a scalar or array-like N_eff; returns the matching type.
    """
    N_eff = np.asarray(N_eff, dtype=float)
    return 3.0 * a ** 2 * f * (1.0 - f) * L ** 2 / N_eff


# --------------------------------------------------------------------------
# Utility: log-log slope estimator (for the N_eff collapse exponent test)
# --------------------------------------------------------------------------
def fit_exponent(x, y) -> float:
    """Least-squares slope of log(y) vs log(x); returns the power-law exponent."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    slope, _ = np.polyfit(np.log(x), np.log(y), 1)
    return float(slope)
