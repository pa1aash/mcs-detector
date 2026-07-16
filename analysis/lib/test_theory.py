"""Limit-tests for theory.py against the analytic results verified in THEORY_CHECK.md.

Run:  python -m pytest analysis/lib/test_theory.py -v
GATE 2 numerical conditions:
  * extruded:           Var_t = f(1-f)L^2  =>  N_eff = 1  =>  gamma2 = 3(1-f)/f
  * exponential L>>lc:  N_eff -> L/(2 lam_c)
  * exponential L<<lc:  N_eff -> 1  (extruded limit)
  * gamma2 vs N_eff log-log slope = -1
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))

import theory as th  # noqa: E402


# ---------------------------------------------------------------------------
# RESULT 1 -- width invariance: kappa2 = a f L, independent of topology
# ---------------------------------------------------------------------------
def test_result1_width_invariance_topology_independent():
    a, f, L = 2.5e-4, 0.40, 4.0
    k2 = th.kappa2(a, f, L)
    assert k2 == pytest.approx(a * f * L)
    # kappa2 depends on p(t) only through the mean -> same for extruded and Markov
    # (both share <t> = f L); verify it is literally a*f*L regardless of Var(t).
    assert th.kappa2(a, f, L) == pytest.approx(a * th.t_mean(f, L))


# ---------------------------------------------------------------------------
# Extruded / binary cell -- EXACT N_eff = 1, gamma2 = 3(1-f)/f
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("f", [0.1, 0.25, 0.40, 0.6, 0.9])
def test_extruded_var_and_neff_exact(f):
    L = 4.0
    Var_t = th.var_t_extruded(f, L)
    assert Var_t == pytest.approx(f * (1 - f) * L ** 2)
    N_eff = th.N_eff_exact(f, L, Var_t)
    assert N_eff == pytest.approx(1.0)


@pytest.mark.parametrize("f", [0.1, 0.25, 0.40, 0.6, 0.9])
def test_extruded_gamma2(f):
    L = 4.0
    Var_t = th.var_t_extruded(f, L)
    N_eff = th.N_eff_exact(f, L, Var_t)
    assert th.gamma2(f, N_eff) == pytest.approx(3 * (1 - f) / f)
    assert th.gamma2(f, N_eff) == pytest.approx(th.gamma2_extruded(f))


def test_extruded_neff1_recovers_floor_aside_shape():
    # N_eff -> 1 reproduces the bare scale-mixture excess kurtosis 3(1-f)/f
    # (the intrinsic kappa_M floor is a separate additive term, cancelled in Delta_kappa4).
    f = 0.40
    assert th.gamma2(f, 1.0) == pytest.approx(3 * (1 - f) / f)


# ---------------------------------------------------------------------------
# Exponential-chord / telegraph -- finite-L Var(t) and both limits
# ---------------------------------------------------------------------------
def test_markov_var_matches_closed_form():
    f, L, lam_c = 0.4, 4.0, 0.3
    expected = 2 * f * (1 - f) * (L * lam_c - lam_c ** 2 * (1 - np.exp(-L / lam_c)))
    assert th.var_t_markov(f, L, lam_c) == pytest.approx(expected)


def test_markov_many_cell_limit_neff_to_L_over_2lamc():
    # L >> lam_c  =>  N_eff -> L/(2 lam_c)
    f, L, lam_c = 0.4, 1000.0, 0.5
    Var_t = th.var_t_markov(f, L, lam_c)
    N_exact = th.N_eff_exact(f, L, Var_t)
    N_asymp = th.N_eff_asymp(L, th.l_int_markov(lam_c))
    assert N_asymp == pytest.approx(L / (2 * lam_c))
    # exact converges to the asymptote in the deep many-cell regime
    assert N_exact == pytest.approx(N_asymp, rel=2e-3)


def test_markov_thin_limit_neff_to_one():
    # L << lam_c  =>  Var(t) -> f(1-f)L^2,  N_eff -> 1 (extruded limit)
    f, L, lam_c = 0.4, 1.0e-3, 1.0e3
    Var_t = th.var_t_markov(f, L, lam_c)
    assert Var_t == pytest.approx(f * (1 - f) * L ** 2, rel=1e-3)
    assert th.N_eff_exact(f, L, Var_t) == pytest.approx(1.0, rel=1e-3)


def test_markov_interpolates_monotonically():
    # N_eff increases monotonically from 1 (thin) toward L/(2lam_c) (thick) as L/lam_c grows
    f, lam_c = 0.4, 1.0
    Ls = np.array([0.01, 0.1, 1.0, 10.0, 100.0])
    neffs = [th.N_eff_exact(f, L, th.var_t_markov(f, L, lam_c)) for L in Ls]
    assert all(np.diff(neffs) > 0)
    assert neffs[0] == pytest.approx(1.0, rel=1e-2)


# ---------------------------------------------------------------------------
# Shape-law exponent: gamma2 vs N_eff log-log slope = -1
# ---------------------------------------------------------------------------
def test_gamma2_vs_neff_loglog_slope_minus_one():
    f = 0.40
    N_eff = np.logspace(0, 3, 40)  # 1 .. 1000
    g = th.gamma2(f, N_eff)
    assert th.fit_exponent(N_eff, g) == pytest.approx(-1.0, abs=1e-9)


def test_delta_kappa4_vs_neff_loglog_slope_minus_one():
    a, f, L = 2.5e-4, 0.40, 4.0
    N_eff = np.logspace(0, 3, 40)
    dk4 = th.delta_kappa4_curve(a, f, L, N_eff)
    assert th.fit_exponent(N_eff, dk4) == pytest.approx(-1.0, abs=1e-9)


# ---------------------------------------------------------------------------
# Cross-consistency: Delta_kappa4 = 3 a^2 Var(t) and the curve agree;
# gamma2 from N_eff agrees with 3 Var(t)/<t>^2
# ---------------------------------------------------------------------------
def test_delta_kappa4_curve_consistent_with_var_form():
    a, f, L, lam_c = 2.5e-4, 0.4, 4.0, 0.3
    Var_t = th.var_t_markov(f, L, lam_c)
    N_eff = th.N_eff_exact(f, L, Var_t)
    assert th.delta_kappa4_curve(a, f, L, N_eff) == pytest.approx(th.Delta_kappa4(a, Var_t))


def test_gamma2_consistent_with_cv_definition():
    # gamma2 = 3 Var(t)/<t>^2 should equal gamma2(f, N_eff_exact)
    f, L, lam_c = 0.4, 4.0, 0.3
    Var_t = th.var_t_markov(f, L, lam_c)
    direct = 3 * Var_t / (f * L) ** 2
    N_eff = th.N_eff_exact(f, L, Var_t)
    assert th.gamma2(f, N_eff) == pytest.approx(direct)


def test_extruded_neff_disagrees_with_asymptote_factor_two():
    # Documents the verified factor-of-2 note: for the extruded cell l_int = L,
    # so the asymptote L/(2 l_int) = 1/2 (WRONG) while the exact N_eff = 1 (RIGHT).
    f, L = 0.4, 4.0
    Var_t = th.var_t_extruded(f, L)
    assert th.N_eff_exact(f, L, Var_t) == pytest.approx(1.0)
    assert th.N_eff_asymp(L, l_int=L) == pytest.approx(0.5)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
