"""transport_raytrace.py -- S5(rebuilt) Phase B: transport-aware ray-tracer.

Integrates a proton trajectory through the ANALYTIC material field chi(x,y,z)
(TPMS implicit functions / analytic struts) in the small-angle MCS diffusion
limit (valid t/X0 <~ 0.05). No facets, no voxels -> cell size unlimited, no
resolution/step/MSC artifact. Per proton, marching in z with step dz << cell:

    dtheta ~ N(0, a * chi(r,z) * dz)     (angular diffusion ONLY in solid)
    dr     = theta * dz                  (transverse drift)
    t_actual += chi(r,z) * dz            (solid path length along the wandering ray)

a = the S3-validated Highland scattering power per unit SOLID path length
(kappa2/t through-origin slope, rad^2/mm). Sampling many protons across the
transverse unit cell builds p(t_actual). The straight-chord p(t_straight) is the
same integral at FIXED transverse position (no drift). The break is where
wandering makes Var(t_actual) < Var(t_straight).

The PER-PATH-LENGTH physics (Gaussian core + Moliere/kappa_M tail) is folded in
separately by the scale mixture over t_actual (see analysis/e0_break.py), so this
tool supplies only the geometric sampling under realistic wandering.
"""
from __future__ import annotations
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import topologies as topo  # noqa: E402

L_DEFAULT = 10.0


def _tuned_param(name, cell, f):
    p, _ = topo.tune_analytic(name, cell, f, n=80, tol=8e-4)
    return p


def transport_trace(name, cell, f, a, L=L_DEFAULT, n_proton=20000,
                    steps_per_cell=20, rng=None, param=None):
    """March n_proton protons through the analytic lattice. Returns dict with
    t_actual (wandering) and t_straight (fixed-transverse) arrays + trajectory
    moments. a in rad^2/mm; transverse wander is 2-D (x,y)."""
    if rng is None:
        rng = np.random.default_rng(0)
    if param is None:
        param = _tuned_param(name, cell, f)
    dz = cell / steps_per_cell
    nz = int(round(L / dz))
    dz = L / nz

    # start positions uniform over one unit cell; angles zero
    x0 = rng.uniform(0.0, cell, n_proton)
    y0 = rng.uniform(0.0, cell, n_proton)
    x = x0.copy(); y = y0.copy()
    tx = np.zeros(n_proton); ty = np.zeros(n_proton)
    t_act = np.zeros(n_proton)
    t_str = np.zeros(n_proton)
    sa = np.sqrt(a * dz)
    z = 0.5 * dz
    for _ in range(nz):
        zc = np.full(n_proton, z)
        chi = topo.chi_analytic(name, x, y, zc, cell, param).astype(np.float64)
        chi_str = topo.chi_analytic(name, x0, y0, zc, cell, param).astype(np.float64)
        t_act += chi * dz
        t_str += chi_str * dz
        # angular diffusion only in solid (variance a*chi*dz per projected plane)
        m = np.sqrt(chi) * sa
        tx += m * rng.standard_normal(n_proton)
        ty += m * rng.standard_normal(n_proton)
        x += tx * dz
        y += ty * dz
        z += dz
    return dict(t_actual=t_act, t_straight=t_str, theta_x=tx, theta_y=ty,
                x=x - x0, y=y - y0, dz=dz, nz=nz, param=param)


def solid_limit_covariance(a, L=L_DEFAULT, steps=400, n=200000, rng=None):
    """Unit test (chi==1): integrate the bare MCS SDE and return the simulated
    Var(theta(L)), Var(y(L)), Cov(theta,y)(L) vs the analytic aL, aL^3/3, aL^2/2."""
    if rng is None:
        rng = np.random.default_rng(1)
    dz = L / steps
    th = np.zeros(n); y = np.zeros(n)
    sa = np.sqrt(a * dz)
    for _ in range(steps):
        th += sa * rng.standard_normal(n)
        y += th * dz
    return dict(
        var_theta=float(np.var(th, ddof=1)), var_theta_analytic=a * L,
        var_y=float(np.var(y, ddof=1)), var_y_analytic=a * L ** 3 / 3.0,
        cov=float(np.cov(th, y)[0, 1]), cov_analytic=a * L ** 2 / 2.0,
        y_rms=float(np.std(y, ddof=1)), y_rms_analytic=np.sqrt(a * L ** 3 / 3.0))


# S3-locked Highland scattering power per unit solid path length (rad^2/mm),
# = kappa2/t through-origin slope (STATE Locked-analysis-parameters / e0_break.py).
A_OF_P = {100: 1.7197e-05, 200: 3.9146e-06, 500: 7.2058e-07, 1000: 2.1943e-07}


if __name__ == "__main__":
    # quick self-check
    a = A_OF_P[200]
    r = transport_trace("rectilinear", 2.5, 0.40, a, n_proton=8000, rng=np.random.default_rng(0))
    va, vs = np.var(r["t_actual"], ddof=1), np.var(r["t_straight"], ddof=1)
    print(f"rect c2.5 200MeV: <t_act>={r['t_actual'].mean():.3f} Var_act={va:.3f} "
          f"Var_str={vs:.3f} ratio={va/vs:.3f}  y_rms={np.std(r['y']):.4f} mm")
