#!/usr/bin/env python3
"""n2_wandering.py -- N2: analytic MCS-blurred chord autocovariance -> Var(t_actual),
validated against the transport ray-tracer across a cell scan (the c_break region).

The straight-chord line-integral variance is a double integral of the along-beam
autocovariance of the material indicator:
    Var(t_str) = int_0^L int_0^L C_z(|z-z'|) dz dz'.
A wandering proton is transversely displaced by MCS, so at depths z, z' it samples the
field at points separated (transversely) by the relative displacement delta = y(z')-y(z),
which DECORRELATES the two samples. Model the field's transverse correlation as a length
l_perp; averaging over the Gaussian relative displacement (variance sigma_delta^2 per axis,
two axes) multiplies the autocovariance by
    R(sigma_delta) = 1 / (1 + 2 sigma_delta^2 / l_perp^2)^2      (2-D transverse blur),
giving the MCS-blurred (wandering) variance
    Var(t_act) = int_0^L int_0^L C_z(|z-z'|) R(sigma_delta(z,z')) dz dz'.
sigma_delta(z,z') is the standard MCS relative-transverse-displacement rms for a proton
diffusing with effective scattering power a_g = a*f (only the solid fraction scatters):
    Var(y(z)) = a_g z^3/3,  Cov(y(z),y(z')) = a_g (z^2 z'/2 - z^3/6)  (z<z'),
    sigma_delta^2(z,z') = Var(y(z)) + Var(y(z')) - 2 Cov(y(z),y(z')).

Inputs C_z(u) (along-beam autocovariance) and l_perp (transverse integral correlation
length) are read off the SAME analytic field the collapse uses (no new simulation). The
tool then compares Var(t_act) against transport_trace's Monte-Carlo Var(t_actual) across a
cell scan; <=15% -> c_break promotes to [PREDICTED-with-theory], else keep empirical (C3).

Run inside g4highland:  python analysis/n2_wandering.py
"""
from __future__ import annotations
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "geom"))
import topologies as topo          # noqa
import raytrace as rt              # noqa
import transport_raytrace as tr    # noqa

OUT = os.path.join(ROOT, "results", "N2")
L = 10.0
A_OF_P = tr.A_OF_P                 # per unit SOLID path length (rad^2/mm), S3-locked


def straight_field(name, cell, f, nxy=80):
    """As-built straight-chord field: chi (nxy^2 rays, nz), dz, param."""
    p, _ = topo.tune_analytic(name, cell, f, n=80, tol=8e-4)
    chi, dz = topo.ray_chi_analytic(name, cell, p, L, nxy=nxy, dz=cell / 40.0)
    return chi, dz, p


def l_perp_of(name, cell, f, p, dz):
    """Transverse integral correlation length l_perp (mm) of a mid-depth z-slice:
    average the 1-D transverse autocovariance over x- and y-lines, integral length."""
    nxy = 96
    g = np.linspace(0.5, nxy - 0.5, nxy) * (cell / nxy)
    z0 = np.full((nxy, nxy), 0.5 * L)
    X, Y = np.meshgrid(g, g, indexing="ij")
    sl = topo.chi_analytic(name, X, Y, z0, cell, p).astype(float)
    f_sl = sl.mean()
    Cx = rt.autocovariance(sl, f_sl, nxy // 2)           # along-x, averaged over y-rows
    Cy = rt.autocovariance(sl.T, f_sl, nxy // 2)         # along-y
    dxy = cell / nxy
    lx = rt.integral_corr_length(Cx, dxy)
    ly = rt.integral_corr_length(Cy, dxy)
    return float(0.5 * (lx + ly))


def var_t_blurred(C_z, dz, a_g, l_perp, blur=True):
    """Double-integral Var(t) with the MCS transverse-decorrelation factor R.
    C_z: along-beam autocovariance array C_z[k], k=0..K (lag k*dz)."""
    nz = int(round(L / dz))
    z = (np.arange(nz) + 0.5) * dz
    K = C_z.size - 1
    lag = np.abs(z[:, None] - z[None, :])                # (nz,nz) |z-z'|
    k = np.clip(np.round(lag / dz).astype(int), 0, K)    # z on a dz grid -> lags exact
    Cij = C_z[k]                                         # autocovariance at each pair
    if blur:
        zi, zj = np.minimum(z[:, None], z[None, :]), np.maximum(z[:, None], z[None, :])
        var_yi = a_g * zi ** 3 / 3.0
        var_yj = a_g * zj ** 3 / 3.0
        cov = a_g * (zi ** 2 * zj / 2.0 - zi ** 3 / 6.0)
        s2 = np.clip(var_yi + var_yj - 2.0 * cov, 0.0, None)   # sigma_delta^2 per axis
        R = 1.0 / (1.0 + 2.0 * s2 / l_perp ** 2) ** 2
    else:
        R = 1.0
    return float(np.sum(Cij * R) * dz * dz)


def run(name, E, cells, f=0.30):
    a = A_OF_P[E]
    a_g = a * f
    rows = []
    rng = np.random.default_rng(7)
    for cell in cells:
        chi, dz, p = straight_field(name, cell, f)
        fbuilt = float(chi.mean())
        C_z = rt.autocovariance(chi, fbuilt, int(round(L / dz)) - 1)
        lp = l_perp_of(name, cell, f, p, dz)
        vs_pred = var_t_blurred(C_z, dz, a_g, lp, blur=False)
        va_pred = var_t_blurred(C_z, dz, a_g, lp, blur=True)
        # transport MC ground truth
        t = tr.transport_trace(name, cell, f, a, L=L, n_proton=30000,
                               steps_per_cell=max(20, int(round(0.05 / (cell / 20)))),
                               rng=rng, param=p)
        vs_mc = float(np.var(t["t_straight"], ddof=1))
        va_mc = float(np.var(t["t_actual"], ddof=1))
        y_rms = float(np.std(t["y"], ddof=1))
        rows.append(dict(topo=name, E=E, cell=cell, f_built=fbuilt, l_perp=lp,
                         y_rms=y_rms, c_over_yrms=cell / y_rms,
                         var_str_mc=vs_mc, var_act_mc=va_mc, ratio_mc=va_mc / vs_mc,
                         var_str_pred=vs_pred, var_act_pred=va_pred,
                         ratio_pred=va_pred / vs_pred,
                         rel_err_act=abs(va_pred / va_mc - 1.0),
                         rel_err_ratio=abs((va_pred / vs_pred) / (va_mc / vs_mc) - 1.0)))
    return rows


def main():
    os.makedirs(OUT, exist_ok=True)
    cells = [0.15, 0.25, 0.4, 0.56, 0.8, 1.2, 2.5]        # brackets c_break (~0.3-0.6 mm)
    allrows = []
    for E in (100, 200):
        for name in ("rectilinear", "gyroid"):
            allrows += run(name, E, cells)
    # verdict: the wandering CORRECTION is the ratio va/vs; validate the predicted ratio
    # against the MC ratio where the correction actually bites (ratio_mc < 0.97).
    bite = [r for r in allrows if r["ratio_mc"] < 0.97]
    worst = max((r["rel_err_ratio"] for r in bite), default=0.0)
    med = float(np.median([r["rel_err_ratio"] for r in bite])) if bite else 0.0
    verdict = "PASS (<=15%) -> c_break [PREDICTED-with-theory]" if worst <= 0.15 \
        else "FAIL (>15%) -> keep empirical c_break (C3)"
    out = dict(model="Var(t_act)=int int C_z(|z-z'|)/(1+2 sigma_delta^2/l_perp^2)^2; "
               "sigma_delta from a_g=a*f MCS position covariance",
               n_bite=len(bite), worst_rel_err_ratio=worst, median_rel_err_ratio=med,
               verdict=verdict, rows=allrows)
    json.dump(out, open(os.path.join(OUT, "wandering.json"), "w"), indent=1)
    print(f"N2 wandering: {len(bite)} configs in the correction regime (ratio_mc<0.97); "
          f"predicted-vs-MC correction ratio: median {100*med:.1f}%, worst {100*worst:.1f}%")
    print(f"  verdict: {verdict}")
    for r in allrows:
        print(f"  {r['topo']:11s} {r['E']:4d} c={r['cell']:.2f} c/yrms={r['c_over_yrms']:5.1f} "
              f"ratio_mc={r['ratio_mc']:.3f} ratio_pred={r['ratio_pred']:.3f} "
              f"rel_err={100*r['rel_err_ratio']:.1f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
