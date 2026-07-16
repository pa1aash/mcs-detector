#!/usr/bin/env python3
"""e0_break.py -- S5 E0 analysis (calibration + line-integral prediction + break scale).

Produces the pieces of E0 that do NOT require lattice Geant4 transport:
  (T2) 100 MeV solid calibration computed IDENTICALLY to S3 (fixed window
       W(E)=5*sigma_core(16 mm)); emits W(100), sigma_core@4mm@100, a(100).
  (T3) The line-integral (straight-chord) prediction, the null the break departs
       from:  Delta_kappa4_LI(topology, E) = 3 a(E)^2 Var(t)  [Result 2/3],
       with the S3 D4 floor amplitude renormalisation carried.  Var(t) is the
       S4 ray-traced, CELL-INDEPENDENT chord variance (data/geom_stats/*.npz),
       so the null is FLAT in cell size -- the break is where transport departs.
  (T5-analytic) The predicted break scale: lateral MCS displacement
       y_rms(E) = (L/sqrt(3)) * theta0(E, <t>=fL) for the areal-matched solid.
       The line-integral picture fails when cell <~ y_rms (proton wanders across
       many cells, averaging out the chord-length fluctuation that drives Dk4).

It does NOT run lattice Geant4 (the production DetectorConstruction is a
homogeneous G4Box; STL/CADMesh loading was never integrated -- see e0_break.md
"BLOCKER").  Delta_kappa4_G4 and the empirical residual are therefore left
unmeasured; this script writes the prediction + the bound.
"""
from __future__ import annotations
import glob, json, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "lib"))
import kink_stats as ks  # noqa: E402

ROOT = os.path.dirname(HERE)
RUNS = os.path.join(ROOT, "data", "runs")
GEOM = os.path.join(ROOT, "data", "geom_stats")
OUT = os.path.join(ROOT, "data", "analysis")
os.makedirs(OUT, exist_ok=True)

L_TARGET = 10.0        # mm, fixed target depth (S4 GUARD2 decision)
F_INFILL = 0.40        # break-probe infill
PRINT_LIMIT_FDM = 0.5  # mm, conservative FDM minimum feature
PRINT_LIMIT_SLA = 0.05 # mm, fine resin / SLA minimum feature
ENERGIES = [100, 200, 500, 1000]
THICKS = [3, 4, 8, 16]


def solid_path(E, T):
    p = os.path.join(RUNS, f"solid_E{E}_t{T}.root")
    return p if os.path.exists(p) else None


def linfit_origin(t, y):
    t = np.asarray(t, float); y = np.asarray(y, float)
    return float(np.sum(t * y) / np.sum(t * t))


def quadfit_origin(t, y):
    t = np.asarray(t, float); y = np.asarray(y, float)
    A = np.vstack([t, t ** 2]).T
    (b, d), *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(b), float(d)


def calibrate():
    """Per-energy: sigma_core@16, W=5*sigma_core@16, sigma_core@4, a (kappa2 slope),
    kappa_M'' (quad), d4frac=kM''/(6 a^2).  Computed exactly as S3."""
    cal = {}
    for E in ENERGIES:
        ref = solid_path(E, 16)
        if ref is None:
            continue
        sc16 = ks.core_sigma(ks.load_run(ref).angles)        # rad
        W = ks.ACCEPT_K * sc16                                # rad (fixed window)
        ts, k2s, k4s, sc4 = [], [], [], None
        for T in THICKS:
            p = solid_path(E, T)
            if p is None:
                continue
            a = ks.load_run(p).angles
            if T == 4:
                sc4 = ks.core_sigma(a)
            k2, k4 = ks.cumulants_in_window(a, W)            # within the FIXED window
            ts.append(T); k2s.append(k2); k4s.append(k4)
        a_pow = linfit_origin(ts, k2s)                        # rad^2/mm  (kappa2 = a t)
        _, d = quadfit_origin(ts, k4s)
        kM2 = 2.0 * d
        d4 = kM2 / (6.0 * a_pow ** 2)
        cal[E] = dict(sigma_core_16=sc16, W=W, sigma_core_4=sc4,
                      a_pow=a_pow, kM2=kM2, d4frac=d4,
                      betacp_MeV=ks.proton_betacp_MeV(E))
    return cal


def main():
    cal = calibrate()

    # ---- (T3) line-integral prediction, per (topology, E) ------------------
    topo = {}
    for name in ("rectilinear", "gyroid"):
        d = np.load(os.path.join(GEOM, f"{name}_f40.npz"))
        topo[name] = dict(var_t=float(d["var_t"]), l_int=float(d["l_int"]),
                          N_eff=float(d["N_eff_exact"]), f=float(d["f"]),
                          t_mean=float(d["t_mean"]))

    li_rows = []
    for name, g in topo.items():
        for E in ENERGIES:
            if E not in cal:
                continue
            a = cal[E]["a_pow"]
            dk4_lead = 3.0 * a ** 2 * g["var_t"]              # rad^4  (Result 3)
            dk4_d4 = dk4_lead * (1.0 + cal[E]["d4frac"])      # D4-renormalised
            li_rows.append(dict(topology=name, E_MeV=E,
                                var_t_mm2=g["var_t"], N_eff=g["N_eff"],
                                a_pow_rad2_per_mm=a,
                                Dk4_LI_leading_rad4=dk4_lead,
                                Dk4_LI_d4corr_rad4=dk4_d4))

    # ---- (T5-analytic) predicted break scale: lateral MCS displacement -----
    # theta0 at the areal-matched solid thickness <t> = f L (=4 mm).  Use the
    # measured sigma_core@4mm (tail-robust) as theta0; y_rms = (L/sqrt3) theta0.
    break_rows = []
    for E in ENERGIES:
        if E not in cal:
            continue
        th0 = cal[E]["sigma_core_4"]                          # rad, <t>=4mm
        y_rms = (L_TARGET / np.sqrt(3.0)) * th0               # mm
        break_rows.append(dict(E_MeV=E, theta0_4mm_mrad=th0 * 1e3,
                               y_rms_lateral_mm=y_rms,
                               cell_break_pred_mm=y_rms,
                               vs_FDM=y_rms / PRINT_LIMIT_FDM,
                               vs_SLA=y_rms / PRINT_LIMIT_SLA))

    # ---- write artifacts ---------------------------------------------------
    out = dict(
        meta=dict(L_mm=L_TARGET, f=F_INFILL, energies_MeV=ENERGIES,
                  print_limit_FDM_mm=PRINT_LIMIT_FDM,
                  print_limit_SLA_mm=PRINT_LIMIT_SLA,
                  physics_list="FTFP_BERT + EmStandard_opt4 (WentzelVI+SS), cut0.05 maxStep0.1",
                  residual_threshold="max(20%, 3 sigma_bootstrap)  [for the empirical pass, pending Geant4]"),
        calibration={str(E): cal[E] for E in cal},
        topology_geometry=topo,
        line_integral_prediction=li_rows,
        break_scale_prediction=break_rows,
        empirical_status="BLOCKED: production DetectorConstruction is homogeneous "
                         "G4Box only; lattice STL loading not integrated; campaign "
                         "compute unprovisioned. Dk4_G4 and the measured residual "
                         "are NOT computed here. See e0_break.md.",
    )
    with open(os.path.join(OUT, "e0_break.json"), "w") as f:
        json.dump(out, f, indent=2)

    # CSV of the LI prediction
    import csv
    with open(os.path.join(OUT, "e0_li_prediction.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(li_rows[0].keys()))
        w.writeheader(); w.writerows(li_rows)

    # console summary
    print("=== 100 MeV calibration (computed identically to S3) ===")
    for E in ENERGIES:
        if E not in cal:
            continue
        c = cal[E]
        print(f"  E={E:4d}  sigma_core@16={c['sigma_core_16']*1e3:6.3f} mrad  "
              f"W=5sig={c['W']*1e3:6.2f} mrad  sigma_core@4={c['sigma_core_4']*1e3:6.3f} "
              f"a={c['a_pow']:.4e} rad^2/mm  d4frac={c['d4frac']*100:+.1f}%")
    print("\n=== break-scale prediction (lateral MCS y_rms vs print limit) ===")
    for r in break_rows:
        print(f"  E={r['E_MeV']:4d}  theta0@4mm={r['theta0_4mm_mrad']:6.3f} mrad  "
              f"y_rms={r['y_rms_lateral_mm']*1e3:6.1f} um  "
              f"= {r['vs_FDM']:.3f}x FDM(0.5mm)  {r['vs_SLA']:.2f}x SLA(0.05mm)")
    print("\n=== line-integral prediction Dk4_LI (cell-INDEPENDENT null) ===")
    for r in li_rows:
        print(f"  {r['topology']:11s} E={r['E_MeV']:4d}  Var(t)={r['var_t_mm2']:6.2f} "
              f"N_eff={r['N_eff']:.2f}  Dk4_LI(lead)={r['Dk4_LI_leading_rad4']:.3e} "
              f"Dk4_LI(D4)={r['Dk4_LI_d4corr_rad4']:.3e} rad^4")
    print(f"\nwrote {OUT}/e0_break.json + e0_li_prediction.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
