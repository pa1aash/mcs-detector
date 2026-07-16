"""kink_stats.py -- S3 analysis helpers for the Geant4 kink-angle output.

Loads the projected MCS kink angles from a run's ROOT ntuple (uproot), computes
robust width estimators and acceptance-defined cumulants, and provides the
Highland projected-angle reference.

Why an acceptance for kappa4
----------------------------
The projected MCS angular distribution has a heavy single-scattering (Rutherford)
tail, so the raw 4th moment is dominated by rare large-angle events and is
acceptance-dependent. A real tracking telescope has finite acceptance; we mirror
that by computing all 4th-cumulant quantities within |theta| < ACCEPT_K * sigma_c
(sigma_c = 98%-central RMS), applied identically to every run. The kappa4
subtraction (Result 3) is then performed at matched acceptance, so the acceptance
cancels.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass

import numpy as np
import uproot

PROTON_MASS_MEV = 938.27208816
ACCEPT_K = 5.0  # angular acceptance for kappa4 in units of the 98%-central RMS


# --------------------------------------------------------------------------
# Proton kinematics and the Highland projected-angle width
# --------------------------------------------------------------------------
def proton_betacp_MeV(Ekin_MeV: float) -> float:
    """beta*(p c) in MeV for a proton of kinetic energy Ekin_MeV."""
    E = Ekin_MeV + PROTON_MASS_MEV
    pc = np.sqrt(E * E - PROTON_MASS_MEV ** 2)
    beta = pc / E
    return beta * pc


def proton_beta(Ekin_MeV: float) -> float:
    E = Ekin_MeV + PROTON_MASS_MEV
    pc = np.sqrt(E * E - PROTON_MASS_MEV ** 2)
    return pc / E


def highland_theta0(Ekin_MeV: float, t_mm: float, X0_mm: float, z: float = 1.0) -> float:
    """Highland projected RMS kink angle [rad] with the log correction.

    theta0 = (13.6/betacp) z sqrt(t/X0) [1 + 0.038 ln(t z^2/(X0 beta^2))].
    """
    betacp = proton_betacp_MeV(Ekin_MeV)
    beta = proton_beta(Ekin_MeV)
    tx0 = t_mm / X0_mm
    log_arg = tx0 * z ** 2 / (beta ** 2)
    return (13.6 / betacp) * z * np.sqrt(tx0) * (1.0 + 0.038 * np.log(log_arg))


# --------------------------------------------------------------------------
# Run loading
# --------------------------------------------------------------------------
@dataclass
class RunData:
    tag: str
    angles: np.ndarray  # pooled projected kink angles (thetax + thetay) [rad]
    meta: dict


def load_run(root_path: str) -> RunData:
    """Load a run's pooled projected kink angles and its metadata sidecar."""
    with uproot.open(root_path) as f:
        t = f["kinks"]
        tx = t["thetax"].array(library="np")
        ty = t["thetay"].array(library="np")
    angles = np.concatenate([tx, ty])
    angles = angles[np.isfinite(angles)]
    meta_path = root_path + ".meta.json"
    meta = {}
    if os.path.exists(meta_path):
        with open(meta_path) as fh:
            meta = json.load(fh)
    tag = os.path.splitext(os.path.basename(root_path))[0]
    return RunData(tag=tag, angles=angles, meta=meta)


# --------------------------------------------------------------------------
# Width estimators
# --------------------------------------------------------------------------
def rms(a: np.ndarray) -> float:
    return float(np.std(a, ddof=1))


def central_rms(a: np.ndarray, frac: float = 0.98) -> float:
    """RMS of the central `frac` of the (zero-centred) projected angles.

    Note: for a Gaussian this is biased *low* by the truncation (e.g. ~0.94 sigma
    at frac=0.98). Use `core_sigma` as the unbiased Highland comparison; this is
    kept only for the kappa4 acceptance window and as a reported cross-check.
    """
    a0 = a - np.median(a)
    lo, hi = np.percentile(a0, [50 * (1 - frac), 100 - 50 * (1 - frac)])
    core = a0[(a0 >= lo) & (a0 <= hi)]
    return float(np.std(core, ddof=1))


def core_sigma(a: np.ndarray) -> float:
    """Tail-robust Gaussian-core width = half the central 68.27% interval.

    sigma = (P84.135 - P15.865)/2. Equals sigma exactly for a Gaussian, is set by
    the core (insensitive to the single-scattering tail and to truncation bias),
    and is therefore the correct comparison to the Highland theta0 (which is a
    Gaussian fit to the central ~98%).
    """
    lo, hi = np.percentile(a, [15.865, 84.135])
    return float((hi - lo) / 2.0)


# --------------------------------------------------------------------------
# Acceptance-defined cumulants
# --------------------------------------------------------------------------
def acceptance_window(a: np.ndarray, k: float = ACCEPT_K) -> float:
    """Angular acceptance half-width = k * (98%-central RMS)."""
    return k * central_rms(a, 0.98)


def cumulants_in_window(a: np.ndarray, theta_acc: float) -> tuple[float, float]:
    """(kappa2, kappa4) of the zero-mean angles within |theta| < theta_acc."""
    a0 = a - np.mean(a)
    w = a0[np.abs(a0) < theta_acc]
    mu2 = np.mean(w ** 2)
    mu4 = np.mean(w ** 4)
    k2 = mu2
    k4 = mu4 - 3.0 * mu2 ** 2
    return float(k2), float(k4)


def bootstrap_kappa4(a: np.ndarray, theta_acc: float, n_boot: int = 400,
                     seed: int = 1234) -> tuple[float, float, float]:
    """Bootstrap (median, lo95, hi95) of kappa4 at fixed acceptance."""
    rng = np.random.default_rng(seed)
    n = a.size
    vals = np.empty(n_boot)
    for i in range(n_boot):
        s = a[rng.integers(0, n, n)]
        _, k4 = cumulants_in_window(s, theta_acc)
        vals[i] = k4
    return (float(np.median(vals)),
            float(np.percentile(vals, 2.5)),
            float(np.percentile(vals, 97.5)))
