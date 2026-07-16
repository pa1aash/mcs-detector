#!/usr/bin/env python3
"""Model robustness of the fourth-cumulant results.

  (a) matched four-configuration periodic-subset collapse slopes under the locked
      WentzelVI and alternative Urban MSC models.  Both models are recomputed from
      the same rows in e5_msc_systematic.json; error bars are honest two-sided 95%
      Student-t OLS intervals (n=4, two residual degrees of freedom).
  (b) the directly simulated foam-scale excess kurtosis at the 0.2 mm rectilinear
      cell under the two MSC models.

The scenario/configuration magnitude budget (Urban-MSC and fabrication shifts of
the deconvolved Delta kappa_4) is tabulated in the discussion (Table tab:systematics),
where the wide dynamic range reads more clearly than on a single linear axis.

Reads data/analysis/e5_msc_systematic.json.
"""
from __future__ import annotations

import re

import numpy as np
from scipy.stats import t as student_t

import figstyle as fs

fs.set_style()
import matplotlib.pyplot as plt


E5 = fs.load("e5_msc_systematic.json")

MSC = "#6f4c8b"          # alternative Urban MSC model (distinct from topology colours)
MATCHED = {
    "rectilinear_f40",
    "gyroid_f20",
    "gyroid_f40",
    "schwarzp_f40",
}
L_MM = 10.0


def fill_from_config(config):
    match = re.search(r"_f(\d+)$", config)
    if not match:
        raise ValueError(f"cannot parse fill fraction from {config}")
    return int(match.group(1)) / 100.0


def matched_periodic_slope(energy, model):
    """OLS slope and 95% Student-t interval for the same four E5 configurations.

    The ordinate is the normalized collapse variable
      C = Delta kappa4 / [3 a_eff^2 f(1-f) L^2].
    The common a_eff factor does not affect a per-energy slope, but retaining it
    here makes the calculation identical to the plotted collapse definition.
    """
    rows = [
        row for row in E5["collapse"]["energies"][str(energy)]
        if row["config"] in MATCHED
    ]
    if {row["config"] for row in rows} != MATCHED:
        missing = MATCHED - {row["config"] for row in rows}
        raise ValueError(f"E={energy}: missing matched E5 configurations: {missing}")

    a_eff = float(E5["energies"][str(energy)][f"a_eff_{model}"])
    x, y = [], []
    for row in rows:
        f = fill_from_config(row["config"])
        dk4 = float(row[f"dk4_{model}"])
        C = dk4 / (3.0 * a_eff**2 * f * (1.0 - f) * L_MM**2)
        x.append(np.log(float(row["N_eff"])))
        y.append(np.log(C))

    x = np.asarray(x)
    y = np.asarray(y)
    A = np.column_stack([x, np.ones_like(x)])
    slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
    residual = y - (slope * x + intercept)
    dof = len(x) - 2
    slope_se = np.sqrt(
        np.sum(residual**2) / dof / np.sum((x - np.mean(x)) ** 2)
    )
    half_width = student_t.ppf(0.975, dof) * slope_se
    return float(slope), float(slope - half_width), float(slope + half_width), len(x)


def main():
    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=(fs.FULL, fs.FULL * 0.44),
        gridspec_kw={"width_ratios": (1.12, 0.88)},
    )

    encodings = {
        "locked": dict(label="WentzelVI", offset=-0.11, marker="o", colour=fs.HERO,
                       mfc=fs.HERO),
        "alt": dict(label="Urban", offset=0.11, marker="s", colour=MSC, mfc="white"),
    }

    # ------------------------------------------------------------------ panel a
    # Same four periodic configurations under both models at each energy.  With
    # n=4, the honest t intervals are necessarily broad; that is part of the result.
    energies = (200, 500)
    xbase = np.arange(len(energies), dtype=float)
    for model, style in encodings.items():
        vals = [matched_periodic_slope(E, model) for E in energies]
        centre = np.array([v[0] for v in vals])
        lower = np.array([v[1] for v in vals])
        upper = np.array([v[2] for v in vals])
        axA.errorbar(
            xbase + style["offset"], centre,
            yerr=[centre - lower, upper - centre],
            fmt=style["marker"], ms=6.2, mfc=style["mfc"],
            mec=style["colour"], mew=1.2, ecolor=style["colour"],
            elinewidth=1.0, capsize=3.0, capthick=0.9, ls="", zorder=4,
            label=style["label"],
        )

    axA.axhline(-1.0, color=fs.MUTE, lw=1.0, ls=(0, (4, 3)), zorder=1)
    axA.text(0.50, -0.955, r"derived $-1$", fontsize=8.5, color=fs.MUTE,
             ha="center", va="bottom")
    axA.set_xlim(-0.42, 1.42)
    axA.set_ylim(-1.82, -0.02)
    axA.set_xticks(xbase)
    axA.set_xticklabels(["200", "500"])
    axA.set_yticks([-1.8, -1.4, -1.0, -0.6, -0.2])
    axA.set_xlabel(r"proton energy [\si{MeV}]")
    axA.set_ylabel(r"collapse exponent")
    axA.legend(
        loc="upper center", bbox_to_anchor=(0.5, 1.00), ncol=2,
        fontsize=8.5, frameon=False,
        handletextpad=0.4, columnspacing=1.0, borderaxespad=0.25,
    )
    fs.despine(axA)
    fs.panel_title(axA, "a", "Matched periodic-subset scaling")

    # ------------------------------------------------------------------ panel b
    gamma_locked = float(E5["foam"]["locked_gamma2_g4"])
    gamma_urban = float(E5["foam"]["measured_gamma2_alt"])
    foam_rows = [
        (1, gamma_locked, encodings["locked"]),
        (0, gamma_urban, encodings["alt"]),
    ]
    for y, value, style in foam_rows:
        line_style = "-" if style["label"] == "WentzelVI" else (0, (3, 2))
        axB.plot([0, value], [y, y], color=style["colour"], lw=1.8,
                 ls=line_style, zorder=2)
        axB.plot(
            value, y, marker=style["marker"], ms=7.2, mfc=style["mfc"],
            mec=style["colour"], mew=1.2, ls="", zorder=4,
        )
        axB.text(value + 0.10, y, f"{value:.2f}", fontsize=9.0,
                 color=style["colour"], ha="left", va="center")

    axB.axvline(0.0, color=fs.MUTE, lw=0.8, zorder=1)
    axB.set_xlim(-0.05, 2.35)
    axB.set_ylim(-0.65, 1.55)
    axB.set_xticks([0, 1, 2])
    axB.set_yticks([1, 0])
    axB.set_yticklabels(["WentzelVI", "Urban"], fontsize=8.7)
    axB.set_xlabel(r"foam-scale excess kurtosis $\gamma_2$")
    axB.text(
        0.03, -0.48, r"homogeneous: $\gamma_2=0$", fontsize=8.5,
        color=fs.MUTE, ha="left", va="center",
    )
    fs.despine(axB)
    fs.panel_title(axB, "b", "Foam-scale shape")

    fs.save(fig, "fig_systematics")
    return 0


if __name__ == "__main__":
    main()
