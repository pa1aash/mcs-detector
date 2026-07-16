"""figstyle.py -- publication figure design system (v4, JINST).

Single source of truth for the figure set. Renders through matplotlib's PGF/LaTeX
backend so text and math are typeset in the SAME font family as the JINST body
(newtx / Times, matching jinstpub.sty) via pdflatex.

Rules encoded here (used by every figure script):
  - Figures are BUILT AT THEIR DISPLAY WIDTH so fonts render at true point size.
    JINST text width = 0.72 * A4 = 151.2 mm = 5.95 in (FULL). Figures shown at
    0.78\\linewidth use W078, at 0.62\\linewidth use W062.
  - One FIXED topology encoding in every figure:
      rectilinear = red circle, Schwarz-P = green square, gyroid = blue triangle,
      Voronoi = orange inverted triangle, diamond = grey thin diamond.
  - Energy encoded by marker fill: 200 MeV open, 500 MeV half-filled (bottom),
    1000 MeV filled.  Use estyle(E, color).
  - Physics-list variants by marker EDGE STYLE: WentzelVI solid edge, Urban dashed
    edge.  Use phys_scatter(ax, ...).
  - All text >= 8.5 pt at final size.  ``save`` enforces this invariant.
  - Multi-panel figures use a short descriptive header as well as a panel letter;
    a reader should understand each panel before opening the caption.
  - ANNOTATION DOCTRINE (professional-clarity pass): panels carry at most two short
    in-axes labels (<= 4 words each); no sentences on the canvas (captions talk);
    no leader lines that cross data; at most ONE legend per figure; direct curve
    labels preferred for <= 3 series; one shaded band per panel, subtle; no
    half-empty axes (ranges follow the data).
  - Colour is never the only curve encoding: topology curves also use fixed dash
    patterns.  Markers retain the paper-wide topology shapes.
  - save() writes the vector PDF (paper) and a 300-dpi PNG preview under
    figs/preview/ (for review).
"""
from __future__ import annotations
import json, os
import matplotlib

matplotlib.use("pgf")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
AOUT = os.path.join(ROOT, "data", "analysis")
FIGS = os.path.join(ROOT, "figs")
PREVIEW = os.path.join(FIGS, "preview")

# ---- print geometry (inches): JINST text width and the used display fractions ----
FULL = 5.95           # \linewidth  (0.72 * A4 width = 151.2 mm)
W078 = 0.78 * FULL    # figures shown at width=0.78\linewidth
W062 = 0.62 * FULL    # figures shown at width=0.62\linewidth
GOLDEN = 0.618
# legacy aliases (older scripts)
HALF = W062
SINGLE, ONEHALF, DOUBLE = W062, W078, FULL

# ---- ONE fixed colour + marker per topology (colourblind- and grayscale-safe) ----
COLORS = {
    "rectilinear": "#b23a48",   # red
    "schwarzp":    "#3f7d5b",   # green
    "gyroid":      "#2f5d8a",   # blue
    "voronoi":     "#c9821f",   # orange
    "diamond":     "#8b8b8b",   # grey (recedes)
}
MARKERS = {"rectilinear": "o", "schwarzp": "s", "gyroid": "^",
           "voronoi": "v", "diamond": "d"}
LINESTYLES = {
    "rectilinear": "-",
    "schwarzp": (0, (5.0, 2.0)),
    "gyroid": (0, (4.0, 1.6, 1.2, 1.6)),
    "voronoi": (0, (1.3, 1.5)),
    "diamond": (0, (7.0, 2.2)),
}
TOPO_LABEL = {"rectilinear": "rectilinear", "schwarzp": "Schwarz-P",
              "gyroid": "gyroid", "voronoi": "Voronoi", "diamond": "diamond"}

HERO = "#1a1a1a"        # the key (derived) line
MUTE = "#7a7a7a"        # supporting / recessed
FAINT = "#c2c2c2"
FOAM = "#efd9a7"        # detector-foam + print pore band fill
FOAM_EDGE = "#8f6d1d"
FOAM_BAND = (0.10, 0.50)  # mm
SOLID_C, VOID_C = "#3a3f4a", "#eef1f5"   # geometry slices: solid phase / void

# physics-list EDGE styles (not colours): WentzelVI solid, Urban dashed
PHYS_EDGE = {"WentzelVI": "solid", "Urban": (0, (2.0, 1.2))}

# ---- shared line weights / marker+error-bar sizes (single source of truth) ----
# Every figure pulls these; do not hardcode weights/sizes per figure.
LW_DATA = 1.7          # primary data lines and curves (the focal result heaviest)
LW_FIT = 1.15          # secondary fit / prediction lines (drawn dashed, subordinate)
LW_REF = 0.9           # reference lines (y=x, unity, zero) -- thin, grey, behind data
LW_SPINE = 0.8
MS = 5.5               # standard data-marker size
MS_SMALL = 4.6
REF = MUTE             # colour for reference lines and guides
BAND = FAINT           # colour for reference/uncertainty bands (drawn faint, behind)
DASH_FIT = (0, (5, 2.6))   # dash pattern for subordinate fit lines
MIN_FONT_PT = 8.5


def ebar(color):
    """Consistent thin error-bar kwargs in the series colour."""
    return dict(ecolor=color, elinewidth=0.9, capsize=2.2, capthick=0.8)


def estyle(E, color):
    """Marker style kwargs encoding energy by fill: 200 open / 500 half / 1000 full."""
    E = int(E)
    if E == 200:
        return dict(mfc="white", mec=color)
    if E == 500:
        return dict(mfc=color, mfcalt="white", fillstyle="bottom", mec=color)
    return dict(mfc=color, mec=color)   # 1000


def phys_scatter(ax, x, y, color, phys, marker="o", size=34, zorder=5, filled=True):
    """Scatter with the physics list encoded in the marker EDGE style."""
    return ax.scatter(x, y, s=size, marker=marker,
                      facecolor=(color if filled else "white"),
                      edgecolor=(HERO if filled else color),
                      linestyle=PHYS_EDGE[phys], linewidth=1.1, zorder=zorder)


def set_style():
    plt.rcParams.update({
        "pgf.texsystem": "pdflatex",
        "text.usetex": True,
        "font.family": "serif",
        "pgf.rcfonts": False,
        "pgf.preamble": r"\usepackage{newtxtext}\usepackage{newtxmath}"
                        r"\usepackage[T1]{fontenc}"
                        r"\usepackage{siunitx}",
        "font.size": 10,
        "axes.labelsize": 10, "axes.titlesize": 10,
        "xtick.labelsize": 9, "ytick.labelsize": 9,
        "legend.fontsize": 9, "legend.frameon": False,
        "legend.handlelength": 1.5, "legend.handletextpad": 0.5,
        "legend.labelspacing": 0.35, "legend.borderaxespad": 0.4,
        "axes.linewidth": 0.8,
        "lines.linewidth": 1.6, "lines.markersize": 5,
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "xtick.major.size": 3.2, "ytick.major.size": 3.2,
        "xtick.minor.size": 1.8, "ytick.minor.size": 1.8,
        "xtick.direction": "out", "ytick.direction": "out",
        "axes.grid": False,
        "figure.constrained_layout.use": True,
        "figure.constrained_layout.h_pad": 0.04,
        "figure.constrained_layout.w_pad": 0.04,
        "savefig.dpi": 600, "figure.dpi": 600,
    })


def despine(ax, keep=("left", "bottom")):
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(s in keep)
    ax.tick_params(top=False, right=False)


def panel(ax, letter, dx=-0.0, dy=1.02):
    """Bold panel letter just outside the top-left of the axes."""
    ax.text(dx, dy, r"\textbf{(%s)}" % letter, transform=ax.transAxes,
            ha="right", va="bottom", fontsize=10.5)


def panel_title(ax, letter, title, y=1.035):
    """Left-aligned panel letter + concise descriptive header above an axis."""
    ax.text(0.0, y, rf"\textbf{{({letter})}}\quad {title}", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=9.5, color=HERO)


def faint_grid(ax, axis="y"):
    ax.grid(True, axis=axis, color=FAINT, lw=0.5, alpha=0.7, zorder=0)
    ax.set_axisbelow(True)


def load(name):
    return json.load(open(os.path.join(AOUT, name)))


def loadr(*parts):
    """Load a JSON from results/<parts>."""
    return json.load(open(os.path.join(ROOT, "results", *parts)))


def save(fig, stem):
    # Publication figures are built at their exact display width.  A hard check is
    # preferable to discovering 6--7 pt inset labels only after manuscript assembly.
    import matplotlib.text as mtext
    undersized = []
    for artist in fig.findobj(match=mtext.Text):
        if artist.get_text().strip() and artist.get_visible() and artist.get_fontsize() < MIN_FONT_PT:
            undersized.append((artist.get_text(), artist.get_fontsize()))
    if undersized:
        details = ", ".join(f"{s!r} ({z:g} pt)" for s, z in undersized[:8])
        raise ValueError(f"{stem}: text below {MIN_FONT_PT:g} pt: {details}")
    out = os.path.join(FIGS, stem + ".pdf")
    fig.savefig(out, facecolor="white", transparent=False)
    os.makedirs(PREVIEW, exist_ok=True)
    fig.savefig(os.path.join(PREVIEW, stem + ".png"), dpi=300,
                facecolor="white", transparent=False)
    print("wrote", out)
    return out
