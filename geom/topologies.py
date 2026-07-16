"""topologies.py -- S4 target topologies as analytic / voxel material fields.

Each topology exposes a binary indicator chi(x,y,z) on a periodic cubic cell of
edge `cell`, parametrised by a single tuning scalar (TPMS level, strut radius, or
foam wall) that is bisected to hit a target solid fraction f.

Topologies:
  rectilinear : cubic strut lattice (struts along the x,y,z cell edges). 3-D, has
                z-structure, l_int tunable via cell -- the printable "rectilinear
                lattice". (The 2-D extruded grid-infill is the analytic N_eff=1
                extruded limit, already covered by theory.py B9.)
  gyroid      : G = sin x cos y + sin y cos z + sin z cos x   (solid = G > level)
  schwarzp    : P = cos x + cos y + cos z                     (solid = P > level)
  diamond     : D = sum of the 4 diamond terms               (solid = D > level)
  voronoi     : stochastic closed-cell foam (Voronoi-face walls) -- the near-Markov
                exponential-chord case. Voxel field (not analytic).

All lengths in mm. Beam axis = +z.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

TPMS = ("gyroid", "schwarzp", "diamond")
ANALYTIC = TPMS + ("rectilinear",)
ALL_TOPOS = ANALYTIC + ("voronoi",)


# --------------------------------------------------------------------------
# implicit fields (period = cell)
# --------------------------------------------------------------------------
def _tpms_field(name, x, y, z, cell):
    k = 2.0 * np.pi / cell
    kx, ky, kz = k * x, k * y, k * z
    if name == "gyroid":
        return np.sin(kx) * np.cos(ky) + np.sin(ky) * np.cos(kz) + np.sin(kz) * np.cos(kx)
    if name == "schwarzp":
        return np.cos(kx) + np.cos(ky) + np.cos(kz)
    if name == "diamond":
        return (np.sin(kx) * np.sin(ky) * np.sin(kz)
                + np.sin(kx) * np.cos(ky) * np.cos(kz)
                + np.cos(kx) * np.sin(ky) * np.cos(kz)
                + np.cos(kx) * np.cos(ky) * np.sin(kz))
    raise ValueError(name)


def _dist_to_grid(v, cell):
    """distance from v to the nearest multiple of cell."""
    m = np.mod(v, cell)
    return np.minimum(m, cell - m)


def chi_analytic(name, x, y, z, cell, p):
    """Binary indicator for an analytic topology. p = level (TPMS) or radius (strut)."""
    if name in TPMS:
        return _tpms_field(name, x, y, z, cell) > p
    if name == "rectilinear":
        dx = _dist_to_grid(x, cell); dy = _dist_to_grid(y, cell); dz = _dist_to_grid(z, cell)
        r = p
        return ((np.hypot(dy, dz) < r) | (np.hypot(dx, dz) < r) | (np.hypot(dx, dy) < r))
    raise ValueError(name)


def _frac_analytic(name, cell, p, n=64):
    """solid fraction over one unit cell for an analytic topology."""
    g = (np.arange(n) + 0.5) / n * cell
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    return float(chi_analytic(name, X, Y, Z, cell, p).mean())


def tune_analytic(name, cell, target_f, n=64, tol=2e-3, iters=60):
    """Bisect the tuning scalar p so the unit-cell fraction == target_f."""
    if name in TPMS:
        lo, hi = -3.2, 3.2            # level: f decreasing in level
        for _ in range(iters):
            mid = 0.5 * (lo + hi)
            f = _frac_analytic(name, cell, mid, n)
            if abs(f - target_f) < tol:
                return mid, f
            if f > target_f:          # too much solid -> raise level
                lo = mid
            else:
                hi = mid
        return mid, f
    if name == "rectilinear":
        lo, hi = 0.0, cell * 0.62      # radius: f increasing in r
        for _ in range(iters):
            mid = 0.5 * (lo + hi)
            f = _frac_analytic(name, cell, mid, n)
            if abs(f - target_f) < tol:
                return mid, f
            if f < target_f:
                lo = mid
            else:
                hi = mid
        return mid, f
    raise ValueError(name)


# --------------------------------------------------------------------------
# voxel fields (voronoi foam, and a grid sampler for STL generation)
# --------------------------------------------------------------------------
def _voronoi_seeds(cell, box, rng, halo=1):
    """jittered-grid seeds at density ~1/cell^3 over box=(Lx,Ly,Lz) (+halo cells)."""
    nx, ny, nz = [int(np.ceil(b / cell)) + 2 * halo for b in box]
    gx = (np.arange(nx) - halo + 0.5) * cell
    gy = (np.arange(ny) - halo + 0.5) * cell
    gz = (np.arange(nz) - halo + 0.5) * cell
    G = np.stack(np.meshgrid(gx, gy, gz, indexing="ij"), -1).reshape(-1, 3)
    G += (rng.random(G.shape) - 0.5) * cell * 0.8   # jitter
    return G


def voronoi_delta(cell, box, grid_shape, rng):
    """(d2-d1) field on a voxel grid (distance to 2nd-nearest minus nearest seed).
    Thin Voronoi-face walls = small (d2-d1)."""
    seeds = _voronoi_seeds(cell, box, rng)
    tree = cKDTree(seeds)
    nx, ny, nz = grid_shape
    xs = (np.arange(nx) + 0.5) / nx * box[0]
    ys = (np.arange(ny) + 0.5) / ny * box[1]
    zs = (np.arange(nz) + 0.5) / nz * box[2]
    P = np.stack(np.meshgrid(xs, ys, zs, indexing="ij"), -1).reshape(-1, 3)
    d, _ = tree.query(P, k=2, workers=-1)
    return (d[:, 1] - d[:, 0]).reshape(grid_shape)


def voronoi_field(cell, box, target_f, grid_shape, rng):
    """Closed-cell foam chi tuned to target_f exactly via a quantile of (d2-d1):
    returns (chi, wall, f_achieved). No bisection -- one distance query."""
    delta = voronoi_delta(cell, box, grid_shape, rng)
    wall = float(np.quantile(delta, target_f))
    chi = delta < wall
    return chi, wall, float(chi.mean())


# --------------------------------------------------------------------------
# voxel sampler for STL generation (all topologies on a common grid)
# --------------------------------------------------------------------------
def ray_chi_analytic(name, cell, p, L, nxy=80, dz=None, n_cells_xy=1, chunk=400):
    """(n_rays, nz) indicator array for an analytic topology: rays at nxy*nxy
    transverse positions spanning n_cells_xy unit cells, sampled along z over
    [0,L]. Chunked over rays to bound memory (small cells -> large nz)."""
    if dz is None:
        dz = cell / 30.0
    nz = int(round(L / dz))
    z = (np.arange(nz) + 0.5) * dz
    span = n_cells_xy * cell
    g = (np.arange(nxy) + 0.5) / nxy * span
    XX, YY = np.meshgrid(g, g, indexing="ij")
    xr = XX.ravel(); yr = YY.ravel()
    n_rays = xr.size
    out = np.empty((n_rays, nz), dtype=np.uint8)
    for i in range(0, n_rays, chunk):
        xc = xr[i:i + chunk][:, None]
        yc = yr[i:i + chunk][:, None]
        out[i:i + chunk] = chi_analytic(name, xc, yc, z[None, :], cell, p)
    return out, dz


def sample_field(name, box, cell, p, grid_shape, rng=None):
    """Float field on a regular grid over box; marching-cubes at level 0.5 of the
    returned (0/1) array gives the surface. For analytic topologies p=level/radius;
    for voronoi p=wall."""
    nx, ny, nz = grid_shape
    xs = (np.arange(nx) + 0.5) / nx * box[0]
    ys = (np.arange(ny) + 0.5) / ny * box[1]
    zs = (np.arange(nz) + 0.5) / nz * box[2]
    if name == "voronoi":
        chi, _, _ = voronoi_field(cell, box, p, grid_shape, rng)  # p = target_f
        return chi.astype(np.float32)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    return chi_analytic(name, X, Y, Z, cell, p).astype(np.float32)
