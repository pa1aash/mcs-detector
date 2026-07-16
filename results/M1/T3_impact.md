# T3 -- M1 reconstruction impact (struct vs matched slab)

Massless telescope, pencil beam, 5e6 protons/config. sigma-var ratio ~1 confirms Result 1 on the reconstructed kinks; tail ratios (bootstrap 68% CI) and track-loss ratio quantify the geometry-induced cost. Numbers written by run_m1_analysis.py.

| target | E [MeV] | sig-ratio | tail>3s | tail>4s | track-loss ratio | pull>3 ratio |
|--|--:|--:|--:|--:|--:|--:|
| rectilinear | 200 | 0.998 | 1.68x [1.67,1.69] | 1.39x | 2.35x (0.048 vs 0.020) | 1.02x |
| gyroid | 200 | 1.001 | 1.33x [1.32,1.33] | 1.15x | 2.10x (0.043 vs 0.020) | 0.97x |
| rectilinear | 1000 | 0.996 | 1.60x [1.60,1.61] | 1.33x | 1.42x (0.029 vs 0.021) | 0.84x |
| gyroid | 1000 | 1.000 | 1.29x [1.29,1.30] | 1.14x | 1.24x (0.026 vs 0.021) | 0.85x |
