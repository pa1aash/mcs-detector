# CAMPAIGN_QC.md (S6 -- 500 MeV stage)

Per-config QC. **thin** = n_events < nevt-min OR deconvolved Delta_kappa4 bootstrap CI / |Delta_kappa4| >= 30% (unresolved). Flagged configs are NOT used in the collapse fit.

| tag | N | n_events | sigma_core[mrad] | gamma2 | dk4_geom | dk4 CI frac | stuck | thin? |
|--|--:|--:|--:|--:|--:|--:|--:|:--:|
| camp_rectilinear_f20_E500 | 5984320 | 3000000 | 0.65 | 17.855 | 1.472e-11 | 0.09 | 0 |  |
| camp_rectilinear_f30_E500 | 5975762 | 3000000 | 0.99 | 11.034 | 2.097e-11 | 0.07 | 0 |  |
| camp_rectilinear_f40_E500 | 5968362 | 3000000 | 1.26 | 8.027 | 2.651e-11 | 0.07 | 0 |  |
| camp_rectilinear_f50_E500 | 5960168 | 3000000 | 1.52 | 6.006 | 3.219e-11 | 0.06 | 0 |  |
| camp_gyroid_f20_E500 | 11967466 | 6000000 | 0.86 | 14.381 | 5.296e-12 | 0.17 | 0 |  |
| camp_gyroid_f30_E500 | 11952034 | 6000000 | 1.15 | 9.527 | 1.007e-11 | 0.11 | 0 |  |
| camp_gyroid_f40_E500 | 11936598 | 6000000 | 1.38 | 7.047 | 1.484e-11 | 0.09 | 0 |  |
| camp_gyroid_f50_E500 | 11920346 | 6000000 | 1.61 | 5.362 | 1.977e-11 | 0.06 | 0 |  |
| camp_schwarzp_f20_E500 | 5983830 | 3000000 | 0.54 | 16.742 | 1.315e-11 | 0.10 | 0 |  |
| camp_schwarzp_f30_E500 | 5976332 | 3000000 | 0.97 | 10.873 | 2.019e-11 | 0.09 | 0 |  |
| camp_schwarzp_f40_E500 | 5968162 | 3000000 | 1.27 | 7.745 | 2.640e-11 | 0.07 | 0 |  |
| camp_schwarzp_f50_E500 | 5960314 | 3000000 | 1.53 | 5.913 | 3.277e-11 | 0.06 | 0 |  |
| camp_diamond_f20_E500 | 1994514 | 1000000 | 0.95 | 13.304 | 1.831e-12 | 1.20 | 0 | YES |
| camp_diamond_f30_E500 | 1991990 | 1000000 | 1.25 | 8.541 | 1.541e-12 | 1.65 | 0 | YES |
| camp_diamond_f40_E500 | 1989278 | 1000000 | 1.49 | 6.089 | 7.635e-13 | 3.82 | 0 | YES |
| camp_diamond_f50_E500 | 1986900 | 1000000 | 1.70 | 4.685 | 5.819e-13 | 5.74 | 0 | YES |
| camp_voronoi_f20_E500 | 59829442 | 30000000 | 1.02 | 12.091 | 1.048e-12 | 0.40 | 0 | YES |
| camp_voronoi_f30_E500 | 59740722 | 30000000 | 1.31 | 7.751 | 1.906e-12 | 0.27 | 0 |  |
| camp_voronoi_f40_E500 | 59672158 | 30000000 | 1.47 | 6.177 | 7.415e-12 | 0.07 | 0 |  |
| camp_voronoi_f50_E500 | 59603578 | 30000000 | 1.65 | 4.974 | 1.042e-11 | 0.06 | 0 |  |

**5/20 configs flagged thin/unstable.**
