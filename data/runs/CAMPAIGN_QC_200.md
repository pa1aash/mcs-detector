# CAMPAIGN_QC.md (S6 -- 200 MeV stage)

Per-config QC. **thin** = n_events < nevt-min OR deconvolved Delta_kappa4 bootstrap CI / |Delta_kappa4| >= 30% (unresolved). Flagged configs are NOT used in the collapse fit.

| tag | N | n_events | sigma_core[mrad] | gamma2 | dk4_geom | dk4 CI frac | stuck | thin? |
|--|--:|--:|--:|--:|--:|--:|--:|:--:|
| camp_rectilinear_f20_E200 | 5985588 | 3000000 | 1.51 | 17.879 | 3.039e-10 | 0.12 | 0 |  |
| camp_rectilinear_f30_E200 | 5978026 | 3000000 | 2.30 | 11.197 | 4.431e-10 | 0.11 | 0 |  |
| camp_rectilinear_f40_E200 | 5971020 | 3000000 | 2.91 | 8.070 | 5.170e-10 | 0.10 | 0 |  |
| camp_rectilinear_f50_E200 | 5963270 | 3000000 | 3.52 | 5.939 | 5.057e-10 | 0.12 | 0 |  |
| camp_gyroid_f20_E200 | 11969882 | 6000000 | 1.99 | 14.344 | 1.238e-10 | 0.22 | 0 |  |
| camp_gyroid_f30_E200 | 11955978 | 6000000 | 2.66 | 9.605 | 2.220e-10 | 0.15 | 0 |  |
| camp_gyroid_f40_E200 | 11942216 | 6000000 | 3.20 | 7.088 | 2.779e-10 | 0.15 | 0 |  |
| camp_gyroid_f50_E200 | 11926998 | 6000000 | 3.73 | 5.416 | 3.056e-10 | 0.14 | 0 |  |
| camp_schwarzp_f20_E200 | 5985372 | 3000000 | 1.26 | 16.897 | 2.955e-10 | 0.13 | 0 |  |
| camp_schwarzp_f30_E200 | 5978286 | 3000000 | 2.25 | 10.956 | 4.134e-10 | 0.10 | 0 |  |
| camp_schwarzp_f40_E200 | 5970578 | 3000000 | 2.96 | 7.759 | 4.858e-10 | 0.11 | 0 |  |
| camp_schwarzp_f50_E200 | 5963714 | 3000000 | 3.54 | 5.916 | 5.417e-10 | 0.11 | 0 |  |
| camp_diamond_f20_E200 | 1995044 | 1000000 | 2.20 | 13.173 | 2.439e-11 | 2.34 | 0 | YES |
| camp_diamond_f30_E200 | 1992690 | 1000000 | 2.89 | 8.424 | 5.785e-12 | 13.80 | 0 | YES |
| camp_diamond_f40_E200 | 1990190 | 1000000 | 3.45 | 6.131 | -4.359e-12 | 18.77 | 0 | YES |
| camp_diamond_f50_E200 | 1987522 | 1000000 | 3.95 | 4.668 | -3.252e-11 | 3.02 | 0 | YES |
| camp_voronoi_f20_E200 | 59842934 | 30000000 | 2.37 | 12.218 | 2.413e-11 | 0.43 | 0 | YES |
| camp_voronoi_f30_E200 | 59761446 | 30000000 | 3.03 | 7.836 | 3.687e-11 | 0.40 | 0 | YES |
| camp_voronoi_f40_E200 | 59698628 | 30000000 | 3.42 | 6.277 | 1.187e-10 | 0.13 | 0 |  |
| camp_voronoi_f50_E200 | 59633628 | 30000000 | 3.83 | 5.016 | 1.087e-10 | 0.15 | 0 |  |

**6/20 configs flagged thin/unstable.**
