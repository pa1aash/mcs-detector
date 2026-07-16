# CAMPAIGN_QC.md (S6 -- 500 MeV stage)

Per-config QC. **thin** = n_events < nevt-min OR deconvolved Delta_kappa4 bootstrap CI / |Delta_kappa4| >= 30% (unresolved). Flagged configs are NOT used in the collapse fit.

| tag | N | n_events | sigma_core[mrad] | gamma2 | dk4_geom | dk4 CI frac | stuck | thin? |
|--|--:|--:|--:|--:|--:|--:|--:|:--:|
| camp_rectilinear_f20_E500 | 1994608 | 1000000 | 0.65 | 18.055 | 1.511e-11 | 0.15 | 0 |  |
| camp_rectilinear_f30_E500 | 1991830 | 1000000 | 0.99 | 11.094 | 2.114e-11 | 0.13 | 0 |  |
| camp_rectilinear_f40_E500 | 1989612 | 1000000 | 1.26 | 8.028 | 2.671e-11 | 0.12 | 0 |  |
| camp_rectilinear_f50_E500 | 1986766 | 1000000 | 1.52 | 5.989 | 3.234e-11 | 0.11 | 0 |  |
| camp_gyroid_f20_E500 | 1994598 | 1000000 | 0.86 | 14.292 | 5.056e-12 | 0.44 | 0 | YES |
| camp_gyroid_f30_E500 | 1991924 | 1000000 | 1.14 | 9.446 | 9.479e-12 | 0.29 | 0 |  |
| camp_gyroid_f40_E500 | 1989440 | 1000000 | 1.38 | 6.999 | 1.449e-11 | 0.20 | 0 |  |
| camp_gyroid_f50_E500 | 1986660 | 1000000 | 1.61 | 5.415 | 2.078e-11 | 0.18 | 0 |  |
| camp_schwarzp_f20_E500 | 1994688 | 1000000 | 0.55 | 16.783 | 1.326e-11 | 0.16 | 0 |  |
| camp_schwarzp_f30_E500 | 1992104 | 1000000 | 0.97 | 10.763 | 1.966e-11 | 0.15 | 0 |  |
| camp_schwarzp_f40_E500 | 1989656 | 1000000 | 1.27 | 7.650 | 2.557e-11 | 0.11 | 0 |  |
| camp_schwarzp_f50_E500 | 1986640 | 1000000 | 1.53 | 5.948 | 3.319e-11 | 0.10 | 0 |  |
| camp_voronoi_f30_E500 | 5973834 | 3000000 | 1.31 | 7.818 | 2.314e-12 | 0.68 | 0 | YES |

**2/13 configs flagged thin/unstable.**
