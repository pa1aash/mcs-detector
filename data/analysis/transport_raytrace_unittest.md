# Transport ray-tracer solid-limit unit-test (Phase B2)

SDE in pure solid (chi==1): dtheta~N(0,a dz), dy=theta dz. Analytic targets:
Var(theta(L))=aL=theta0^2, Var(y(L))=aL^3/3, Cov(theta,y)=aL^2/2.

| p [MeV] | a [rad^2/mm] | quantity | simulated | analytic | rel.err |
|--:|--:|--|--:|--:|--:|
| 100 | 1.720e-05 | Var(theta) | 1.7175e-04 | 1.7197e-04 | -0.13% |
| 100 | 1.720e-05 | Var(y) | 5.7229e-03 | 5.7323e-03 | -0.16% |
| 100 | 1.720e-05 | Cov(theta,y) | 8.5877e-04 | 8.5985e-04 | -0.13% |
| 100 | 1.720e-05 | y_rms [mm] | 7.5650e-02 | 7.5712e-02 | -0.08% |
| 200 | 3.915e-06 | Var(theta) | 3.9095e-05 | 3.9146e-05 | -0.13% |
| 200 | 3.915e-06 | Var(y) | 1.3027e-03 | 1.3049e-03 | -0.16% |
| 200 | 3.915e-06 | Cov(theta,y) | 1.9548e-04 | 1.9573e-04 | -0.13% |
| 200 | 3.915e-06 | y_rms [mm] | 3.6093e-02 | 3.6123e-02 | -0.08% |
| 500 | 7.206e-07 | Var(theta) | 7.1964e-06 | 7.2058e-06 | -0.13% |
| 500 | 7.206e-07 | Var(y) | 2.3980e-04 | 2.4019e-04 | -0.16% |
| 500 | 7.206e-07 | Cov(theta,y) | 3.5984e-05 | 3.6029e-05 | -0.13% |
| 500 | 7.206e-07 | y_rms [mm] | 1.5485e-02 | 1.5498e-02 | -0.08% |
| 1000 | 2.194e-07 | Var(theta) | 2.1914e-06 | 2.1943e-06 | -0.13% |
| 1000 | 2.194e-07 | Var(y) | 7.3023e-05 | 7.3143e-05 | -0.16% |
| 1000 | 2.194e-07 | Cov(theta,y) | 1.0958e-05 | 1.0971e-05 | -0.13% |
| 1000 | 2.194e-07 | y_rms [mm] | 8.5454e-03 | 8.5524e-03 | -0.08% |

**Solid-limit covariance unit-test: PASS** (tol 2%). The trajectory model reproduces the analytic correlated-MCS moments -> the SDE is correct.
