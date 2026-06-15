# Live Rcpp fidelity report

- R version: `13.0`
- R commit: `3c03a9c9d183fee06e2e8b718d2a0f4d7f45c153`
- tolerance: `1e-09`
- cases: 28  |  failures: 0

| case | status | C++ | R | abs diff |
| --- | --- | --- | --- | --- |
| lpm_d0_t3 | ok | 0.6 | 0.6 | 0 |
| upm_d0_t3 | ok | 0.4 | 0.4 | 0 |
| lpm_d0_t3 | ok | 0.6 | 0.6 | 0 |
| upm_d0_t3 | ok | 0.4 | 0.4 | 0 |
| lpm_d1_t3 | ok | 0.6 | 0.6 | 0 |
| upm_d1_t3 | ok | 0.6 | 0.6 | 0 |
| lpm_d1_t3 | ok | 0.6 | 0.6 | 0 |
| upm_d1_t3 | ok | 0.6 | 0.6 | 0 |
| lpm_d2_t3 | ok | 1 | 1 | 0 |
| upm_d2_t3 | ok | 1 | 1 | 0 |
| lpm_d2_t3 | ok | 1 | 1 | 0 |
| upm_d2_t3 | ok | 1 | 1 | 0 |
| lpm_x2_d1 | ok | 2.05102040816 | 2.05102040816 | 0 |
| upm_x2_d1 | ok | 2.05102040816 | 2.05102040816 | 0 |
| lpm_x2_d2 | ok | 10.263483965 | 10.263483965 | 0 |
| upm_x2_d2 | ok | 10.6727405248 | 10.6727405248 | 0 |
| co_lpm_d0 | ok | 0.3 | 0.3 | 0 |
| co_upm_d0 | ok | 0.275 | 0.275 | 0 |
| co_lpm_d1 | ok | 0.104773701343 | 0.104773701343 | 0 |
| co_upm_d1 | ok | 0.191756606133 | 0.191756606133 | 0 |
| d_lpm_l0_u0 | ok | 0.25 | 0.25 | 0 |
| d_upm_l0_u0 | ok | 0.175 | 0.175 | 0 |
| d_lpm_l0_u1 | ok | 0.18053466875 | 0.18053466875 | 0 |
| d_upm_l0_u1 | ok | 0.085356753125 | 0.085356753125 | 0 |
| d_lpm_l1_u0 | ok | 0.16293538125 | 0.16293538125 | 0 |
| d_upm_l1_u0 | ok | 0.191029031875 | 0.191029031875 | 0 |
| d_lpm_l1_u1 | ok | 0.152298706835 | 0.152298706835 | 0 |
| d_upm_l1_u1 | ok | 0.100336892208 | 0.100336892208 | 0 |

## Result

All fidelity cases passed. Core matches live R within tolerance.
