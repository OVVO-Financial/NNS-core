# R src sync inspection

- R repo: OVVO-Financial/NNS
- R commit: `01b4325735f14352c8e955a3a27165a4e0ce377b`
- R version: `13.1`
- Previous src hash: `03ac013e5cdaa606300eea2b5c0ac16aa23a61362b555d5522e42eccd069d581`
- Incoming src hash: `b0552b0439d1bf4ac07f8f52ffd790144180d2b1d4dbadbe1aa393a8dd3ad362`
- Changed files:
  - `R/Multivariate_Regression.R`
  - `R/RcppExports.R`
  - `R/Stack.R`
  - `src/NNS_mreg_predict.cpp`
  - `src/RcppExports.cpp`
  - `tests/testthat/test-regression-audit-repairs.R`

## Result

`src/**` changed. Review and port the affected R C++ files into NNS-core.

This inspection does not auto-port arbitrary C++ changes. After the core source has actually been updated and C++ build and tests pass, update `sync/nns_r_source.json` with the new `r_commit`, `r_version`, `r_src_tree_hash`, and `core_commit`.
