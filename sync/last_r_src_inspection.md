# R src sync inspection

- R repo: OVVO-Financial/NNS
- R commit: `3c03a9c9d183fee06e2e8b718d2a0f4d7f45c153`
- R version: `13.0`
- Previous src hash: `03ac013e5cdaa606300eea2b5c0ac16aa23a61362b555d5522e42eccd069d581`
- Incoming src hash: `1b3ba783a489734c69ebcc2a767486187cb51fb869cf2bce91be41e829c8d9ec`
- Changed files:
  - `R/ARMA.R`
  - `R/RcppExports.R`
  - `src/RcppExports.cpp`
  - `src/internal_functions.cpp`

## Result

`src/**` changed. Review and port the affected R C++ files into NNS-core.

This inspection does not auto-port arbitrary C++ changes. After the core source has actually been updated and C++ build and tests pass, update `sync/nns_r_source.json` with the new `r_commit`, `r_version`, `r_src_tree_hash`, and `core_commit`.
