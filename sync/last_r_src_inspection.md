# R src sync inspection

- R repo: OVVO-Financial/NNS
- R commit: `7dec0a0e39ec9e604b9e9c47c7ecf2d99673abf0`
- R version: `13.1`
- Previous src hash: `03ac013e5cdaa606300eea2b5c0ac16aa23a61362b555d5522e42eccd069d581`
- Incoming src hash: `381ac35d86fc3554a4b486adf44d45979387e1d2d389e22e33e568003ddad507`
- Changed files:
  - `R/Copula.R`
  - `src/NNS_dep.cpp`
  - `tests/testthat/test_Copula.R`

## Result

`src/**` changed. Review and port the affected R C++ files into NNS-core.

This inspection does not auto-port arbitrary C++ changes. After the core source has actually been updated and C++ build and tests pass, update `sync/nns_r_source.json` with the new `r_commit`, `r_version`, `r_src_tree_hash`, and `core_commit`.
