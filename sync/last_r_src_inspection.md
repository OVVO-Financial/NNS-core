# R src sync inspection

- R repo: OVVO-Financial/NNS
- R commit: `34d153938e9ba473e3e6c687687ef2d56f54e66c`
- R version: `13.0`
- Previous src hash: `03ac013e5cdaa606300eea2b5c0ac16aa23a61362b555d5522e42eccd069d581`
- Incoming src hash: `926c56bbf6161b4bc80b363c5464ce9f0404711e3e7634847b5096b3b9214fbd`
- Changed files:
  - `R/ARMA.R`
  - `R/Dependence.R`
  - `R/RcppExports.R`
  - `R/Regression.R`
  - `src/RcppExports.cpp`
  - `src/central_tendencies.cpp`
  - `src/central_tendencies.h`
  - `src/nns_reg_points.cpp`

## Result

`src/**` changed. Review and port the affected R C++ files into NNS-core.

This inspection does not auto-port arbitrary C++ changes. After the core source has actually been updated and C++ build and tests pass, update `sync/nns_r_source.json` with the new `r_commit`, `r_version`, `r_src_tree_hash`, and `core_commit`.
