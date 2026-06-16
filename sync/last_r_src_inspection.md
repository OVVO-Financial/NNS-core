# R src sync inspection

- R repo: OVVO-Financial/NNS
- R commit: `a73b3b17847c38d15223e3fcd073075942b58106`
- R version: `13.0`
- Previous src hash: `03ac013e5cdaa606300eea2b5c0ac16aa23a61362b555d5522e42eccd069d581`
- Incoming src hash: `2ba77ee3cd0f8192bba0e6621b2b179cd5a5afb2bf243c5ec85d13ac2d10631f`
- Changed files:
  - `R/Regression.R`
  - `src/nns_reg_points.cpp`

## Result

`src/**` changed. Review and port the affected R C++ files into NNS-core.

This inspection does not auto-port arbitrary C++ changes. After the core source has actually been updated and C++ build and tests pass, update `sync/nns_r_source.json` with the new `r_commit`, `r_version`, `r_src_tree_hash`, and `core_commit`.
