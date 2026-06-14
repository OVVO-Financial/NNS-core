# R src sync inspection

- R repo: OVVO-Financial/NNS
- R commit: `3250bc985883e21eeb8fa4bf97656d8469de2836`
- R version: `13.0`
- Previous src hash: `03ac013e5cdaa606300eea2b5c0ac16aa23a61362b555d5522e42eccd069d581`
- Incoming src hash: `705094f1c69478faf2cb89347730a0dc48c9340152ad42a606d74cc33c4d1e66`
- Changed files:
  - `NNS_13.0.tar.gz`
  - `NNS_13.0.zip`
  - `R/ARMA.R`
  - `src/NNS.dll`
  - `tests/testthat/Rplots.pdf`

## Result

`src/**` changed. Review and port the affected R C++ files into NNS-core.

This inspection does not auto-port arbitrary C++ changes. After the core source has actually been updated and C++ build and tests pass, update `sync/nns_r_source.json` with the new `r_commit`, `r_version`, `r_src_tree_hash`, and `core_commit`.
