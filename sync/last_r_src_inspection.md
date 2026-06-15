# R src sync inspection

- R repo: OVVO-Financial/NNS
- R commit: `36b760589ca0a08639de65d923bac5a82b0d2ab9`
- R version: `13.0`
- Previous src hash: `03ac013e5cdaa606300eea2b5c0ac16aa23a61362b555d5522e42eccd069d581`
- Incoming src hash: `ae5bdfa9cc458ac1f2472efb7a4d37190fb878b957604911ffe24bb8e8850179`
- Changed files:
  - `DESCRIPTION`
  - `NNS_13.0.tar.gz`
  - `NNS_13.0.zip`
  - `R/ARMA.R`
  - `src/NNS.dll`

## Result

`src/**` changed. Review and port the affected R C++ files into NNS-core.

This inspection does not auto-port arbitrary C++ changes. After the core source has actually been updated and C++ build and tests pass, update `sync/nns_r_source.json` with the new `r_commit`, `r_version`, `r_src_tree_hash`, and `core_commit`.
