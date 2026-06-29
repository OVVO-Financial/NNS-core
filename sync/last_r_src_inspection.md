# R src sync inspection

- R repo: OVVO-Financial/NNS
- R commit: `459164e2bb14557e44f0f2b4b448432c317afbbf`
- R version: `13.0`
- Previous src hash: `03ac013e5cdaa606300eea2b5c0ac16aa23a61362b555d5522e42eccd069d581`
- Incoming src hash: `189f103c9c054e13912e2e8ef0a9e6cb60a9fe31f518db358ce9dd80f73925e9`
- Changed files:
  - `NNS_13.0.tar.gz`
  - `NNS_13.0.zip`
  - `R/NNS_meboot.R`
  - `src/NNS.dll`

## Result

`src/**` changed. Review and port the affected R C++ files into NNS-core.

This inspection does not auto-port arbitrary C++ changes. After the core source has actually been updated and C++ build and tests pass, update `sync/nns_r_source.json` with the new `r_commit`, `r_version`, `r_src_tree_hash`, and `core_commit`.
