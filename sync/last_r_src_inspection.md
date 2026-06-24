# R src sync inspection

- R repo: OVVO-Financial/NNS
- R commit: `8183e964c941d9981e19b23dfbc9fc8336903d89`
- R version: `13.0`
- Previous src hash: `03ac013e5cdaa606300eea2b5c0ac16aa23a61362b555d5522e42eccd069d581`
- Incoming src hash: `648e2e70b146fdc37d8fd694c618c2cfc328e8fc94a4947e9b84727e1b8a959a`
- Changed files:
  - `DESCRIPTION`
  - `NNS_13.0.tar.gz`
  - `NNS_13.0.zip`
  - `R/Multivariate_Regression.R`
  - `src/NNS.dll`
  - `tests/testthat/Rplots.pdf`

## Result

`src/**` changed. Review and port the affected R C++ files into NNS-core.

This inspection does not auto-port arbitrary C++ changes. After the core source has actually been updated and C++ build and tests pass, update `sync/nns_r_source.json` with the new `r_commit`, `r_version`, `r_src_tree_hash`, and `core_commit`.
