# R src sync inspection

- R repo: OVVO-Financial/NNS
- R commit: `9012ec236283e6584323cd482b9f0141b8b32a8b`
- R version: `13.1`
- Previous src hash: `03ac013e5cdaa606300eea2b5c0ac16aa23a61362b555d5522e42eccd069d581`
- Incoming src hash: `12d7788da1e2dd6f3cc666a8703a732089aabcbc8adfd283e04911e494fa78eb`
- Changed files:
  - `DESCRIPTION`
  - `NNS_13.0.tar.gz`
  - `NNS_13.0.zip`
  - `NNS_13.1.tar.gz`
  - `NNS_13.1.zip`
  - `R/ARMA_optim.R`
  - `R/NNS_VAR.R`
  - `README.md`
  - `man/NNS.ARMA.optim.Rd`
  - `man/NNS.VAR.Rd`
  - `src/NNS.dll`

## Result

`src/**` changed. Review and port the affected R C++ files into NNS-core.

This inspection does not auto-port arbitrary C++ changes. After the core source has actually been updated and C++ build and tests pass, update `sync/nns_r_source.json` with the new `r_commit`, `r_version`, `r_src_tree_hash`, and `core_commit`.
