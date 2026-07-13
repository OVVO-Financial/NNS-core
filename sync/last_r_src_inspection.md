# R src sync inspection

- R repo: OVVO-Financial/NNS
- R commit: `022be470d765600e09c13a7e53cd0de03d72e391`
- R version: `13.1`
- Previous src hash: `03ac013e5cdaa606300eea2b5c0ac16aa23a61362b555d5522e42eccd069d581`
- Incoming src hash: `80d42b9bf79bd0283b9e819da29db1c46aaa1db5cf16e8dd529eab744242eec6`
- Changed files:
  - `NNS_13.1.tar.gz`
  - `NNS_13.1.zip`
  - `R/Stack.R`
  - `man/NNS.boost.Rd`
  - `man/NNS.stack.Rd`
  - `src/NNS.dll`

## Result

`src/**` changed. Review and port the affected R C++ files into NNS-core.

This inspection does not auto-port arbitrary C++ changes. After the core source has actually been updated and C++ build and tests pass, update `sync/nns_r_source.json` with the new `r_commit`, `r_version`, `r_src_tree_hash`, and `core_commit`.
