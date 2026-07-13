# R src sync inspection

- R repo: OVVO-Financial/NNS
- R commit: `4a5c74fe0ea1546fcbe25b359055d87ff69ba5e5`
- R version: `13.1`
- Previous src hash: `03ac013e5cdaa606300eea2b5c0ac16aa23a61362b555d5522e42eccd069d581`
- Incoming src hash: `61c94af1ea9fc432b7cf92b7124f57e2eb9f9ce2370834ed851ca93aad7eda7d`
- Changed files:
  - `NNS_13.1.tar.gz`
  - `NNS_13.1.zip`
  - `R/Boost.R`
  - `R/Stack.R`
  - `man/NNS.boost.Rd`
  - `man/NNS.stack.Rd`
  - `src/NNS.dll`
  - `tests/testthat/Rplots.pdf`

## Result

`src/**` changed. Review and port the affected R C++ files into NNS-core.

This inspection does not auto-port arbitrary C++ changes. After the core source has actually been updated and C++ build and tests pass, update `sync/nns_r_source.json` with the new `r_commit`, `r_version`, `r_src_tree_hash`, and `core_commit`.
