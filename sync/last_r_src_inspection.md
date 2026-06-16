# R src sync inspection

- R repo: OVVO-Financial/NNS
- R commit: `7a3e54b144a0968d1b4ed61b6bf84be7a8d7f935`
- R version: `13.0`
- Previous src hash: `03ac013e5cdaa606300eea2b5c0ac16aa23a61362b555d5522e42eccd069d581`
- Incoming src hash: `40e5c89778a7b9dd8f553b60bda54eb95201476d0dafbb92f92549b3c627898f`
- Changed files:
  - `DESCRIPTION`
  - `NNS_13.0.tar.gz`
  - `NNS_13.0.zip`
  - `man/NNS.part.Rd`
  - `src/NNS.dll`
  - `tests/testthat/Rplots.pdf`

## Result

`src/**` changed. Review and port the affected R C++ files into NNS-core.

This inspection does not auto-port arbitrary C++ changes. After the core source has actually been updated and C++ build and tests pass, update `sync/nns_r_source.json` with the new `r_commit`, `r_version`, `r_src_tree_hash`, and `core_commit`.
