# R src sync inspection

- R repo: OVVO-Financial/NNS
- R commit: `cf18e19a8ab2584b5369857fdcda59ef3a864ec5`
- R version: `13.0`
- Previous src hash: `03ac013e5cdaa606300eea2b5c0ac16aa23a61362b555d5522e42eccd069d581`
- Incoming src hash: `0d6458c3893a2936a368c9e130de68fccf81f7b31d7369a98f8b8355547fc731`
- Changed files:
  - `R/Partition_Map.R`
  - `src/NNS_part.cpp`
  - `tests/testthat/test_NNS_reg_infinite.R`

## Result

`src/**` changed. Review and port the affected R C++ files into NNS-core.

This inspection does not auto-port arbitrary C++ changes. After the core source has actually been updated and C++ build and tests pass, update `sync/nns_r_source.json` with the new `r_commit`, `r_version`, `r_src_tree_hash`, and `core_commit`.
