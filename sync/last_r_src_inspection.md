# R src sync inspection

- R repo: OVVO-Financial/NNS
- R commit: `5452831b23c77e2e0d0730623c57597809e07c1f`
- R version: `13.1`
- Previous src hash: `03ac013e5cdaa606300eea2b5c0ac16aa23a61362b555d5522e42eccd069d581`
- Incoming src hash: `a1f0b090eaa624f1ee4ad8c6c479a1fb96bcfc1942b49c8c89671743ae1f467d`
- Changed files:
  - `R/Boost.R`
  - `R/Partition_Map.R`
  - `R/Stack.R`
  - `src/NNS_distance.cpp`

## Result

`src/**` changed. Review and port the affected R C++ files into NNS-core.

This inspection does not auto-port arbitrary C++ changes. After the core source has actually been updated and C++ build and tests pass, update `sync/nns_r_source.json` with the new `r_commit`, `r_version`, `r_src_tree_hash`, and `core_commit`.
