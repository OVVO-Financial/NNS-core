# R src sync inspection

- R repo: OVVO-Financial/NNS
- R commit: `cc7b5f971db01170f99660c822d0f35ec683375d`
- R version: `13.0`
- Previous src hash: `03ac013e5cdaa606300eea2b5c0ac16aa23a61362b555d5522e42eccd069d581`
- Incoming src hash: `6e83ca62509cfff844921107cd7ff03b3789205913357465126c72508d8fb05b`
- Changed files:
  - `DESCRIPTION`
  - `NAMESPACE`
  - `R/Copula.R`
  - `R/Multivariate_Regression.R`
  - `R/Partial_Moments.R`
  - `R/gvload.R`
  - `doc/NNSvignette_01_Overview.html`
  - `src/NNS.dll`
  - `vignettes/NNSvignette_01_Overview.html`

## Result

`src/**` changed. Review and port the affected R C++ files into NNS-core.

This inspection does not auto-port arbitrary C++ changes. After the core source has actually been updated and C++ build and tests pass, update `sync/nns_r_source.json` with the new `r_commit`, `r_version`, `r_src_tree_hash`, and `core_commit`.
