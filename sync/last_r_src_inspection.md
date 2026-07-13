# R src sync inspection

- R repo: OVVO-Financial/NNS
- R commit: `c13fb4f8162ba97559f87de03db0a3de34317ec6`
- R version: `13.1`
- Previous src hash: `03ac013e5cdaa606300eea2b5c0ac16aa23a61362b555d5522e42eccd069d581`
- Incoming src hash: `841dc0d5f74b6ffd6542b6515d43369121a7349a0ff2eed60ba006b221ede3ad`
- Changed files:
  - `DESCRIPTION`
  - `NAMESPACE`
  - `NNS_13.1.tar.gz`
  - `NNS_13.1.zip`
  - `R/Multivariate_Regression.R`
  - `R/Partition_Map.R`
  - `R/RcppExports.R`
  - `R/Regression.R`
  - `R/Stack.R`
  - `R/gvload.R`
  - `doc/NNSvignette_01_Overview.R`
  - `doc/NNSvignette_01_Overview.Rmd`
  - `doc/NNSvignette_01_Overview.html`
  - `doc/NNSvignette_03_Correlation_and_Dependence.html`
  - `doc/NNSvignette_07_Clustering_and_Regression.R`
  - `doc/NNSvignette_07_Clustering_and_Regression.Rmd`
  - `doc/NNSvignette_07_Clustering_and_Regression.html`
  - `doc/NNSvignette_08_Classification.R`
  - `doc/NNSvignette_08_Classification.Rmd`
  - `doc/NNSvignette_08_Classification.html`
  - `doc/NNSvignette_09_Forecasting.html`
  - `man/NNS.boost.Rd`
  - `man/NNS.part.Rd`
  - `man/NNS.stack.Rd`
  - `src/NNS.dll`
  - `src/RcppExports.cpp`
  - `tests/testthat/Rplots.pdf`
  - `vignettes/NNSvignette_01_Overview.R`
  - `vignettes/NNSvignette_01_Overview.Rmd`
  - `vignettes/NNSvignette_01_Overview.html`
  - `vignettes/NNSvignette_03_Correlation_and_Dependence.html`
  - `vignettes/NNSvignette_07_Clustering_and_Regression.R`
  - `vignettes/NNSvignette_07_Clustering_and_Regression.Rmd`
  - `vignettes/NNSvignette_07_Clustering_and_Regression.html`
  - `vignettes/NNSvignette_08_Classification.R`
  - `vignettes/NNSvignette_08_Classification.Rmd`
  - `vignettes/NNSvignette_08_Classification.html`
  - `vignettes/NNSvignette_09_Forecasting.html`

## Result

`src/**` changed. Review and port the affected R C++ files into NNS-core.

This inspection does not auto-port arbitrary C++ changes. After the core source has actually been updated and C++ build and tests pass, update `sync/nns_r_source.json` with the new `r_commit`, `r_version`, `r_src_tree_hash`, and `core_commit`.
