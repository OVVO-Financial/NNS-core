# R src sync inspection

- R repo: OVVO-Financial/NNS
- R commit: `79f3255b1a4629edfc930d25423c6c46d676ed23`
- R version: `13.1`
- Previous src hash: `03ac013e5cdaa606300eea2b5c0ac16aa23a61362b555d5522e42eccd069d581`
- Incoming src hash: `3722fd299c3f4e530a59b75793fdfafa4699c58599521de7d3e55e93f0605bcd`
- Changed files:
  - `DESCRIPTION`
  - `NNS_13.1.tar.gz`
  - `NNS_13.1.zip`
  - `book/chapter-03-distribution-theory-from-partial-moments.Rmd`
  - `book/chapter-09-why-correlation-fails.Rmd`
  - `book/chapter-10-directional-dependence.Rmd`
  - `book/chapter-11-directional-spectral-decomposition.Rmd`
  - `book/chapter-12-copula-interpretation.Rmd`
  - `book/chapter-15-distribution-comparison.Rmd`
  - `book/chapter-17-prediction-intervals.Rmd`
  - `book/chapter-18-recursive-mean-split-estimation.Rmd`
  - `book/chapter-20-synthetic-data-and-maximum-entropy-bootstrap.Rmd`
  - `book/chapter-21-clustering.Rmd`
  - `book/chapter-22-nonparametric-regression.Rmd`
  - `book/chapter-23-classification.Rmd`
  - `book/chapter-24-ensemble-methods.Rmd`
  - `book/chapter-25-nonparametric-time-series-models.Rmd`
  - `book/chapter-26-multivariate-forecasting.Rmd`
  - `book/chapter-28-appendix-notation-and-function-reference.Rmd`
  - `book/index.Rmd`
  - `src/NNS.dll`

## Result

`src/**` changed. Review and port the affected R C++ files into NNS-core.

This inspection does not auto-port arbitrary C++ changes. After the core source has actually been updated and C++ build and tests pass, update `sync/nns_r_source.json` with the new `r_commit`, `r_version`, `r_src_tree_hash`, and `core_commit`.
