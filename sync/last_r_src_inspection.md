# R src sync inspection

- R repo: OVVO-Financial/NNS
- R commit: `15d74ef047d5c60cda53dfa8ca6db13c96ad3f6d`
- R version: `13.1`
- Previous src hash: `03ac013e5cdaa606300eea2b5c0ac16aa23a61362b555d5522e42eccd069d581`
- Incoming src hash: `cbf8abd8426c7a5d9de3fe63959fa8adb453790521c0742e832b693d2a962457`
- Changed files:
  - `R/Multivariate_Regression.R`
  - `R/RcppExports.R`
  - `R/Regression.R`
  - `R/Stack.R`
  - `benchmarks/benchmark_stack_native.R`
  - `benchmarks/stack_baseline_report.md`
  - `benchmarks/stack_baseline_results.csv`
  - `benchmarks/stack_native_report.md`
  - `benchmarks/stack_native_results.csv`
  - `src/NNS_mreg_predict.cpp`
  - `src/NNS_mreg_setup.cpp`
  - `src/NNS_stack_fast.cpp`
  - `src/RcppExports.cpp`
  - `tests/testthat/helper-stack-parity.R`
  - `tests/testthat/test-mreg-path-native.R`
  - `tests/testthat/test-mreg-setup-native.R`
  - `tests/testthat/test-stack-method1-native-parity.R`
  - `tests/testthat/test-stack-method2-native-parity.R`
  - `tests/testthat/test-stack-native-parity.R`
  - `tests/testthat/test-univariate-fast-parity.R`
  - `tests/testthat/test-xstar-path-native.R`

## Result

`src/**` changed. Review and port the affected R C++ files into NNS-core.

This inspection does not auto-port arbitrary C++ changes. After the core source has actually been updated and C++ build and tests pass, update `sync/nns_r_source.json` with the new `r_commit`, `r_version`, `r_src_tree_hash`, and `core_commit`.
