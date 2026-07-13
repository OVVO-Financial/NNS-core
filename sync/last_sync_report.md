# R src conservative sync report

- R repo: `OVVO-Financial/NNS`
- R commit: `15d74ef047d5c60cda53dfa8ca6db13c96ad3f6d`
- mode: dispatch

## Auto-updated provenance (fidelity-covered)
- none

## Manual review required
- none

## Unmapped changed R src (no port invented)
- `src/NNS_mreg_predict.cpp`
- `src/NNS_mreg_setup.cpp`
- `src/NNS_stack_fast.cpp`

## Ignored Rcpp glue / build files
- `src/RcppExports.cpp`

## Out of scope (outside source_root)
- `R/Multivariate_Regression.R`
- `R/RcppExports.R`
- `R/Regression.R`
- `R/Stack.R`
- `benchmarks/benchmark_stack_native.R`
- `benchmarks/stack_baseline_report.md`
- `benchmarks/stack_baseline_results.csv`
- `benchmarks/stack_native_report.md`
- `benchmarks/stack_native_results.csv`
- `tests/testthat/helper-stack-parity.R`
- `tests/testthat/test-mreg-path-native.R`
- `tests/testthat/test-mreg-setup-native.R`
- `tests/testthat/test-stack-method1-native-parity.R`
- `tests/testthat/test-stack-method2-native-parity.R`
- `tests/testthat/test-stack-native-parity.R`
- `tests/testthat/test-univariate-fast-parity.R`
- `tests/testthat/test-xstar-path-native.R`

## Result

Manual attention is required. This sync must NOT be marked ready for merge until the listed files are ported by a human and the live Rcpp fidelity suite passes.
