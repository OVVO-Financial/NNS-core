# R src conservative sync report

- R repo: `OVVO-Financial/NNS`
- R commit: `3c03a9c9d183fee06e2e8b718d2a0f4d7f45c153`
- mode: dispatch

## Auto-updated provenance (fidelity-covered)
- none

## Manual review required
- `src/internal_functions.cpp`: mapped but auto=false (not covered by the fidelity suite); requires a human port and a fidelity case

## Unmapped changed R src (no port invented)
- none

## Ignored Rcpp glue / build files
- `src/RcppExports.cpp`

## Out of scope (outside source_root)
- `R/ARMA.R`
- `R/RcppExports.R`

## Result

Manual attention is required. This sync must NOT be marked ready for merge until the listed files are ported by a human and the live Rcpp fidelity suite passes.
