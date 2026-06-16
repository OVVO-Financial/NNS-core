# R src conservative sync report

- R repo: `OVVO-Financial/NNS`
- R commit: `34d153938e9ba473e3e6c687687ef2d56f54e66c`
- mode: dispatch

## Auto-updated provenance (fidelity-covered)
- none

## Manual review required
- `src/central_tendencies.cpp`: mapped but auto=false (not covered by the fidelity suite); requires a human port and a fidelity case
- `src/central_tendencies.h`: mapped but auto=false (not covered by the fidelity suite); requires a human port and a fidelity case

## Unmapped changed R src (no port invented)
- `src/nns_reg_points.cpp`

## Ignored Rcpp glue / build files
- `src/RcppExports.cpp`

## Out of scope (outside source_root)
- `R/ARMA.R`
- `R/Dependence.R`
- `R/RcppExports.R`
- `R/Regression.R`

## Result

Manual attention is required. This sync must NOT be marked ready for merge until the listed files are ported by a human and the live Rcpp fidelity suite passes.
