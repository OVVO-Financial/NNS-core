# R src conservative sync report

- R repo: `OVVO-Financial/NNS`
- R commit: `7dec0a0e39ec9e604b9e9c47c7ecf2d99673abf0`
- mode: dispatch

## Auto-updated provenance (fidelity-covered)
- none

## Manual review required
- `src/NNS_dep.cpp`: mapped but auto=false (not covered by the fidelity suite); requires a human port and a fidelity case

## Unmapped changed R src (no port invented)
- none

## Out of scope (outside source_root)
- `R/Copula.R`
- `tests/testthat/test_Copula.R`

## Result

Manual attention is required. This sync must NOT be marked ready for merge until the listed files are ported by a human and the live Rcpp fidelity suite passes.
