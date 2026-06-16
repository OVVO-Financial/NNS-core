# R src conservative sync report

- R repo: `OVVO-Financial/NNS`
- R commit: `cf18e19a8ab2584b5369857fdcda59ef3a864ec5`
- mode: dispatch

## Auto-updated provenance (fidelity-covered)
- none

## Manual review required
- `src/NNS_part.cpp`: mapped but auto=false (not covered by the fidelity suite); requires a human port and a fidelity case

## Unmapped changed R src (no port invented)
- none

## Out of scope (outside source_root)
- `R/Partition_Map.R`
- `tests/testthat/test_NNS_reg_infinite.R`

## Result

Manual attention is required. This sync must NOT be marked ready for merge until the listed files are ported by a human and the live Rcpp fidelity suite passes.
