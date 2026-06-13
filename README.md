# NNS-core
Core NNS functions for C++ port

## Role in the NNS project

`NNS-core` is the portable C++ core extracted from the R package
[`OVVO-Financial/NNS`](https://github.com/OVVO-Financial/NNS).

The R package is authoritative. This repository records the exact R commit and
`src/**` tree hash in `sync/nns_r_source.json`.

`OVVO-Financial/NNS` exposes two truths. `NNS-core` owns only the **`src/**`
native truth** (the C++/Rcpp numerics), verified by its live Rcpp fidelity
suite. The separate **R API behavior truth** is verified directly and is
`NNS-python`'s responsibility, not this repository's.

Downstream Python integration happens in
[`OVVO-Financial/NNS-python`](https://github.com/OVVO-Financial/NNS-python),
where the two truths converge.

See [`docs/r-source-sync.md`](docs/r-source-sync.md) for the full source sync
contract.
