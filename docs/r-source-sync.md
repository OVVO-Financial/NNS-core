# R source sync contract

`OVVO-Financial/NNS-core` is the public portable C++ extraction layer for the R
`OVVO-Financial/NNS` package.

The R repo is authoritative. This repository records the exact R commit and
`src/**` tree hash used to produce the current C++ core.

Flow:

```text
OVVO-Financial/NNS
  -> OVVO-Financial/NNS-core
  -> OVVO-Financial/NNS-python
```

## Responsibilities

`NNS-core` is responsible for:

* receiving `nns-r-src-updated` dispatch events from `OVVO-Financial/NNS`
* fetching the requested R commit
* detecting whether `src/**` changed relative to `sync/nns_r_source.json`
* updating the portable C++ core if needed
* preserving numerical semantics from R
* running C++ build and tests
* opening a reviewable sync PR
* notifying `NNS-python` after a merged core update

`NNS-core` is not responsible for:

* Python wrapper behavior
* Python parity cache regeneration
* R package runtime changes outside the C++ source layer

It is the R repository (`OVVO-Financial/NNS`) that remains authoritative, and it
is the Python repository (`OVVO-Financial/NNS-python`) that decides Python API
parity and regenerates the Python parity cache. `NNS-core` does not make those
decisions.

## Required metadata

Every accepted sync must update `sync/nns_r_source.json` with:

* `r_repo`
* `r_commit`
* `r_version`
* `r_src_tree_hash`
* `core_repo`
* `core_commit`
* `source_paths`

The `r_src_tree_hash` is a SHA-256 over the R `src/**` tree computed with
`scripts/hash_r_src.py`. It is stable across machines because it hashes sorted
relative file paths plus file bytes, never timestamps.

## Downstream notification

After a core sync merges to `main`, this repository dispatches
`nns-core-updated` to `OVVO-Financial/NNS-python`.

The Python repository then updates its vendored `extern/NNS-core`, builds the
native backend, and runs parity gates.
