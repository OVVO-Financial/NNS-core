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

## Conservative sync engine

The sync workflow may automatically update `NNS-core`, but only within strict
limits encoded in `sync/r_src_map.json` and enforced by `scripts/sync_adapter.py`.

The adapter's **only** transformation to the C++ core is rewriting a managed
provenance comment block that records the incoming R commit:

```cpp
// ---- NNS-CORE PROVENANCE (managed by scripts/sync_adapter.py; do not hand-edit) ----
// NNS-R-SOURCE: src/partial_moments.cpp @ <r_commit>
// ---- end NNS-CORE PROVENANCE ----
```

It never edits C++ logic. For each changed R `src` file the adapter:

* fetches the incoming R source file from the requested commit
* verifies the expected function signatures (markers) are still present, and
  **refuses** (routing to manual review) if any are missing — a changed public
  signature means a real port is required, not a provenance bump
* rewrites the provenance block **only** for `auto: true` mappings, which are
  limited to files whose numerical behavior is covered by the live Rcpp fidelity
  suite
* writes `sync/last_sync_report.md` explaining exactly what changed

For `auto: false` mappings and for **unmapped** changed R `src` files, the
workflow does not invent a port. It opens a draft PR / manual-review request
enumerating the files a human must port (and add fidelity coverage for).
Generated Rcpp glue and build files listed under `ignore` do not trigger a
manual-port request.

## Live Rcpp fidelity gate

Before a sync PR may be marked ready, the workflow:

1. installs the incoming R checkout with `R CMD INSTALL`
2. generates JSON truth from live R calls (`tests/fidelity/generate_truth.R`)
3. builds `NNS-core` (with `-DNNSCORE_BUILD_FIDELITY=ON`)
4. compares C++ outputs against the live R truth (`scripts/run_fidelity.py`)

A sync PR is opened as a **draft** unless every changed `src/**` file was cleanly
auto-handled by the adapter **and** the live Rcpp fidelity suite passed. No sync
is marked ready on fidelity failure or when fidelity could not be run.

## Downstream notification

After a core sync merges to `main`, this repository dispatches
`nns-core-updated` to `OVVO-Financial/NNS-python`.

The Python repository then updates its vendored `extern/NNS-core`, builds the
native backend, and runs parity gates.
