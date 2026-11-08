# NNS formal proofs (Lean 4 + Mathlib)

Machine-checked statements about the NNS partial-moment kernels in
`src/partial_moments.cpp`. Definitions copy the C++ conventions exactly
(LPM counts ties, UPM is strict, degree 0 counts observations), so a theorem
here is a theorem about what the code computes, for natural-number degrees and
exact real arithmetic.

| File | Contents |
|---|---|
| `NNSProofs/PartialMoments.lean` | `LPM`, `UPM`, `mean`, `centralMoment`; master identity `UPM(d) + (-1)^d LPM(d) = E[(x-t)^d]`; degree 0/1/2 corollaries; LPM(0) is the empirical CDF; monotonicity; nonnegativity |
| `NNSProofs/CoPartialMoments.lean` | `CoLPM`, `CoUPM`, `DLPM`, `DUPM`; covariance decomposition `CoUPM + CoLPM - DLPM - DUPM = Cov`; degree-0 quadrant partition sums to 1; symmetry and self-moment lemmas |

## Build

```sh
cd proofs
lake exe cache get   # prebuilt Mathlib, ~5 GB, avoids hours of compiling
lake build
```

Toolchain and Mathlib are pinned together (`lean-toolchain`, `lakefile.toml`,
`lake-manifest.json`). Bump them together.

## Scope and gaps

* Natural-number degrees only. Real degrees (`std::pow`) need `Real.rpow`.
* Exact reals, not IEEE doubles: floating-point error is out of scope.
* Equal-length samples for co-moments (the C++ also handles unequal lengths
  by truncating to the shorter sample and dividing by the longer).
* `sorry` is rejected in CI.
