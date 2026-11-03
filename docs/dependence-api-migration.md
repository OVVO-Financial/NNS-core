# Dependence API after the R de167f08 sync

R redesigned `NNS.dep` (`src/NNS_dep.cpp`). Dependence is now read from every
level of three recursive mean-split partitions: on x, on y, and the joint
quadrants of (x, y). NNS-core follows R, so `nns::dep_pair` changed signature
and three functions are new. Bindings built against the old API won't compile
until they are updated.

## dep_pair

Old:

```cpp
DepResult dep_pair(const double* x, const double* y, std::size_t n,
                   const uint64_t* quad_xy, const uint64_t* quad_yx, bool asym);
```

New:

```cpp
DepResult dep_pair(const double* x, const double* y, std::size_t n,
                   const std::vector<std::string>& quad_xy,
                   const std::vector<std::string>& quad_yx,
                   const std::vector<std::string>& quad_joint, bool asym = false);
```

- The labels are the prefix-coded quadrant strings that `nns::partition`
  returns ("q", "q1", "q12", ...). The levels are read from the label
  prefixes, so hashed or renumbered cell ids give wrong results.
- `quad_xy` is the x-only partition of (x, y), `quad_yx` the x-only partition
  of (y, x), and `quad_joint` the joint partition of (x, y). With `asym = true`
  only `quad_xy` is read, and the other two may be empty.
- Non-finite x or y, or label vectors of the wrong length, throw
  `std::invalid_argument`.

Most callers should use `nns::dep` instead, which builds the partitions itself.

## New functions

- `DepResult dep(x, y, n, asym = false)` is R's `NNS.dep(x, y, asym)`. It
  builds the three partitions with `nns::partition` using R's settings: order
  `max(1, ceiling(log2 n))`, `obs_req = max(8, n / 8)`, `min_obs_stop = false`,
  noise reduction "off", quadrants only. The joint partition uses no `type`;
  the other two use "XONLY".
- `std::vector<double> dep_levels(x, y, n, quad_xy)` is R's
  `NNS_dep_levels_cpp`: the unfloored dependence of y on x at each distinct
  level of the partition on x. R's `NNS.caus` averages it.
- `double dep_local(x, y, n, quad_xy, quad_yx, asym = true)` is R's
  `NNS_dep_local_cpp`: the finest-cell measure from before the redesign.
  R's `NNS.reg` default order and the `dy.d_` step rule use it, not `NNS.dep`.

## dep_matrix

The signature is unchanged. Each off-diagonal entry [i, j] is now
`dep(column i, column j, asym)`, and the diagonal is fixed at 1,
as R's `NNS.dep(X, asym)` loops over the bivariate path. The old kernel used its
own partition (floor(log2 n) levels, no joint partition), which R has removed.

## Binding checklist

1. Build the three partitions with `nns::partition`, or call `nns::dep`.
2. Pass the label strings as they are, with no hashing.
3. Expose `dep_levels` and `dep_local` where the causation and regression
   ports need them.
4. Regenerate parity values: dependence results change everywhere, not only on
   edge cases.

`nns::seasonality` also gained an optional trailing `max_period` argument.
Existing calls still compile.
