// tests/cpp/smoke_partial_moments.cpp
//
// Minimal dependency-free smoke test for the core partial moment functions.
// Exits non-zero on failure so it can be driven by ctest.
//
// SPDX-License-Identifier: GPL-3.0-only
#include "nns/partial_moments.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

namespace {

int g_failures = 0;

void check_close(const char* label, double got, double want, double tol) {
  if (std::fabs(got - want) > tol) {
    std::fprintf(stderr, "FAIL %s: got %.12g, want %.12g\n", label, got, want);
    ++g_failures;
  } else {
    std::fprintf(stderr, "ok   %s: %.12g\n", label, got);
  }
}

} // namespace

int main() {
  const std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0};
  const double target = 3.0;
  const double tol = 1e-12;

  // Degree-1 lower partial moment: mean(max(target - x, 0)).
  // (2 + 1 + 0 + 0 + 0) / 5 = 0.6
  const double lpm1 = nns::lpm(1.0, target, x.data(), x.size());
  check_close("lpm(1, 3)", lpm1, 0.6, tol);

  // Degree-1 upper partial moment: mean(max(x - target, 0)).
  // (0 + 0 + 0 + 1 + 2) / 5 = 0.6
  const double upm1 = nns::upm(1.0, target, x.data(), x.size());
  check_close("upm(1, 3)", upm1, 0.6, tol);

  // For a symmetric series about the target the two sides are equal.
  check_close("lpm == upm (symmetric)", lpm1, upm1, tol);

  // Degree 0 with points tied with the target: CLPM (all <=), CUPM (all >) and
  // DPM (neither) partition the sample. R 13.3 with the degree-0 tie patch gives
  // 1/2, 1/6, 1/3; unpatched R 13.3 gives 1/2, 1/2, 2/3.
  const double tied[] = {0, -1, 0, 1, 1, -1,   // x
                         0, 0, 1, 1, -1, -1};  // y
  const double origin[] = {0, 0};
  for (bool norm : {true, false}) {
    check_close("degree-0 tied Co.LPM_nD", nns::clpm_nd(tied, 6, 2, origin, 0.0, norm), 0.5, 0.0);
    check_close("degree-0 tied Co.UPM_nD", nns::cupm_nd(tied, 6, 2, origin, 0.0, norm), 1.0 / 6.0, 0.0);
    check_close("degree-0 tied DPM_nD", nns::dpm_nd(tied, 6, 2, origin, 0.0, norm), 1.0 / 3.0, 0.0);
  }

  if (g_failures != 0) {
    std::fprintf(stderr, "%d partial moment smoke check(s) failed\n",
                 g_failures);
    return 1;
  }
  std::fprintf(stderr, "all partial moment smoke checks passed\n");
  return 0;
}
