// tests/cpp/test_fast_lm.cpp
//
// fast_lm uses the same long-double two-pass mean as Rcpp's mean(), so the
// coefficients match NNS:::fast_lm bit for bit on data with a large offset.
// Expected values: R 13.3 on x86-64 Linux,
//   i <- 1:500; x <- 1e9 + i * 1e-3; y <- ((i * 7919) %% 1000) / 1000
//   NNS:::fast_lm(x, y)
//
// SPDX-License-Identifier: GPL-3.0-only
#include "nns/fast_lm.hpp"

#include <vector>

#include "check.hpp"

using nns_test::check;
using nns_test::check_close;

int main() {
  std::vector<double> x(500), y(500);
  for (int i = 1; i <= 500; ++i) {
    x[i - 1] = 1e9 + i * 1e-3;
    y[i - 1] = ((i * 7919) % 1000) / 1000.0;
  }
  const nns::FastLmResult r = nns::fast_lm(x.data(), y.data(), x.size());
  // Bit-exact where long double is wider than double (as on the reference
  // platform); otherwise the Rcpp mean itself runs in double precision.
  const double tol = sizeof(long double) > sizeof(double) ? 0.0 : 1e-13;
  check_close("intercept", r.coef[0], -4356025.5750295557, tol);
  check_close("slope", r.coef[1], 0.0043560260754383715, tol);
  check_close("fitted[1]", r.fitted_values[0], 0.50041317194700241, tol);
  check_close("fitted[250]", r.fitted_values[249], 0.50149782188236713, tol);
  check_close("fitted[500]", r.fitted_values[499], 0.50258682854473591, tol);
  check(r.df_residual == 498, "df.residual");

  // Constant x: slope 0, intercept mean(y).
  std::vector<double> cx(4, 2.0), cy = {1.0, 2.0, 3.0, 6.0};
  const nns::FastLmResult c = nns::fast_lm(cx.data(), cy.data(), 4);
  check_close("constant x intercept", c.coef[0], 3.0, 0.0);
  check_close("constant x slope", c.coef[1], 0.0, 0.0);
  return nns_test::check_summary("test_fast_lm");
}
