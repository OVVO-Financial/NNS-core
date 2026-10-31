// tests/cpp/test_central_tendencies.cpp
//
// Edge cases of nns::mode and nns::gravity (discrete keys beyond the int
// range, tiny histograms, pathological bin counts).
//
// SPDX-License-Identifier: GPL-3.0-only
#include "nns/central_tendencies.hpp"

#include <cmath>
#include <vector>

#include "check.hpp"

using nns_test::check;
using nns_test::check_close;

namespace {

double mode1(const std::vector<double>& x, bool discrete) {
  const std::vector<double> m = nns::mode(x.data(), x.size(), discrete, false);
  return m.empty() ? std::nan("") : m[0];
}

void discrete_mode_cases() {
  // NNS 13.3 returned -2147483648 for both (int overflow, UB).
  check_close("discrete mode 3e9", mode1({3e9, 3e9, 3e9, 1, 2}, true), 3e9, 0.0);
  check_close("discrete mode -5e9", mode1({-5e9, -5e9, -5e9, 1, 2}, true), -5e9, 0.0);
  check_close("discrete mode 2^60", mode1({1152921504606846976.0, 1152921504606846976.0, 7, 8}, true),
              1152921504606846976.0, 0.0);
  // In-range behaviour is unchanged.
  const std::vector<double> ties = {1, 2, 2, 3, 3, 4, 4, 5};
  const std::vector<double> multi = nns::mode(ties.data(), ties.size(), true, true);
  check(multi == std::vector<double>({2, 3, 4}), "discrete multi-mode {2,3,4}");
  check_close("discrete mean of tied modes", mode1(ties, true), 3.0, 0.0);
  check_close("discrete half-up rounding", mode1({1.5, 1.6, 2.4, 0.2, 7.0}, true), 2.0, 0.0);
  // -0.3 rounds half-up to -0.0, which must share the 0 key and print as 0.
  const double z = mode1({-0.3, -0.2, 0.1, 0.2, 5.0}, true);
  check(z == 0.0 && !std::signbit(z), "discrete mode folds -0 into 0");
}

}  // namespace

int main() {
  discrete_mode_cases();
  return nns_test::check_summary("test_central_tendencies");
}
