// tests/cpp/test_stochastic_dominance.cpp
//
// Empty samples are an error in every dominance test, as in R's
// "Each sample needs at least one observation." (13.3 crashed on them).
// The equal-length _uni forms forward to the two-length ones.
//
// SPDX-License-Identifier: GPL-3.0-only
#include "nns/stochastic_dominance.hpp"

#include <vector>

#include "check.hpp"

using nns_test::check;
using nns_test::check_throws;

int main() {
  const std::vector<double> y = {0.0, 1.0};
  const double* none = nullptr;
  check_throws("fsd_uni empty x", [&] { nns::fsd_uni(none, 0, y.data(), y.size(), true); });
  check_throws("fsd_uni empty y", [&] { nns::fsd_uni(y.data(), y.size(), none, 0, false); });
  check_throws("ssd_uni empty x", [&] { nns::ssd_uni(none, 0, y.data(), y.size()); });
  check_throws("ssd_uni empty y", [&] { nns::ssd_uni(y.data(), y.size(), none, 0); });
  check_throws("tsd_uni empty x", [&] { nns::tsd_uni(none, 0, y.data(), y.size()); });
  check_throws("tsd_uni empty y", [&] { nns::tsd_uni(y.data(), y.size(), none, 0); });
  // 0 x 2 matrices
  check_throws("sd_dom_matrix no rows", [&] { nns::sd_dom_matrix(none, 0, 2, 2, true); });
  check_throws("fsd no rows", [&] { nns::fsd(none, 0, 2, true); });
  check_throws("tsd no rows", [&] { nns::tsd(none, 0, 2); });
  const double a[2] = {2.0, 3.0}, b[2] = {0.0, 1.0};
  check(nns::fsd_uni(a, b, 2, true) == 1 && nns::fsd_uni(b, a, 2, false) == 0, "fsd_uni equal length");
  check(nns::ssd_uni(a, b, 2) == 1 && nns::ssd_uni(b, a, 2) == 0, "ssd_uni equal length");
  check(nns::tsd_uni(a, b, 2) == 1 && nns::tsd_uni(b, a, 2) == 0, "tsd_uni equal length");
  return nns_test::check_summary("test_stochastic_dominance");
}
