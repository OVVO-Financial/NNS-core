// tests/cpp/test_degree_validation.cpp
//
// Negative (and NaN) degrees are rejected by every partial moment kernel.
// NNS 13.3 silently returned the degree-0 count for LPM(-1, 0, c(-2, 1)).
//
// SPDX-License-Identifier: GPL-3.0-only
#include "nns/partial_moments.hpp"

#include <cmath>
#include <limits>
#include <vector>

#include "check.hpp"

using nns_test::check_close;
using nns_test::check_throws;

int main() {
  const std::vector<double> x = {-2.0, 1.0, 0.5, 3.0};
  const std::vector<double> y = {1.0, -1.0, 2.0, 0.0};
  std::vector<double> t(40, 0.25);
  std::vector<double> out(40);
  std::vector<double> data(x);
  data.insert(data.end(), y.begin(), y.end());  // 4 x 2 column-major
  const double tgt[2] = {0.0, 0.0};
  const double nan = std::numeric_limits<double>::quiet_NaN();

  for (double bad : {-1.0, -0.5, nan}) {
    check_throws("lpm", [&] { nns::lpm(bad, 0.0, x.data(), x.size()); });
    check_throws("upm", [&] { nns::upm(bad, 0.0, x.data(), x.size()); });
    check_throws("lpm_v direct", [&] { nns::lpm_v(bad, t.data(), 2, x.data(), x.size(), out.data()); });
    check_throws("upm_v prefix", [&] { nns::upm_v(bad, t.data(), 40, x.data(), x.size(), out.data()); });
    check_throws("lpm_ratio_v", [&] { nns::lpm_ratio_v(bad, t.data(), 40, x.data(), x.size(), out.data()); });
    check_throws("upm_ratio_v", [&] { nns::upm_ratio_v(bad, t.data(), 3, x.data(), x.size(), out.data()); });
    check_throws("co_lpm x", [&] { nns::co_lpm(bad, 1.0, x.data(), y.data(), 4, 4, 0.0, 0.0); });
    check_throws("co_upm y", [&] { nns::co_upm(1.0, bad, x.data(), y.data(), 4, 4, 0.0, 0.0); });
    check_throws("d_lpm", [&] { nns::d_lpm(bad, 1.0, x.data(), y.data(), 4, 4, 0.0, 0.0); });
    check_throws("d_upm", [&] { nns::d_upm(1.0, bad, x.data(), y.data(), 4, 4, 0.0, 0.0); });
    check_throws("co_lpm_v", [&] {
      nns::co_lpm_v(bad, 1.0, x.data(), y.data(), 4, 4, t.data(), 40, t.data(), 40, out.data());
    });
    check_throws("co_upm_v", [&] {
      nns::co_upm_v(1.0, bad, x.data(), y.data(), 4, 4, t.data(), 40, t.data(), 40, out.data());
    });
    check_throws("d_lpm_v", [&] {
      nns::d_lpm_v(bad, 1.0, x.data(), y.data(), 4, 4, t.data(), 40, t.data(), 40, out.data());
    });
    check_throws("d_upm_v", [&] {
      nns::d_upm_v(1.0, bad, x.data(), y.data(), 4, 4, t.data(), 40, t.data(), 40, out.data());
    });
    check_throws("clpm_nd", [&] { nns::clpm_nd(data.data(), 4, 2, tgt, bad, true); });
    check_throws("cupm_nd", [&] { nns::cupm_nd(data.data(), 4, 2, tgt, bad, false); });
    check_throws("dpm_nd", [&] { nns::dpm_nd(data.data(), 4, 2, tgt, bad, true); });
    check_throws("clpm_nd_batch", [&] {
      nns::clpm_nd_batch(data.data(), 4, 2, tgt, 1, bad, true, out.data());
    });
    check_throws("pm_matrix lpm", [&] { nns::pm_matrix(bad, 1.0, tgt, data.data(), 4, 2, false, false); });
    check_throws("pm_matrix upm", [&] { nns::pm_matrix(1.0, bad, tgt, data.data(), 4, 2, false, false); });
  }

  // Valid degrees still work, including 0 and fractional.
  check_close("lpm degree 0", nns::lpm(0.0, 0.0, x.data(), x.size()), 0.25, 0.0);
  check_close("lpm degree 0.5", nns::lpm(0.5, 0.0, x.data(), x.size()), std::sqrt(2.0) / 4.0,
              1e-15);
  return nns_test::check_summary("test_degree_validation");
}
