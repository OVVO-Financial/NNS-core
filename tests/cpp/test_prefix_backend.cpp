// tests/cpp/test_prefix_backend.cpp
//
// The vectorised LPM/UPM kernels switch to a sorted prefix backend when more
// than 32 targets are requested. Its result must not depend on how many
// targets are passed, so every value is checked against the direct scan
// (nns::lpm / nns::upm) on data where the NNS 13.0-13.3 binomial expansion
// cancelled catastrophically (price levels, heavy tails, targets in the tails
// and at order statistics, ties).
//
// SPDX-License-Identifier: GPL-3.0-only
#include "nns/partial_moments.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include "check.hpp"

using nns_test::check;
using nns_test::check_close;

namespace {

// Platform-independent uniform(0, 1) from mt19937_64 raw output.
double unif(std::mt19937_64& g) {
  return (static_cast<double>(g() >> 11) + 0.5) / 9007199254740992.0;
}

std::vector<double> price_path(std::size_t n, std::uint64_t seed) {
  std::mt19937_64 g(seed);
  std::vector<double> x(n);
  double p = 700.0;
  for (std::size_t i = 0; i < n; ++i) {
    p *= 1.0 + 0.02 * (unif(g) - 0.5);
    x[i] = p;
  }
  return x;
}

std::vector<double> heavy_tailed(std::size_t n, std::uint64_t seed) {
  std::mt19937_64 g(seed);
  std::vector<double> x(n);
  for (std::size_t i = 0; i < n; ++i) {
    const double u = unif(g) - 0.5;
    x[i] = 0.01 * u / std::pow(unif(g), 0.4);  // Pareto-like tails
  }
  return x;
}

std::vector<double> targets_for(const std::vector<double>& x) {
  std::vector<double> s = x;
  std::sort(s.begin(), s.end());
  const double lo = s.front(), hi = s.back(), span = hi - lo;
  std::vector<double> t;
  for (int i = 0; i < 120; ++i) t.push_back(lo - 0.1 * span + 1.2 * span * i / 119.0);
  for (std::size_t i = 0; i < s.size(); i += 7) t.push_back(s[i]);  // order stats
  const std::size_t tail = std::min<std::size_t>(40, s.size());
  for (std::size_t i = 1; i <= tail; ++i) {                         // just inside both tails
    t.push_back(std::nextafter(s[s.size() - i], -1e300));
    t.push_back(std::nextafter(s[i - 1], 1e300));
  }
  return t;
}

void check_against_direct(const std::string& name, const std::vector<double>& x,
                          double degree) {
  const std::vector<double> t = targets_for(x);
  const std::size_t m = t.size();
  std::vector<double> lo(m), up(m), lr(m), ur(m);
  nns::lpm_v(degree, t.data(), m, x.data(), x.size(), lo.data());
  nns::upm_v(degree, t.data(), m, x.data(), x.size(), up.data());
  nns::lpm_ratio_v(degree, t.data(), m, x.data(), x.size(), lr.data());
  nns::upm_ratio_v(degree, t.data(), m, x.data(), x.size(), ur.data());
  int negatives = 0, bad_ratio = 0;
  const std::string tag = name + " deg " + std::to_string(degree);
  for (std::size_t i = 0; i < m; ++i) {
    const double dl = nns::lpm(degree, t[i], x.data(), x.size());
    const double du = nns::upm(degree, t[i], x.data(), x.size());
    // The direct scan itself carries ~n*eps summation error.
    check_close(tag + " lpm[" + std::to_string(i) + "]", lo[i], dl, 1e-11);
    check_close(tag + " upm[" + std::to_string(i) + "]", up[i], du, 1e-11);
    if (lo[i] < 0.0 || up[i] < 0.0) ++negatives;
    if (!(lr[i] >= 0.0 && lr[i] <= 1.0) && !std::isnan(lr[i])) ++bad_ratio;
    if (!(ur[i] >= 0.0 && ur[i] <= 1.0) && !std::isnan(ur[i])) ++bad_ratio;
  }
  check(negatives == 0, tag + " no negative partial moments");
  check(bad_ratio == 0, tag + " ratios within [0, 1]");
}

// The answer for one target must be the same whether it is requested alone,
// in a 32-target call (direct path) or in a 33-target call (prefix path).
void check_call_size_invariance(const std::vector<double>& x) {
  std::vector<double> t = targets_for(x);
  t.resize(33);
  for (double degree : {1.0, 2.0, 4.0, 8.0}) {
    std::vector<double> o32(32), o33(33);
    nns::upm_v(degree, t.data() + 1, 32, x.data(), x.size(), o32.data());
    nns::upm_v(degree, t.data(), 33, x.data(), x.size(), o33.data());
    for (std::size_t i = 0; i < 32; ++i)
      check_close("call-size invariance deg " + std::to_string(degree), o33[i + 1], o32[i],
                  1e-11);
  }
}

}  // namespace

int main() {
  const std::vector<double> prices = price_path(3000, 7);
  const std::vector<double> tails = heavy_tailed(3000, 11);
  std::vector<double> ties = heavy_tailed(400, 3);
  for (double& v : ties) v = std::round(v * 1000.0);  // heavy ties
  const std::vector<double> tiny = {2.0, 2.0, 5.0};

  for (double degree : {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 12.0, 20.0, 32.0}) {
    check_against_direct("prices", prices, degree);
    check_against_direct("tails", tails, degree);
    check_against_direct("ties", ties, degree);
    check_against_direct("tiny", tiny, degree);
  }
  // Fractional and > 32 degrees take the direct path inside the vectorised
  // kernel; they must agree as well.
  check_against_direct("prices", prices, 2.5);
  check_against_direct("tails", tails, 33.0);
  check_call_size_invariance(prices);

  // Degree 2 upper tail just below the maximum (NNS 13.3 returned -1.35e-17).
  std::vector<double> near_max(10000);
  std::mt19937_64 g(5);
  for (double& v : near_max) v = 100.0 + (unif(g) - 0.5) * 1e-2;
  std::vector<double> s = near_max;
  std::sort(s.begin(), s.end());
  std::vector<double> t(s.end() - 100, s.end());
  for (double& v : t) v -= 1e-9;
  std::vector<double> u(t.size());
  nns::upm_v(2.0, t.data(), t.size(), near_max.data(), near_max.size(), u.data());
  for (std::size_t i = 0; i < t.size(); ++i) {
    const double d = nns::upm(2.0, t[i], near_max.data(), near_max.size());
    check(u[i] >= 0.0, "near-max degree-2 UPM non-negative");
    check_close("near-max degree-2 UPM", u[i], d, 1e-11);
  }

  // NaN / Inf targets fall back to the direct scan, as upstream.
  std::vector<double> odd(40, 1.0);
  odd[0] = std::numeric_limits<double>::quiet_NaN();
  odd[1] = std::numeric_limits<double>::infinity();
  std::vector<double> o(odd.size());
  nns::lpm_v(1.0, odd.data(), odd.size(), prices.data(), prices.size(), o.data());
  check_close("NaN target matches direct scan", o[0],
              nns::lpm(1.0, odd[0], prices.data(), prices.size()), 0.0);
  check(std::isinf(o[1]), "+Inf target gives +Inf LPM at degree 1");

  return nns_test::check_summary("test_prefix_backend");
}
