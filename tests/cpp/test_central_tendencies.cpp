// tests/cpp/test_central_tendencies.cpp
//
// Edge cases of nns::mode and nns::gravity (discrete keys beyond the int
// range, tiny histograms, pathological bin counts). Expected values are the
// NNS 13.3 results for NNS.mode(x, multi = FALSE) / NNS.gravity(x).
//
// SPDX-License-Identifier: GPL-3.0-only
#include "nns/central_tendencies.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <random>
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
  // Tied discrete modes: the one nearest the gravity, then the lower; never their mean.
  check_close("discrete tied modes {2,3,4}", mode1(ties, true), 3.0, 0.0);
  check_close("discrete tied modes {1,1,3,3}", mode1({1, 1, 3, 3}, true), 1.0, 0.0);
  check_close("discrete tied modes {1,3}", mode1({1, 3}, true), 1.0, 0.0);
  check_close("discrete tied modes {1,2}", mode1({1, 2}, true), 1.0, 0.0);
  check_close("discrete mode {2,2,7}", mode1({2, 2, 7}, true), 2.0, 0.0);
  std::vector<double> eq20;
  for (int v = 0; v < 20; ++v) eq20.insert(eq20.end(), 400, static_cast<double>(v));
  check_close("discrete twenty tied modes", mode1(eq20, true), 9.0, 0.0);
  const std::vector<double> two = {1, 3};
  check(nns::mode(two.data(), 2, true, true) == std::vector<double>({1, 3}),
        "discrete multi-mode of two points");
  check_close("discrete half-up rounding", mode1({1.5, 1.6, 2.4, 0.2, 7.0}, true), 2.0, 0.0);
  // -0.3 rounds half-up to -0.0, which must share the 0 key and print as 0.
  const double z = mode1({-0.3, -0.2, 0.1, 0.2, 5.0}, true);
  check(z == 0.0 && !std::signbit(z), "discrete mode folds -0 into 0");
}

// <= 3 histogram bins: the 7-tap smoother used to read past the end of the
// count vector (ASan heap-buffer-overflow). Values are the NNS 13.3 results.
void tiny_histogram_cases() {
  const std::vector<double> a = {0, 0, 1, 1};
  const std::vector<double> ma = nns::mode(a.data(), a.size(), false, true);
  check(ma == std::vector<double>({0, 1}), "mode multi {0,0,1,1}");
  // Two equally high bins: the one nearer the gravity (0.375).
  check_close("mode {0,0,1,1}", mode1(a, false), 0.0, 0.0);
  check_close("mode {-1,-1,1,1}", mode1({-1, -1, 1, 1}, false), -1.0, 0.0);
  check_close("mode {0,1,1,2}", mode1({0, 1, 1, 2}, false), 1.0, 0.0);
  check_close("mode {1,2,2,5}", mode1({1, 2, 2, 5}, false), 2.0, 0.0);
  check_close("mode {1,1,1,2,9}", mode1({1, 1, 1, 2, 9}, false), 1.0, 0.0);
  check_close("gravity {0,0,1,1}", nns::gravity(a.data(), a.size(), false), 0.375, 0.0);
  const std::vector<double> b = {1, 1, 1, 2, 9};
  check_close("gravity {1,1,1,2,9}", nns::gravity(b.data(), b.size(), false), 1.54375, 1e-15);
}

// Reference: the NNS 13.3 dense-histogram gravity (int64 indices), used to
// show that shortening long runs of empty bins leaves gravity bit-identical.
double dense_gravity(std::vector<double> x) {
  std::sort(x.begin(), x.end());
  const int l = static_cast<int>(x.size());
  auto quart = [&](double& q1, double& q2, double& q3) {
    const double l25 = l * 0.25, l50 = l * 0.5, l75 = l * 0.75;
    if (l % 2 == 0) {
      q1 = x[std::max(1, static_cast<int>(std::floor(l25))) - 1];
      q2 = x[std::max(1, static_cast<int>(std::floor(l50))) - 1];
      q3 = x[std::max(1, static_cast<int>(std::floor(l75))) - 1];
      return;
    }
    auto interp = [&](double p) {
      const int f = std::min(std::max(1, static_cast<int>(std::floor(p))), l);
      const int c = std::min(std::max(1, static_cast<int>(std::ceil(p))), l);
      return x[f - 1] + (p - std::floor(p)) * (x[c - 1] - x[f - 1]);
    };
    q1 = interp(l25);
    const int f50 = std::min(std::max(1, static_cast<int>(std::floor(l50))), l);
    const int c50 = std::min(std::max(1, static_cast<int>(std::ceil(l50))), l);
    q2 = 0.5 * (x[f50 - 1] + x[c50 - 1]);
    q3 = interp(l75);
  };
  double q1, q2, q3;
  quart(q1, q2, q3);
  const double width = (q3 - q1) * std::pow(static_cast<double>(l), -0.5);
  const double origin = x.front();
  const std::int64_t nb =
      static_cast<std::int64_t>(std::floor((x.back() - origin) / width + 1e-12)) + 1;
  std::vector<int> counts(static_cast<std::size_t>(nb), 0);
  for (double v : x) {
    std::int64_t k = static_cast<std::int64_t>(std::floor((v - origin) / width));
    k = std::min<std::int64_t>(std::max<std::int64_t>(k, 0), nb - 1);
    ++counts[static_cast<std::size_t>(k)];
  }
  const int maxc = *std::max_element(counts.begin(), counts.end());
  const std::int64_t ties = std::count(counts.begin(), counts.end(), maxc);
  std::int64_t lo = 0, hi = nb - 1;
  if (ties == 1) {
    const std::int64_t zc = std::find(counts.begin(), counts.end(), maxc) - counts.begin();
    lo = std::max<std::int64_t>(0, zc - 1);
    hi = std::min<std::int64_t>(nb - 1, zc + 1);
  }
  long double num = 0.0L, den = 0.0L;
  for (std::int64_t i = lo; i <= hi; ++i) {
    num += static_cast<long double>(origin + static_cast<double>(i) * width) * counts[i];
    den += counts[i];
  }
  long double s = 0.0L;
  for (double v : x) s += v;
  const double mu = static_cast<double>(s / x.size());
  return 0.25 * (q2 + static_cast<double>(num / den) + mu + 0.5 * (q1 + q3));
}

// Tight IQR plus far outliers: the full histogram has millions of bins, most
// of them empty. NNS 13.3 took 12.8 s and > 2.7 GB for the 1e-7 case.
void pathological_bin_counts() {
  std::mt19937_64 g(9);
  std::normal_distribution<double> nd;
  for (double scale : {1e-3, 1e-4, 1e-5}) {
    std::vector<double> x(1000);
    for (double& v : x) v = 5.0 + scale * nd(g);
    x.push_back(6.0);
    x.push_back(-40.0);
    check_close("gravity long empty runs, scale " + std::to_string(scale),
                nns::gravity(x.data(), x.size(), false), dense_gravity(x), 0.0);
  }
  for (double scale : {1e-7, 1e-12, 1e-300}) {
    std::vector<double> x(1000);
    for (double& v : x) v = 5.0 + scale * nd(g);
    x.push_back(6.0);
    const auto t0 = std::chrono::steady_clock::now();
    const double gv = nns::gravity(x.data(), x.size(), false);
    const double mv = mode1(x, false);
    const double secs =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    check(std::fabs(gv - 5.0) < 0.5 && std::fabs(mv - 5.0) < 1e-3,
          "gravity/mode near 5 for scale " + std::to_string(scale));
    check(secs < 1.0, "pathological binning finishes quickly");
  }
}

}  // namespace

// Twenty equally frequent values give equally high smoothed peaks; the single
// mode is the one nearest the gravity (9.2354492587482149), not their average.
void equal_peak_cases() {
  std::vector<double> x;
  for (int v = 0; v < 20; ++v) x.insert(x.end(), 400, static_cast<double>(v));
  const std::vector<double> all = nns::mode(x.data(), x.size(), false, true);
  check(all.size() == 18 && all[8] == mode1(x, false), "equal peaks, the single mode is a peak");
  check_close("equal peaks, R value", mode1(x, false), 8.9442719099991592, 1e-14);
}

int main() {
  discrete_mode_cases();
  equal_peak_cases();
  tiny_histogram_cases();
  pathological_bin_counts();
  return nns_test::check_summary("test_central_tendencies");
}
