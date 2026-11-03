// tests/cpp/test_distance.cpp
//
// The NNS.distance kernels against a direct computation of their weight blend,
// with the exponential rank weight at rate 1/k (dexp(rank, rate = 1/k)) and the
// log-normal one at sd(1:k) in every path.
//
// SPDX-License-Identifier: GPL-3.0-only
#include "nns/distance.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>

#include "check.hpp"

using nns_test::check_close;

namespace {

void normalise(std::vector<double>& v) {
  const double s = std::accumulate(v.begin(), v.end(), 0.0);
  for (double& e : v) e = s > 0 ? e / s : 0.0;
}

double sd(const std::vector<double>& v) {
  const double m = std::accumulate(v.begin(), v.end(), 0.0) / static_cast<double>(v.size());
  double a = 0.0;
  for (double e : v) a += (e - m) * (e - m);
  return std::sqrt(a / static_cast<double>(v.size() - 1));
}

// X is l x 2 column-major; the ranges include the query point.
double reference(const std::vector<double>& X, const std::vector<double>& y, const double* q, int k) {
  const std::size_t l = y.size();
  double lo[2], hi[2];
  for (int j = 0; j < 2; ++j) {
    lo[j] = hi[j] = q[j];
    for (std::size_t i = 0; i < l; ++i) {
      lo[j] = std::min(lo[j], X[j * l + i]);
      hi[j] = std::max(hi[j], X[j * l + i]);
    }
  }
  std::vector<double> d(l);
  for (std::size_t i = 0; i < l; ++i) {
    d[i] = 0.0;
    for (int j = 0; j < 2; ++j) {
      const double z = (X[j * l + i] - q[j]) / (hi[j] - lo[j]);
      d[i] += std::fabs(z) + z * z;
    }
  }
  std::vector<std::size_t> ord(l);
  std::iota(ord.begin(), ord.end(), 0);
  std::stable_sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b) { return d[a] < d[b]; });

  std::vector<double> s(k), r(k), w(k, 0.0);
  for (int i = 0; i < k; ++i) { s[i] = d[ord[i]]; r[i] = i + 1.0; }
  const double sds = sd(s);
  const double sdr = sd(r);
  std::vector<std::vector<double>> parts(8, std::vector<double>(k));
  for (int i = 0; i < k; ++i) {
    parts[0][i] = 1.0 / k;
    parts[1][i] = std::pow(1.0 + s[i] * s[i] / k, -(k + 1.0) / 2.0);
    parts[2][i] = 1.0 / s[i];
    parts[3][i] = std::exp(-r[i] / k) / k;
    const double lz = std::log(r[k - 1 - i]) / sdr;
    parts[4][i] = std::fabs(-std::log(r[k - 1 - i] * sdr * std::sqrt(2.0 * std::acos(-1.0))) - 0.5 * lz * lz);
    parts[5][i] = 1.0 / (r[i] * r[i]);
    parts[6][i] = std::exp(-0.5 * (s[i] / sds) * (s[i] / sds));
    parts[7][i] = std::exp(-s[i] / (2.0 * sds * sds));
  }
  for (auto& p : parts) {
    normalise(p);
    for (int i = 0; i < k; ++i) w[i] += p[i];
  }
  normalise(w);
  double out = 0.0;
  for (int i = 0; i < k; ++i) out += w[i] * y[ord[i]];
  return out;
}

}  // namespace

int main() {
  const std::size_t l = 80, m = 3;
  std::vector<double> X(2 * l), y(l);
  for (std::size_t i = 0; i < l; ++i) {
    const double t = static_cast<double>(i);
    X[i] = 0.5 + 0.5 * std::sin(1.3 * t + 0.2);
    X[l + i] = 0.5 + 0.5 * std::cos(0.7 * t + 1.1);
    y[i] = std::sin(4.0 * X[i]) + X[l + i] * X[l + i];
  }
  const std::vector<double> Q = {0.37, 0.12, 0.81, 0.61, 0.44, 0.05};  // 3 x 2 column-major

  const int kmax = 10;
  const std::vector<double> path =
      nns::distance_path_parallel(X.data(), l, 2, y.data(), Q.data(), m, kmax, false, 2);
  for (int k : {2, 3, 5, 10}) {
    const std::vector<double> single =
        nns::distance_path_single_parallel(X.data(), l, 2, y.data(), Q.data(), m, k, false, 2);
    for (std::size_t r = 0; r < m; ++r) {
      const double q[2] = {Q[r], Q[m + r]};
      const double want = reference(X, y, q, k);
      const std::string at = " k=" + std::to_string(k) + " row " + std::to_string(r);
      check_close("distance" + at, nns::distance(X.data(), l, 2, y.data(), q, k, false), want, 1e-12);
      check_close("distance_path_single_parallel" + at, single[r], want, 1e-12);
      check_close("distance_path_parallel" + at, path[(k - 1) * m + r], want, 1e-12);
    }
  }
  return nns_test::check_summary("test_distance");
}
