// tests/cpp/test_internal_functions.cpp
//
// generate_vectors / generate_lin_vectors: in-range results equal the
// NNS 13.3 values, and out-of-range lags no longer read out of bounds,
// underflow size_t or divide by zero. arma_seas_weighting falls back to
// equal weights.
//
// SPDX-License-Identifier: GPL-3.0-only
#include "nns/internal_functions.hpp"

#include <cmath>
#include <vector>

#include "check.hpp"

using nns_test::check;
using nns_test::check_throws;

int main() {
  std::vector<double> x(10);
  for (int i = 0; i < 10; ++i) x[i] = i + 1;

  // NNS:::generate.vectors(as.numeric(1:10), c(3L, 4L, 10L))
  const int lags[3] = {3, 4, 10};
  const nns::TimeSeriesVectors gv = nns::generate_vectors(x.data(), x.size(), lags, 3);
  check(gv.series[0] == std::vector<double>({2, 5, 8}), "lag 3 series");
  check(gv.index[0] == std::vector<int>({1, 2, 3}), "lag 3 index");
  check(gv.series[1] == std::vector<double>({3, 7}), "lag 4 series");
  check(gv.series[2] == std::vector<double>({1}), "lag n series");

  // lag > n: empty component (13.3 read x[n]; the size_t port underflowed).
  const int big[2] = {7, 1000000};
  const nns::TimeSeriesVectors ge = nns::generate_vectors(x.data(), 3, big, 2);
  check(ge.series[0].empty() && ge.index[0].empty(), "lag 7 > n = 3 gives an empty component");
  check(ge.series[1].empty(), "lag 1e6 > n gives an empty component");
  check(nns::generate_vectors(x.data(), 0, lags, 1).series[0].empty(), "n = 0");

  // len (R 13.3): only the first len elements are the active series.
  const nns::TimeSeriesVectors gl = nns::generate_vectors(x.data(), x.size(), lags, 1, 7);
  const nns::TimeSeriesVectors g7 = nns::generate_vectors(x.data(), 7, lags, 1);
  check(gl.series[0] == g7.series[0] && gl.series[0] == std::vector<double>({2, 5}),
        "len = 7 equals the first 7 elements");
  check(nns::generate_vectors(x.data(), x.size(), lags, 1, -1).series[0] == gv.series[0],
        "len = -1 uses the full length");
  check_throws("len > n", [&] { nns::generate_vectors(x.data(), 3, lags, 1, 4); });

  // NNS:::generate.lin.vectors(as.numeric(1:10), 3L, 5L)
  const nns::ForecastVectors lv = nns::generate_lin_vectors(x.data(), x.size(), 3, 5);
  check(lv.series[0] == std::vector<double>({2, 5, 8}), "lin component 1");
  check(lv.series[1] == std::vector<double>({3, 6, 9}), "lin component 2");
  check(lv.series[2] == std::vector<double>({1, 4, 7, 10}), "lin component 3");
  check(lv.forecast_values[0] == std::vector<double>({4, 5}), "lin forecast values 1");
  check(lv.forecast_values[2] == std::vector<double>({5}), "lin forecast values 3");
  check(lv.forecast_index[0] == std::vector<int>({1, 4}), "lin forecast index 1");

  // l > n: empty components, forecast positions ceil(i / l).
  const nns::ForecastVectors le = nns::generate_lin_vectors(x.data(), 3, 7, 2);
  check(le.series[0].empty() && le.series[1].empty(), "lin l > n gives empty components");
  check(le.forecast_values[0] == std::vector<double>({1}), "lin l > n forecast value");

  // l = 0 was SIGFPE in R 13.3 (exit code 136); negative h resized to SIZE_MAX.
  check_throws("lin l = 0", [&] { nns::generate_lin_vectors(x.data(), x.size(), 0, 1); });
  check_throws("lin l < 0", [&] { nns::generate_lin_vectors(x.data(), x.size(), -2, 1); });
  check_throws("lin h < 0", [&] { nns::generate_lin_vectors(x.data(), x.size(), 3, -1); });
  check(nns::generate_lin_vectors(x.data(), x.size(), 3, 0).series.empty(), "lin h = 0");

  // (VarCV - CV) / sqrt(lag), normalised; equal weights when that sums to 0
  // (13.3 returned zeros, so NNS.ARMA forecast 0).
  const double per[3] = {7, 14, 21}, cv[3] = {0.1, 0.2, 0.3}, vcv[3] = {0.3, 0.3, 0.3};
  const nns::ARMAWeights w = nns::arma_seas_weighting(per, cv, vcv, 3);
  const double w0 = 0.2 / std::sqrt(7.0), w1 = 0.1 / std::sqrt(14.0);
  nns_test::check_close("weight 7", w.weights[0], w0 / (w0 + w1), 1e-15);
  nns_test::check_close("weight 14", w.weights[1], w1 / (w0 + w1), 1e-15);
  check(w.weights[2] == 0.0, "weight 21 (CV == VarCV)");
  const nns::ARMAWeights eq = nns::arma_seas_weighting(per, vcv, vcv, 3);
  check(eq.lags == std::vector<double>({7, 14, 21}), "zero-sum lags");
  check(eq.weights == std::vector<double>(3, 1.0 / 3.0), "zero-sum weights are equal");

  return nns_test::check_summary("test_internal_functions");
}
