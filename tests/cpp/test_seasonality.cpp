// tests/cpp/test_seasonality.cpp
//
// Fixed values from R's NNS_seas_cpp (e1068620): the dispersion of each
// component series is its Gini mean difference, a CV only for a series on one
// side of 0, and max_period caps the scan.
//
// SPDX-License-Identifier: GPL-3.0-only
#include "nns/seasonality.hpp"

#include <cmath>
#include <optional>
#include <string>
#include <vector>

#include "check.hpp"

using nns_test::check;
using nns_test::check_close;

namespace {

constexpr double kPi = 3.14159265358979323846;

nns::SeasonalityResult seas(const std::vector<double>& x, std::vector<int> mods = {},
                            bool mod_only = true, std::optional<int> max_period = std::nullopt) {
  return nns::seasonality(x.data(), x.size(), mods.empty() ? nullptr : mods.data(), mods.size(),
                          mod_only, max_period);
}

void check_cvs(const std::string& label, const nns::SeasonalityResult& r,
               const std::vector<double>& want, double vcv) {
  check(r.all_coef_var.size() == want.size(), label + " CV count");
  for (std::size_t i = 0; i < want.size() && i < r.all_coef_var.size(); ++i)
    check_close(label + " CV " + std::to_string(i), r.all_coef_var[i], want[i], 1e-13);
  for (double v : r.all_var_coef_var) check_close(label + " CV(X)", v, vcv, 1e-13);
}

std::vector<double> repeat2(double a, double b, int times) {
  std::vector<double> x;
  for (int i = 0; i < times; ++i) {
    x.push_back(a);
    x.push_back(b);
  }
  return x;
}

}  // namespace

int main() {
  // 10 + sin(2 pi t / 12), t = 1..96.
  std::vector<double> sine(96);
  for (int t = 1; t <= 96; ++t) sine[t - 1] = 10.0 + std::sin(2.0 * kPi * t / 12.0);
  const nns::SeasonalityResult s = seas(sine);
  check(s.best_period == 12, "sine best period 12");
  check(s.periods == std::vector<int>({12, 6, 18, 4, 13, 8, 2, 3}), "sine periods");
  check_cvs("sine", s,
            {1.1462454605814141e-16, 0.066347569912335647, 0.07465376541686071,
             0.074946105682756572, 0.075743740498660603, 0.078352746850154642,
             0.079025266554754767, 0.079106687507153842},
            0.081000439187511694);

  // max_period stops the scan at lag 5; a modulo period beyond it is still scored.
  const nns::SeasonalityResult m5 = seas(sine, {}, true, 5);
  check(m5.periods == std::vector<int>({4, 2, 3}), "max_period 5 periods");
  const nns::SeasonalityResult m5m = seas(sine, {12}, true, 5);
  check(m5m.periods == std::vector<int>({12}), "max_period 5, modulo 12");
  check_cvs("max_period 5, modulo 12", m5m, {1.1462454605814141e-16}, 0.081000439187511694);
  const nns::SeasonalityResult m5a = seas(sine, {12}, false, 5);
  check(m5a.periods == std::vector<int>({12, 4, 2, 3, 1}), "max_period 5, modulo 12, all");
  check(seas(sine, {}, true, 0).periods == std::vector<int>({1}), "max_period 0 scans lag 1");

  // Rounded noisy sine.
  const std::vector<double> noisy = {
      10.9574, 10.6267, 10.8611, 10.7836, 10.3059, 9.8105,  9.6496,  9.1106,  9.0305,  9.572,
      9.5714,  10.5434, 10.9563, 10.9308, 11.3792, 10.9596, 10.3212, 9.9385,  9.499,   9.3316,
      9.168,   9.275,   9.7612,  9.7224,  10.7546, 10.9029, 11.1505, 10.9844, 10.3034, 9.9448,
      9.3258,  9.2777,  9.0221,  9.1183,  9.4159,  9.8876,  10.6995, 10.645,  10.9715, 10.929,
      10.7437, 9.8601,  9.4429,  8.8717,  8.9218,  9.0537,  9.7701,  10.1182, 10.5201, 11.0522,
      10.9475, 10.8645, 10.5734, 10.3414, 9.6447,  9.2302,  8.6864,  9.1976,  9.5332,  9.82,
      10.5153, 10.8979, 11.1087, 11.007,  10.5638, 10.2218, 9.6538,  9.3647,  9.2521,  9.2741,
      9.5865,  9.8155,  10.3769, 10.6927, 10.6721, 10.6009, 10.3222, 9.8885,  9.4875,  9.6185,
      9.0685,  9.1348,  9.5058,  9.9213,  10.3415, 10.8037, 10.9308, 10.8051, 10.1428, 10.1175,
      9.8272,  9.0049,  9.1238,  9.1813,  9.6693,  9.8853};
  const nns::SeasonalityResult ns = seas(noisy);
  check(ns.periods == std::vector<int>({12, 13, 6, 4}), "noisy sine periods");
  check_cvs("noisy sine", ns,
            {0.021470698078126623, 0.070400573968646879, 0.074796433826565797,
             0.08018168620800166},
            0.080879033767251154);

  // Constant series and components: Gini 0, so no lag is strictly below CV(X) = 0.
  for (double v : {5.0, 0.1}) {
    const nns::SeasonalityResult c = seas(std::vector<double>(20, v));
    check(c.best_period == 1 && c.periods == std::vector<int>({1}), "constant series has no period");
    check_cvs("constant", c, {0.0}, 0.0);
  }
  // Every even lag has constant components; n / 5 = 4 caps the scan.
  const nns::SeasonalityResult a1 = seas(repeat2(0.0, 1.0, 10));
  check(a1.periods == std::vector<int>({2, 4}), "alternating periods");
  check_cvs("alternating", a1, {0.0, 0.0}, 0.52631578947368418);
  check_cvs("alternating x10", seas(repeat2(0.0, 10.0, 10)), {0.0, 0.0}, 5.2631578947368425);
  check_cvs("alternating decimal", seas(repeat2(0.1, 1.0, 10)), {0.0, 0.0}, 0.86124401913875581);

  // A linear ramp has no period; across 0 the Gini replaces the CV.
  std::vector<double> up(40), down(40), cross(40);
  for (int t = 1; t <= 40; ++t) {
    up[t - 1] = t;
    down[t - 1] = -t;
    cross[t - 1] = t - 20.5;
  }
  check_cvs("linear", seas(up), {0.66666666666666663}, 0.66666666666666663);
  check_cvs("linear decreasing", seas(down), {0.66666666666666663}, 0.66666666666666663);
  const nns::SeasonalityResult lc = seas(cross);
  check(lc.periods == std::vector<int>({1}), "linear across 0 has no period");
  check_cvs("linear across 0", lc, {13.666666666666666}, 13.666666666666666);

  // n < 5: the R placeholder.
  const nns::SeasonalityResult sh = seas({1, 2, 3, 4});
  check(sh.best_period == 0 && sh.periods == std::vector<int>({0}), "short series");

  return nns_test::check_summary("test_seasonality");
}
