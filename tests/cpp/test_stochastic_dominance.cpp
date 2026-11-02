// tests/cpp/test_stochastic_dominance.cpp
//
// Empty samples are an error in every dominance test, as in R's
// "Each sample needs at least one observation." (13.3 crashed on them).
// The equal-length _uni forms forward to the two-length ones. The verdicts
// are R's NNS.FSD.uni / NNS.SSD.uni / NNS.TSD.uni after its SD.cpp fixes.
//
// SPDX-License-Identifier: GPL-3.0-only
#include "nns/stochastic_dominance.hpp"

#include <limits>
#include <vector>

#include "check.hpp"

using nns_test::check;
using nns_test::check_throws;

namespace {

using V = std::vector<double>;

int fsd(const V& x, const V& y, bool discrete) {
  return nns::fsd_uni(x.data(), x.size(), y.data(), y.size(), discrete);
}
int ssd(const V& x, const V& y) { return nns::ssd_uni(x.data(), x.size(), y.data(), y.size()); }
int tsd(const V& x, const V& y) { return nns::tsd_uni(x.data(), x.size(), y.data(), y.size()); }

// round(qnorm(ppoints(1e4)), 1) as {10 * value, count}: sorted, so the
// samples match R's even though C++ has no qnorm.
V qnorm_1e4(double scale, double shift) {
  static const int bins[][2] = {
      {-39, 1}, {-36, 1}, {-35, 1}, {-34, 1}, {-33, 2}, {-32, 2}, {-31, 3}, {-30, 5},
      {-29, 6}, {-28, 8}, {-27, 10}, {-26, 14}, {-25, 17}, {-24, 23}, {-23, 28},
      {-22, 36}, {-21, 44}, {-20, 54}, {-19, 66}, {-18, 79}, {-17, 94}, {-16, 111},
      {-15, 129}, {-14, 150}, {-13, 171}, {-12, 195}, {-11, 218}, {-10, 242},
      {-9, 266}, {-8, 289}, {-7, 312}, {-6, 334}, {-5, 352}, {-4, 368}, {-3, 381},
      {-2, 391}, {-1, 397}, {0, 398}, {1, 397}, {2, 391}, {3, 381}, {4, 368},
      {5, 352}, {6, 334}, {7, 312}, {8, 289}, {9, 266}, {10, 242}, {11, 218},
      {12, 195}, {13, 171}, {14, 150}, {15, 129}, {16, 111}, {17, 94}, {18, 79},
      {19, 66}, {20, 54}, {21, 44}, {22, 36}, {23, 28}, {24, 23}, {25, 17},
      {26, 14}, {27, 10}, {28, 8}, {29, 6}, {30, 5}, {31, 3}, {32, 2}, {33, 2},
      {34, 1}, {35, 1}, {36, 1}, {39, 1}};
  V out;
  for (const auto& b : bins)
    for (int k = 0; k < b[1]; ++k) out.push_back(scale * (b[0] / 10.0) + shift);
  return out;
}

}  // namespace

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

  // Each sample is walked with its own length (13.3 used x's for both: 0 0 0).
  check(fsd({4}, {3, 0}, true) == 1, "unequal lengths FSD");
  check(ssd({1, 4}, {3, 1, 3, 0, 2}) == 1, "unequal lengths SSD");
  check(tsd({4, 1}, {0, 3, 3, 4}) == 1, "unequal lengths TSD");

  // Coinciding curves never dominate (13.3: 1 in both directions).
  check(fsd({1, 2}, {1, 1, 2, 2}, false) == 0, "coincide continuous FSD");
  check(ssd({1, 2}, {1, 1, 2, 2}) == 0 && ssd({1, 1, 2, 2}, {1, 2}) == 0, "coincide SSD");
  check(tsd({1, 2}, {1, 1, 2, 2}) == 0, "coincide TSD");
  // Continuous FSD ratios equal at every knot but crossing between them.
  check(fsd({0, 0, 3, 3}, {1, 1, 3, 0}, false) == 0, "knot ratios equal x-y");
  check(fsd({1, 1, 3, 0}, {0, 0, 3, 3}, false) == 0, "knot ratios equal y-x");
  check(fsd({0, 1, 1}, {0, 0, 1}, false) == 1, "midpoint x-y");
  check(fsd({0, 0, 1}, {0, 1, 1}, false) == 0, "midpoint y-x");
  check(fsd({1, 1}, {0, 1}, false) == 1, "midpoint point mass");

  // 0/0 ratio at a point mass is the CDF value k/m, not 0.
  check(fsd({4, 4, 4, 4, 4, 4}, {5, 5, 0, 4, 0, 2}, false) == 0, "0/0 ratio");

  // Curves that cross between sample points.
  check(tsd({-0.85, -1.52, 0.68}, {-0.21, -0.33, -1.67}) == 0, "TSD crossing between knots");
  check(fsd({0.15, 1.55, -0.86, -0.8, 0.25, 0.33}, {-0.15, 1.17, 0.26, -0.04, -2.44, 1.19},
            false) == 0,
        "continuous FSD crossing between knots");

  // Rounding-level differences are ties; the 1e4 samples also need
  // compensated prefix sums to stay inside the tolerance.
  check(ssd({0.1, 0.7}, {0, 0.8}) == 1, "SSD tie tolerance");
  const V zx = qnorm_1e4(0.1, 0.3), zy = qnorm_1e4(0.2, 0.3);
  check(zx.size() == 10000, "qnorm_1e4 size");
  check(ssd(zx, zy) == 1, "SSD 1e4 compensated sums");
  check(tsd(zx, zy) == 1, "TSD 1e4 compensated sums");

  // A common shift keeps the verdict: the tolerances used to grow with the
  // level squared, so +20000 turned "no TSD" into TSD.
  auto shift = [](V v, double s) {
    for (double& e : v) e += s;
    return v;
  };
  const V ax = {-.0085, -.0152, .0068}, ay = {-.0021, -.0033, -.0167};
  check(tsd(shift(ax, 20000), shift(ay, 20000)) == tsd(ax, ay), "TSD +20000 unchanged");
  const V bx = {.01, .07}, by = {0, .08};
  check(tsd(shift(bx, 1e6), shift(by, 1e6)) == 1, "TSD +1e6");
  const V bm = {.01 + 1e6, .07 + 1e6, 0 + 1e6, .08 + 1e6};  // cbind(x, y) + 1e6
  check(nns::tsd(bm.data(), 2, 2) == std::vector<int>({0}), "TSD efficient set +1e6");
  const V cx = {.01, .03}, cy = {0, .04};
  check(ssd(shift(cx, 20000), shift(cy, 20000)) == 1, "SSD +20000");
  check(tsd(shift(cx, 20000), shift(cy, 20000)) == 1, "TSD +20000");
  V dx = {0, .01, .03}, dy = {0, .02, .02}, dz = {0, .02, .02};
  for (int i = 1; i <= 996; ++i) {
    dx.push_back(i);
    dy.push_back(i);
    dz.push_back(i);
  }
  dx.push_back(1000);
  dy.push_back(999);
  dz.push_back(1000);
  for (double s : {20000.0, 1e6}) {
    check(tsd(shift(dx, s), shift(dy, s)) == 0, "TSD long tails, no dominance");
    check(tsd(shift(dz, s), shift(dx, s)) == 1, "TSD long tails, dominance");
  }
  // An infinite value has no midrange to centre on; the tests must still return.
  const double inf = std::numeric_limits<double>::infinity();
  const V ix = {0, 1, inf}, iy = {0, 1, 2};
  check(fsd(ix, iy, false) == 0 && fsd(iy, ix, false) == 0, "continuous FSD with +Inf");
  check(ssd(ix, iy) == 0 && tsd(ix, iy) == 0, "SSD/TSD with +Inf");

  return nns_test::check_summary("test_stochastic_dominance");
}
