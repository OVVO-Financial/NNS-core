// tests/cpp/test_dependence.cpp
//
// NNS.dep and its kernels against R e1068620: dep / dep_matrix (NNS.dep),
// dep_pair (NNS_dep_pair_cpp) on hand-made cells, dep_levels and dep_local.
//
// SPDX-License-Identifier: GPL-3.0-only
#include "nns/dependence.hpp"
#include "nns/partition.hpp"

#include <cmath>
#include <string>
#include <vector>

#include "check.hpp"

using nns_test::check;
using nns_test::check_close;

namespace {

constexpr double kTol = 1e-13;

void check_result(const std::string& tag, const nns::DepResult& r, double correlation,
                  double dependence) {
  check_close(tag + " correlation", r.correlation, correlation, kTol);
  check_close(tag + " dependence", r.dependence, dependence, kTol);
}

// dep_pair with every partition given the same labels.
void check_cells(const std::string& name, const std::vector<double>& x,
                 const std::vector<double>& y, const std::vector<std::string>& cells,
                 double correlation, double dependence) {
  for (bool asym : {false, true}) {
    check_result(name + (asym ? " asym" : ""),
                 nns::dep_pair(x.data(), y.data(), x.size(), cells, cells, cells, asym),
                 correlation, dependence);
  }
}

}  // namespace

int main() {
  // x = cos(0.37 i) + 0.01 i, y = x^2 + 0.3 sin(1.3 i), i = 1..80.
  std::vector<double> x(80), y(80), z(80), neg_x(80), neg_y(80);
  for (int i = 1; i <= 80; ++i) {
    x[i - 1] = std::cos(0.37 * i) + 0.01 * i;
    y[i - 1] = x[i - 1] * x[i - 1] + 0.3 * std::sin(1.3 * i);
    z[i - 1] = std::sin(2.1 * i);
    neg_x[i - 1] = -x[i - 1];
    neg_y[i - 1] = -y[i - 1];
  }
  const double r_xy = 0.20162838991094528, d_xy = 0.77708974213677307;
  check_result("dep(x, y)", nns::dep(x.data(), y.data(), 80), r_xy, d_xy);
  check_result("dep(x, y) asym", nns::dep(x.data(), y.data(), 80, true), r_xy, d_xy);
  // Symmetric: the same three partitions whichever way round; -y flips the sign.
  check_result("dep(y, x)", nns::dep(y.data(), x.data(), 80), r_xy, d_xy);
  check_result("dep(x, -y)", nns::dep(x.data(), neg_y.data(), 80), -r_xy, d_xy);
  // -x moves the mean splits of x, so R gives a different value.
  check_result("dep(-x, y)", nns::dep(neg_x.data(), y.data(), 80), -0.17944931042143974,
               0.75014070202643035);

  std::vector<double> X = x;
  X.insert(X.end(), y.begin(), y.end());
  X.insert(X.end(), z.begin(), z.end());
  const nns::DepMatrixResult m = nns::dep_matrix(X.data(), 80, 3, false, 1);
  const std::vector<double> mc = {1, r_xy, 0.15198594271078439, r_xy, 1, 0.13104291789955438,
                                  0.15198594271078439, 0.13104291789955438, 1};
  const std::vector<double> md = {1, d_xy, 0.28860025468055317, d_xy, 1, 0.13104291789955438,
                                  0.28860025468055317, 0.13104291789955438, 1};
  for (std::size_t k = 0; k < 9; ++k) {
    check_close("matrix correlation " + std::to_string(k), m.correlation[k], mc[k], kTol);
    check_close("matrix dependence " + std::to_string(k), m.dependence[k], md[k], kTol);
  }
  // asym: entry [i, j] is dep(column i, column j, asym = TRUE).
  const nns::DepMatrixResult ma = nns::dep_matrix(X.data(), 80, 3, true, 2);
  const std::vector<double> ac = {1, 0.74654375385431393, 0.15198594271078439, r_xy, 1,
                                  0.13104291789955438, -0.097324036522073221,
                                  -0.03592012003487987, 1};
  const std::vector<double> ad = {1, 0.74654375385431393, 0.28860025468055317, d_xy, 1,
                                  0.13104291789955438, 0.14257708926470505, 0.03592012003487987, 1};
  for (std::size_t k = 0; k < 9; ++k) {
    check_close("asym matrix correlation " + std::to_string(k), ma.correlation[k], ac[k], kTol);
    check_close("asym matrix dependence " + std::to_string(k), ma.dependence[k], ad[k], kTol);
  }

  // A step function of a discrete x: each x-constant cell has y constant.
  std::vector<double> xd, yd;
  for (int r = 0; r < 8; ++r)
    for (int v = 1; v <= 10; ++v) {
      xd.push_back(v);
      yd.push_back(v > 5 ? 1.0 : 0.0);
    }
  check_result("discrete step", nns::dep(xd.data(), yd.data(), 80), 1.0, 1.0);
  check_result("discrete step asym", nns::dep(xd.data(), yd.data(), 80, true), 1.0, 1.0);
  check_result("constant y", nns::dep(x.data(), std::vector<double>(80, 2.0).data(), 80), 0.0, 0.0);
  std::vector<double> bad = x;
  bad[3] = std::nan("");
  nns_test::check_throws("NaN input", [&] { nns::dep(bad.data(), y.data(), 80); });

  // Every point discordant in one cell: perfect negative dependence.
  std::vector<double> c(100), neg(100);
  for (int i = 0; i < 100; ++i) {
    c[i] = i + 1;
    neg[i] = -(i + 1);
  }
  check_cells("x vs -x", c, neg, std::vector<std::string>(100, "q"), -1.0, 1.0);
  // Binary self-dependence in cells of 51 and 49 points: 49 * (1/49) is not 1.
  std::vector<double> b(100, 0.0);
  std::vector<std::string> bcells(100, "q1");
  for (int i = 51; i < 100; ++i) {
    b[i] = 1.0;
    bcells[i] = "q2";
  }
  check_cells("binary self", b, b, bcells, 1.0, 1.0);
  // Points on the cell mean count as co-lower only.
  const double pattern[9] = {0, 1, -1, 0, -1, 1, 1, 0, -1};
  std::vector<double> tx(99), ty(99);
  for (int i = 0; i < 99; ++i) {
    tx[i] = (i % 3) - 1.0;
    ty[i] = pattern[i % 9];
  }
  check_cells("ties at the mean", tx, ty, std::vector<std::string>(99, "q"), -0.21046138763212088,
              0.21046138763212088);

  // Levels and the local measure on the partition NNS.dep builds for (x, y).
  const std::vector<std::string> qa =
      nns::partition(x.data(), y.data(), 80, std::string("XONLY"), 7, 10, false, "off", true).quadrant;
  const std::vector<std::string> qb =
      nns::partition(y.data(), x.data(), 80, std::string("XONLY"), 7, 10, false, "off", true).quadrant;
  const std::vector<double> lv = nns::dep_levels(x.data(), y.data(), 80, qa);
  const std::vector<double> want = {0.74654375385431404, 0.77708974213677307, 0.49841564287980855,
                                    0.23240655071524363, 0.18649740880281743};
  check(lv.size() == want.size(), "dep_levels count");
  for (std::size_t k = 0; k < want.size() && k < lv.size(); ++k)
    check_close("dep_levels " + std::to_string(k), lv[k], want[k], kTol);
  for (bool asym : {false, true})
    check_close(std::string("dep_local") + (asym ? " asym" : ""),
                nns::dep_local(x.data(), y.data(), 80, qa, qb, asym), 0.11550010474265658, kTol);
  nns_test::check_throws("dep_pair label length", [&] {
    nns::dep_pair(x.data(), y.data(), 80, qa, {}, {}, false);
  });

  return nns_test::check_summary("test_dependence");
}
