// tests/cpp/test_dependence.cpp
//
// NNS.dep and its kernels against R 499f82b1: dep / dep_matrix (NNS.dep),
// dep_pair (NNS_dep_pair_cpp) on hand-made cells, dep_levels and dep_local.
//
// SPDX-License-Identifier: GPL-3.0-only
#include "nns/dependence.hpp"
#include "nns/partition.hpp"

#include <cmath>
#include <limits>
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
  const double r_xy = 0.74654375385431426, d_xy = 0.74654375385431426;
  check_result("dep(x, y)", nns::dep(x.data(), y.data(), 80), r_xy, d_xy);
  check_result("dep(x, y) asym", nns::dep(x.data(), y.data(), 80, true), r_xy, d_xy);
  // Symmetric: the same three partitions whichever way round; -y flips the sign.
  check_result("dep(y, x)", nns::dep(y.data(), x.data(), 80), r_xy, d_xy);
  check_result("dep(x, -y)", nns::dep(x.data(), neg_y.data(), 80), -r_xy, d_xy);
  // -x moves the gravity splits of x, so R gives a different value.
  check_result("dep(-x, y)", nns::dep(neg_x.data(), y.data(), 80), -0.17944931042143963,
               0.75014070202643024);

  std::vector<double> X = x;
  X.insert(X.end(), y.begin(), y.end());
  X.insert(X.end(), z.begin(), z.end());
  const nns::DepMatrixResult m = nns::dep_matrix(X.data(), 80, 3, false, 1);
  const std::vector<double> mc = {1, r_xy, 0.075618087484626714, r_xy, 1, 0.0069948492176366058,
                                  0.075618087484626714, 0.0069948492176366058, 1};
  const std::vector<double> md = {1, d_xy, 0.20434399892988306, d_xy, 1, 0.10254377270875774,
                                  0.20434399892988306, 0.10254377270875774, 1};
  for (std::size_t k = 0; k < 9; ++k) {
    check_close("matrix correlation " + std::to_string(k), m.correlation[k], mc[k], kTol);
    check_close("matrix dependence " + std::to_string(k), m.dependence[k], md[k], kTol);
  }
  // asym: entry [i, j] is dep(column i, column j, asym = TRUE).
  const nns::DepMatrixResult ma = nns::dep_matrix(X.data(), 80, 3, true, 2);
  const std::vector<double> ac = {1, 0.74654375385431393, -0.086388021856362518, r_xy, 1, 0,
                                  -0.11435954175833232, 0, 1};
  const std::vector<double> ad = {1, 0.74654375385431393, 0.16803718110241508, d_xy, 1, 0,
                                  0.1716769006516293, 0, 1};
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

  // Tied data: argument order, row order and the decimal image don't move the
  // result (R set.seed(208); n <- sample(10:80, 1); sample(0:5, n, TRUE) twice).
  const std::vector<double> tx8 = {2, 1, 0, 0, 4, 2, 2, 2, 0, 4, 2, 1, 2, 4, 2, 0, 1, 4, 3, 1, 4, 4};
  const std::vector<double> ty8 = {3, 4, 4, 3, 2, 3, 5, 0, 1, 3, 1, 3, 2, 4, 0, 3, 2, 2, 4, 2, 3, 3};
  const std::size_t n8 = tx8.size();
  const double r8 = -0.18989026648936866;
  std::vector<double> rx8(tx8.rbegin(), tx8.rend()), ry8(ty8.rbegin(), ty8.rend());
  std::vector<double> dx8(n8), dy8(n8), hx8(n8), hy8(n8);
  for (std::size_t i = 0; i < n8; ++i) {
    dx8[i] = tx8[i] / 10.0;
    dy8[i] = ty8[i] / 10.0;
    hx8[i] = 1e9 + tx8[i] / 10.0;
    hy8[i] = 1e9 + ty8[i] / 10.0;
  }
  check_result("tied", nns::dep(tx8.data(), ty8.data(), n8), r8, -r8);
  check_result("tied, swapped", nns::dep(ty8.data(), tx8.data(), n8), r8, -r8);
  check_result("tied, rows reversed", nns::dep(rx8.data(), ry8.data(), n8), r8, -r8);
  check_result("tied / 10", nns::dep(dx8.data(), dy8.data(), n8), r8, -r8);
  check_result("tied / 10 + 1e9", nns::dep(hx8.data(), hy8.data(), n8), -0.18989022573686976,
               0.18989022573686976);
  // Distinct values at a high level are not ties.
  std::vector<double> ts(100), nts(100);
  for (int i = 0; i < 100; ++i) {
    ts[i] = 1.7e9 + (i + 1) / 1e6;
    nts[i] = -ts[i];
  }
  check_result("1.7e9 + t vs its negation", nns::dep(ts.data(), nts.data(), 100), -1.0, 1.0);

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
  const std::vector<double> want = {0.74654375385431404, 0.69123054502417314, 0.48542447845329173,
                                    0.37834727124745682, 0.2541022165228275};
  check(lv.size() == want.size(), "dep_levels count");
  for (std::size_t k = 0; k < want.size() && k < lv.size(); ++k)
    check_close("dep_levels " + std::to_string(k), lv[k], want[k], kTol);
  for (bool asym : {false, true})
    check_close(std::string("dep_local") + (asym ? " asym" : ""),
                nns::dep_local(x.data(), y.data(), 80, qa, qb, asym), 0.26017793115176047, kTol);
  // NaN used to loop forever in the tied-value scans; every entry point rejects it.
  std::vector<double> tied_nan = xd, tied_inf = xd;
  tied_nan[9] = std::nan("");
  tied_inf[9] = std::numeric_limits<double>::infinity();
  const std::vector<std::string> one(80, "q");
  for (const std::vector<double>* v : {&tied_nan, &tied_inf}) {
    const std::string what = (v == &tied_nan) ? " NaN" : " Inf";
    for (bool asym : {false, true}) {
      nns_test::check_throws("dep_pair x" + what, [&] {
        nns::dep_pair(v->data(), yd.data(), 80, one, one, one, asym);
      });
      nns_test::check_throws("dep_pair y" + what, [&] {
        nns::dep_pair(xd.data(), v->data(), 80, one, one, one, asym);
      });
      nns_test::check_throws("dep_local" + what, [&] {
        nns::dep_local(xd.data(), v->data(), 80, one, one, asym);
      });
    }
    nns_test::check_throws("dep_levels" + what, [&] { nns::dep_levels(v->data(), yd.data(), 80, one); });
    nns_test::check_throws("dep_matrix" + what, [&] {
      std::vector<double> M = xd;
      M.insert(M.end(), v->begin(), v->end());
      nns::dep_matrix(M.data(), 80, 2, false, 1);
    });
  }
  nns_test::check_throws("dep_pair label length", [&] {
    nns::dep_pair(x.data(), y.data(), 80, qa, {}, {}, false);
  });

  return nns_test::check_summary("test_dependence");
}
