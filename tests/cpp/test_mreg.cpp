// tests/cpp/test_mreg.cpp
//
// Self-consistency checks for the ported M.reg / stack kernels (the numbers
// themselves are checked against live R by the fidelity suite).
//
// SPDX-License-Identifier: GPL-3.0-only
#include "nns/mreg.hpp"

#include <cmath>
#include <random>
#include <vector>

#include "check.hpp"

using nns_test::check;
using nns_test::check_close;
using nns_test::check_throws;

namespace {

double unif(std::mt19937_64& g) {
  return (static_cast<double>(g() >> 11) + 0.5) / 9007199254740992.0;
}

}  // namespace

int main() {
  std::mt19937_64 g(17);
  const std::size_t n = 80, p = 3, m = 25;
  std::vector<double> rpm(n * p), yhat(n), Xt(m * p), mins(p, 0.0), maxs(p, 1.0);
  for (double& v : rpm) v = std::round(unif(g) * 1000.0) / 1000.0;
  for (double& v : yhat) v = unif(g) * 4.0 - 2.0;
  for (double& v : Xt) v = std::round(unif(g) * 1000.0) / 1000.0;
  for (std::size_t j = 0; j < p; ++j) Xt[0 + j * m] = rpm[3 + j * n];  // exact hit on row 3

  // Every column k of the path equals the single-k kernel, for every metric,
  // and the result does not depend on the thread count.
  for (int dist = 0; dist <= 3; ++dist) {
    const int kmax = 10;
    const std::vector<double> path = nns::mreg_predict_path(
        rpm.data(), n, p, yhat.data(), Xt.data(), m, kmax, dist, mins.data(), maxs.data(), false, 1);
    const std::vector<double> path4 = nns::mreg_predict_path(
        rpm.data(), n, p, yhat.data(), Xt.data(), m, kmax, dist, mins.data(), maxs.data(), false, 4);
    check(path == path4, "path independent of thread count, dist " + std::to_string(dist));
    for (int k = 1; k <= kmax; ++k) {
      const std::vector<double> one = nns::mreg_predict(
          rpm.data(), n, p, yhat.data(), Xt.data(), m, k, dist, mins.data(), maxs.data(), false, 1);
      for (std::size_t r = 0; r < m; ++r)
        check_close("path column k == single k, dist " + std::to_string(dist), path[r + (k - 1) * m],
                    one[r], 0.0);
    }
  }
  // k = 1 on an exact hit returns that regression point's value.
  const std::vector<double> k1 = nns::mreg_predict(rpm.data(), n, p, yhat.data(), Xt.data(), m, 1,
                                                   nns::kMregDistNNS, mins.data(), maxs.data(), false);
  check_close("k = 1 exact hit", k1[0], yhat[3], 0.0);
  // k larger than n is clipped to n.
  const std::vector<double> big = nns::mreg_predict(rpm.data(), n, p, yhat.data(), Xt.data(), m,
                                                    1000, 0, mins.data(), maxs.data(), false);
  const std::vector<double> kn = nns::mreg_predict(rpm.data(), n, p, yhat.data(), Xt.data(), m,
                                                   static_cast<int>(n), 0, mins.data(), maxs.data(), false);
  check(big == kn, "k > n clipped to n");
  // Weighted class mode returns one of the classes.
  std::vector<double> cls(n);
  for (std::size_t i = 0; i < n; ++i) cls[i] = static_cast<double>(i % 3);
  const std::vector<double> c5 = nns::mreg_predict(rpm.data(), n, p, cls.data(), Xt.data(), m, 5, 0,
                                                   mins.data(), maxs.data(), true);
  bool classes_ok = true;
  for (double v : c5) classes_ok = classes_ok && (v == 0.0 || v == 1.0 || v == 2.0);
  check(classes_ok, "class predictions are class labels");
  check_throws("k < 1", [&] {
    nns::mreg_predict(rpm.data(), n, p, yhat.data(), Xt.data(), m, 0, 0, mins.data(), maxs.data(), false);
  });
  check_throws("bad distance code", [&] {
    nns::mreg_predict(rpm.data(), n, p, yhat.data(), Xt.data(), m, 3, 4, mins.data(), maxs.data(), false);
  });

  // 1-D fast path (duplicated coordinate): path column k == single k too.
  std::vector<double> r1(n * 2), q1(m * 2);
  for (std::size_t i = 0; i < n; ++i) r1[i] = r1[i + n] = unif(g);
  for (std::size_t i = 0; i < m; ++i) q1[i] = q1[i + m] = unif(g);
  const double mn2[2] = {0.0, 0.0}, mx2[2] = {1.0, 2.0};
  const std::vector<double> p1 =
      nns::mreg_predict_path(r1.data(), n, 2, yhat.data(), q1.data(), m, 6, 0, mn2, mx2, false);
  for (int k = 1; k <= 6; ++k) {
    const std::vector<double> o =
        nns::mreg_predict(r1.data(), n, 2, yhat.data(), q1.data(), m, k, 0, mn2, mx2, false);
    for (std::size_t r = 0; r < m; ++r) check_close("1-D path vs single k", p1[r + (k - 1) * m], o[r], 0.0);
  }

  // A non-duplicated query turns the 1-D shortcut off for the whole batch, so
  // the shortcut must pick the same neighbours and ties as the generic path.
  // Under NNS, (.07, .07) is equidistant from .01 and .13 only after range scaling.
  {
    const double tx[] = {.01, .13, 1, .01, .13, 1}, ty[] = {0, 100, 999};
    const double tq1[] = {.07, .07}, tq2[] = {.07, .2, .07, .3};
    const double lo[] = {.01, .01}, hi[] = {1, 1};
    std::vector<double> gx(n * 2), one(m * 2), two((m + 1) * 2);
    for (std::size_t i = 0; i < n; ++i) gx[i] = gx[i + n] = std::round(unif(g) * 50.0) / 50.0 + .01;
    for (std::size_t i = 0; i < m; ++i)
      one[i] = one[i + m] = two[i] = two[i + m + 1] = std::round(unif(g) * 100.0) / 100.0;
    two[m] = .2;
    two[2 * m + 1] = .3;
    for (int dist = 0; dist <= 2; ++dist) {
      const std::string tag = "1-D vs batched generic, dist " + std::to_string(dist);
      const std::vector<double> a = nns::mreg_predict(tx, 3, 2, ty, tq1, 1, 1, dist, lo, hi, false);
      const std::vector<double> b = nns::mreg_predict(tx, 3, 2, ty, tq2, 2, 1, dist, lo, hi, false);
      check_close(tag + " tie", a[0], b[0], 0.0);
      if (dist == nns::kMregDistNNS) check_close("R reference helper", a[0], 50.0, 0.0);
      // Dividing by the range, as R's reference does, puts .03 strictly nearer
      // (.26, .26) than .49; multiplying by 1/range ties them (50).
      const double dx[] = {.03, .49, 1, .03, .49, 1}, dq[] = {.26, .26}, dlo[] = {.03, .03};
      const std::vector<double> dv = nns::mreg_predict(dx, 3, 2, ty, dq, 1, 1, dist, dlo, hi, false);
      if (dist == nns::kMregDistNNS) check_close("divide by the range", dv[0], 100.0, 0.0);
      for (bool is_cls : {false, true}) {
        const std::vector<double>& yy = is_cls ? cls : yhat;
        const std::vector<double> pa = nns::mreg_predict_path(gx.data(), n, 2, yy.data(), one.data(), m,
                                                              8, dist, lo, hi, is_cls);
        const std::vector<double> pb = nns::mreg_predict_path(gx.data(), n, 2, yy.data(), two.data(),
                                                              m + 1, 8, dist, lo, hi, is_cls);
        for (int k = 1; k <= 8; ++k) {
          const std::vector<double> sa = nns::mreg_predict(gx.data(), n, 2, yy.data(), one.data(), m, k,
                                                           dist, lo, hi, is_cls);
          for (std::size_t r = 0; r < m; ++r) {
            check_close(tag + " path", pa[r + (k - 1) * m], pb[r + (k - 1) * (m + 1)], 0.0);
            check_close(tag + " single", sa[r], pa[r + (k - 1) * m], 0.0);
          }
        }
      }
    }
  }

  // A query about 58 from every regression point, where the Student-t weights
  // underflow. R 13.3 NNS:::NNS_mreg_predict_v2_cpp(x, y, matrix(0, 1, 60),
  // 1000L, 0L, rep(0, 60), rep(1, 60), FALSE, 1L) on the same x, y returns
  // 7.7618903035977533.
  {
    const std::size_t fn = 1000, fp = 60;
    std::vector<double> fx(fn * fp, 0.0), fy(fn, 0.0), fq(fp, 0.0), lo(fp, 0.0), hi(fp, 1.0);
    for (std::size_t i = 0; i < fn; ++i) {
      for (std::size_t j = 0; j < 29; ++j) fx[i + ((i + j) % fp) * fn] = 1.0;
      const double d = 59.57 - static_cast<double>(i) / 999.0;
      fx[i + ((i + 29) % fp) * fn] = (std::sqrt(1.0 + 4.0 * (d - 58.0)) - 1.0) / 2.0;
    }
    fy[fn - 1] = 100.0;
    const std::vector<double> far = nns::mreg_predict(fx.data(), fn, fp, fy.data(), fq.data(), 1,
                                                      static_cast<int>(fn), nns::kMregDistNNS,
                                                      lo.data(), hi.data(), false, 1);
    check_close("far query, underflowing t weights", far[0], 7.7618903035977533, 1e-12);
  }

  // mreg_setup: rows grouped by findInterval ids, groups sorted as strings.
  const std::vector<double> X = {0.1, 0.2, 0.9, 0.95, 0.5,   // column 1
                                 1.0, 1.0, 3.0, 3.0, 2.0};   // column 2
  const std::vector<double> y = {1, 3, 10, 20, 7};
  const std::vector<std::vector<double>> bounds = {{0.0, 0.5, 1.0}, {1.0, 2.0, 3.0}};
  const nns::MregSetupResult s = nns::mreg_setup(X.data(), 5, 2, y.data(), bounds, nns::kMregMean, false);
  check(s.ids == std::vector<std::string>({"1.1", "1.1", "2.2", "2.2", "2.2"}), "setup ids");
  check(s.row_ids == std::vector<std::string>({"1.1", "2.2"}), "setup row ids");
  check(s.n_groups == 2 && s.rpm.size() == 6, "setup shape");
  check_close("setup mean x1 group 1", s.rpm[0], 0.15, 1e-15);
  check_close("setup mean y group 2", s.rpm[2 * 2 + 1], (10.0 + 20.0 + 7.0) / 3.0, 1e-15);
  const nns::MregSetupResult sm = nns::mreg_setup(X.data(), 5, 2, y.data(), bounds, nns::kMregMedian, false);
  check_close("setup median y group 2", sm.rpm[2 * 2 + 1], 10.0, 0.0);

  // duplicate_column_map: first identical column (0-based), NaN never matches.
  const double nan = std::nan("");
  const std::vector<double> D = {1, 2, 3, 4, 1, 2, nan, 1, nan, 1, 2, 3};
  const std::vector<std::size_t> rep = nns::duplicate_column_map(D.data(), 3, 4);
  check(rep == std::vector<std::size_t>({0, 1, 2, 0}), "duplicate column map");

  // xstar_path: running weighted mean over the ordered columns.
  const std::vector<double> tr = {1, 2, 10, 20}, te = {3, 30};  // 2x2 and 1x2
  const double coef[2] = {2.0, 0.0};
  const std::size_t order[2] = {0, 1};
  const nns::XstarPathResult xs = nns::xstar_path(tr.data(), 2, te.data(), 1, 2, coef, order);
  check(xs.denominator == std::vector<double>({1, 1}), "xstar denominators");
  check(xs.train == std::vector<double>({2, 4, 2, 4}), "xstar train");
  check(xs.representative == std::vector<std::size_t>({0, 0}), "xstar representative");

  return nns_test::check_summary("test_mreg");
}
