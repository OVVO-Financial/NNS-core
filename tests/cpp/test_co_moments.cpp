// tests/cpp/test_co_moments.cpp
//
// Bivariate co- and divergent partial moments against a direct mean of the
// product of their legs, for integer, fractional and zero degrees.
//
// SPDX-License-Identifier: GPL-3.0-only
#include "nns/partial_moments.hpp"

#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include "check.hpp"

using nns_test::check_close;
using nns_test::check_throws;

namespace {

double lower_leg(double v, double t, double d) {
  if (d == 0) return v <= t ? 1.0 : 0.0;
  return t > v ? std::pow(t - v, d) : 0.0;
}

double upper_leg(double v, double t, double d) {
  if (d == 0) return v > t ? 1.0 : 0.0;
  return v > t ? std::pow(v - t, d) : 0.0;
}

template <typename Leg1, typename Leg2>
double direct(const std::vector<double>& x, const std::vector<double>& y, Leg1 lx, Leg2 ly) {
  double s = 0.0;
  for (std::size_t i = 0; i < x.size(); ++i) s += lx(x[i]) * ly(y[i]);
  return s / static_cast<double>(x.size());
}

}  // namespace

int main() {
  // D.UPM(.5, 1, ...) took the .5 to y and truncated it; R NNS 13.3 returned 2.5.
  const std::vector<double> a = {-4.0, -9.0}, b = {2.0, 3.0};
  check_close("d_upm lpm .5 upm 1", nns::d_upm(0.5, 1.0, a.data(), b.data(), 2, 2, 0.0, 0.0), 6.5,
              1e-15);
  const std::vector<double> c = {-2.0, -3.0}, e = {4.0, 9.0};
  check_close("d_upm lpm 1 upm .5", nns::d_upm(1.0, 0.5, c.data(), e.data(), 2, 2, 0.0, 0.0), 6.5,
              1e-15);

  std::vector<double> x(60), y(60);
  for (std::size_t i = 0; i < x.size(); ++i) {
    x[i] = 2.0 * std::sin(1.7 * static_cast<double>(i) + 0.4);
    y[i] = 1.5 * std::cos(0.9 * static_cast<double>(i) + 0.3);
  }
  const double degrees[] = {0.0, 0.5, 1.0, 1.5, 2.0, 3.0};
  const std::size_t n = x.size();
  for (double tx : {0.1, -0.3}) {
    for (double ty : {-0.2, 0.4}) {
      for (double p : degrees) {
        for (double q : degrees) {
          const std::string at = " " + std::to_string(p) + "," + std::to_string(q) + " at " +
                                 std::to_string(tx) + "," + std::to_string(ty);
          auto lo = [](double t, double d) { return [=](double v) { return lower_leg(v, t, d); }; };
          auto up = [](double t, double d) { return [=](double v) { return upper_leg(v, t, d); }; };
          check_close("co_lpm" + at, nns::co_lpm(p, q, x.data(), y.data(), n, n, tx, ty),
                      direct(x, y, lo(tx, p), lo(ty, q)), 1e-12, 1e-15);
          check_close("co_upm" + at, nns::co_upm(p, q, x.data(), y.data(), n, n, tx, ty),
                      direct(x, y, up(tx, p), up(ty, q)), 1e-12, 1e-15);
          check_close("d_lpm" + at, nns::d_lpm(p, q, x.data(), y.data(), n, n, tx, ty),
                      direct(x, y, up(tx, q), lo(ty, p)), 1e-12, 1e-15);
          check_close("d_upm" + at, nns::d_upm(p, q, x.data(), y.data(), n, n, tx, ty),
                      direct(x, y, lo(tx, p), up(ty, q)), 1e-12, 1e-15);
        }
      }
    }
  }

  std::vector<double> targets = {-1.0, 0.0, 1.0}, out(3);
  const double ty = 0.5;
  nns::d_upm_v(0.5, 2.0, x.data(), y.data(), n, n, targets.data(), 3, &ty, 1, out.data());
  for (std::size_t i = 0; i < 3; ++i)
    check_close("d_upm_v " + std::to_string(i), out[i],
                nns::d_upm(0.5, 2.0, x.data(), y.data(), n, n, targets[i], ty), 0.0);

  // Unequal lengths summed min(n_x, n_y) pairs over max(n_x, n_y); an empty
  // target made the workers take i % 0.
  using Scalar = double (*)(double, double, const double*, const double*, std::size_t,
                            std::size_t, double, double);
  using Vector = void (*)(double, double, const double*, const double*, std::size_t, std::size_t,
                          const double*, std::size_t, const double*, std::size_t, double*, int);
  const std::pair<const char*, Scalar> scalars[] = {
      {"co_lpm", nns::co_lpm}, {"co_upm", nns::co_upm}, {"d_lpm", nns::d_lpm}, {"d_upm", nns::d_upm}};
  const std::pair<const char*, Vector> vectors[] = {{"co_lpm_v", nns::co_lpm_v},
                                                    {"co_upm_v", nns::co_upm_v},
                                                    {"d_lpm_v", nns::d_lpm_v},
                                                    {"d_upm_v", nns::d_upm_v}};
  const std::vector<double> t200(200, 0.0);
  std::vector<double> buf(200);
  for (const auto& entry : scalars) {
    const std::string s = entry.first;
    const auto fn = entry.second;
    check_throws(s + " n_x != n_y", [&] { fn(1, 1, x.data(), y.data(), n, n - 1, 0, 0); });
    check_throws(s + " empty", [&] { fn(1, 1, x.data(), y.data(), 0, 0, 0, 0); });
  }
  for (const auto& entry : vectors) {
    const std::string s = entry.first;
    const auto fn = entry.second;
    const double* t = t200.data();
    double* o = buf.data();
    check_throws(s + " n_x != n_y", [&] { fn(1, 1, x.data(), y.data(), n, n - 1, t, 1, t, 1, o, -1); });
    check_throws(s + " n_x != n_y, 200 targets",
                 [&] { fn(1, 1, x.data(), y.data(), n, n - 1, t, 200, t, 1, o, -1); });
    check_throws(s + " empty x, y", [&] { fn(1, 1, x.data(), y.data(), 0, 0, t, 1, t, 1, o, -1); });
    check_throws(s + " empty target_x", [&] { fn(1, 1, x.data(), y.data(), n, n, t, 0, t, 1, o, -1); });
    check_throws(s + " empty target_y", [&] { fn(1, 1, x.data(), y.data(), n, n, t, 1, t, 0, o, -1); });
    check_throws(s + " empty targets", [&] { fn(1, 1, x.data(), y.data(), n, n, t, 0, t, 0, o, -1); });
    check_throws(s + " 3 and 2 targets", [&] { fn(1, 1, x.data(), y.data(), n, n, t, 3, t, 2, o, -1); });
  }
  return nns_test::check_summary("test_co_moments");
}
