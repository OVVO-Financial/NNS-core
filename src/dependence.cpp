// src/dependence.cpp
//
// Implementation extracted from NNS NNS_dep.cpp (synced to R 13.3, 79c60272),
// plus the NNS.dep / NNS.dep.matrix drivers that build the partitions with
// nns::partition as R's native path does. Decoupled from Rcpp.
//
// SPDX-License-Identifier: GPL-3.0-only
#include "nns/dependence.hpp"
#include "nns/parallel.hpp"
#include "nns/partition.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace nns {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kEps = std::numeric_limits<double>::epsilon();

// Rmath pnorm(q, 0, 1, lower = TRUE).
inline double pnorm_std(double q) {
  return 0.5 * std::erfc(-q / std::sqrt(2.0));
}

// Rmath dhyper(k_lo, a, m - a, b): at k_lo = max(0, a + b - m) it is one product of ratios.
double dhyper_lowest(int m, int a, int b) {
  double p = 1.0;
  if (a + b <= m) {
    for (int i = 0; i < b; ++i) p *= static_cast<double>(m - a - i) / static_cast<double>(m - i);
  } else {
    for (int i = 0; i < m - b; ++i) p *= static_cast<double>(a - i) / static_cast<double>(m - i);
  }
  return p;
}

struct CellCopula {
  double signed_value;  // debiased dependence times the sign of the cell covariance
  double magnitude;     // debiased dependence; may be negative below chance
};

// E0|2 p0 - 1| under random pairing: the co-lower count K is hypergeometric.
double null_T0(int m, int a, int b) {
  const int k_lo = std::max(0, a + b - m);
  const int k_hi = std::min(a, b);
  if (k_lo >= k_hi) return 1.0;  // a single attainable count carries no information
  const double md = static_cast<double>(m);
  if (m > 400) {
    // Normal approximation with the exact hypergeometric mean and variance.
    const double ad = a, bd = b;
    const double mean_k = ad * bd / md;
    const double var_k = bd * (ad / md) * ((md - ad) / md) * ((md - bd) / (md - 1.0));
    const double mu = (2.0 * (md - ad - bd + 2.0 * mean_k) - md) / md;
    const double sd = 4.0 * std::sqrt(std::max(var_k, 0.0)) / md;
    if (sd <= 0.0) return std::abs(mu);
    const double z = mu / sd;
    return sd * std::sqrt(2.0 / kPi) * std::exp(-0.5 * z * z) + mu * (1.0 - 2.0 * pnorm_std(-z));
  }
  double pk = dhyper_lowest(m, a, b);
  double e = 0.0;
  for (int k = k_lo; k <= k_hi; ++k) {
    const double conc = static_cast<double>(m - a - b + 2 * k);
    e += pk * std::abs(2.0 * conc / md - 1.0);
    if (k < k_hi) {
      pk *= (static_cast<double>(a - k) * static_cast<double>(b - k)) /
            (static_cast<double>(k + 1) * static_cast<double>(m - a - b + k + 1));
    }
  }
  return e;
}

// log(k!) without std::lgamma, which writes the global signgam and so races
// across dep_matrix workers.
double log_factorial(int k) {
  if (k < 2) return 0.0;
  if (k <= 170) {
    double f = 1.0;
    for (int i = 2; i <= k; ++i) f *= static_cast<double>(i);
    return std::log(f);
  }
  const double x = static_cast<double>(k), x2 = x * x;
  return x * std::log(x) - x + 0.5 * std::log(2.0 * kPi * x) +
         (1.0 / 12.0 - (1.0 / 360.0 - 1.0 / (1260.0 * x2)) / x2) / x;
}

// log(m! / prod(multiplicity!)), the number of distinct arrangements of v.
double log_arrangements(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  const int m = static_cast<int>(v.size());
  double out = log_factorial(m);
  for (int i = 0; i < m;) {
    int j = i;
    while (j < m && v[j] == v[i]) ++j;
    out -= log_factorial(j - i);
    i = j;
  }
  return out;
}

// Up to 6! = 720 distinct arrangements are enumerated exactly.
const double kLogExactArrangements = std::log(720.0) + 1e-9;
constexpr int kExactCellSize = 6;

// E0 of |sum dx dy| / sum |dx dy| under random pairing; legacy is the rule dep_local was tuned on.
double null_T1(const std::vector<double>& dx, const std::vector<double>& dy,
               bool x_ties, bool y_ties, bool legacy = false) {
  const int m = static_cast<int>(dx.size());
  // Permute whichever side has fewer distinct arrangements; the null is the same.
  double lx = 1e300, ly = 1e300;
  if (m <= kExactCellSize || (x_ties && !legacy)) lx = log_arrangements(dx);
  if (m <= kExactCellSize || (y_ties && !legacy)) ly = log_arrangements(dy);
  if (std::min(lx, ly) <= kLogExactArrangements) {
    std::vector<double> v = (lx <= ly) ? dx : dy;
    const std::vector<double>& w = (lx <= ly) ? dy : dx;
    std::sort(v.begin(), v.end());
    double total = 0.0;
    double count = 0.0;
    do {
      double s = 0.0, d = 0.0;
      for (int i = 0; i < m; ++i) {
        const double p = w[i] * v[i];
        s += p;
        d += std::abs(p);
      }
      total += d > 0.0 ? std::abs(s) / d : 0.0;
      count += 1.0;
    } while (std::next_permutation(v.begin(), v.end()));
    return total / count;
  }
  double sx2 = 0.0, sy2 = 0.0, sax = 0.0, say = 0.0;
  for (int i = 0; i < m; ++i) {
    sx2 += dx[i] * dx[i];
    sy2 += dy[i] * dy[i];
    sax += std::abs(dx[i]);
    say += std::abs(dy[i]);
  }
  if (sax <= 0.0 || say <= 0.0) return 1.0;
  // sum dx dy has permutation variance sx2 sy2 / (m - 1); sum |dx dy| has mean sax say / m.
  const double md = static_cast<double>(m);
  const double e = std::sqrt(2.0 / kPi) * std::sqrt(sx2 * sy2 / (md - 1.0)) * md / (sax * say);
  return std::min(1.0, e);
}

inline double debias_term(double t, double e0) {
  if (!(e0 < 1.0 - 1e-12)) return 0.0;
  return (t - e0) / (1.0 - e0);
}

enum CellKind { CELL_VARYING = 0, CELL_X_CONST = 1, CELL_BOTH_CONST = 2 };

// Degree-0 and degree-1 concordance, each debiased against its within-cell permutation null.
CellCopula cell_dependence(const std::vector<double>& xv, const std::vector<double>& yv,
                           bool x_ties, bool y_ties, CellKind& kind, bool legacy_null = false) {
  const int n = static_cast<int>(xv.size());
  kind = CELL_VARYING;
  if (n < 2) return {0.0, 0.0};

  // Measure from the cell minimum, so deviations round at the cell's spread rather than its level.
  double rx = 0.0, ry = 0.0;
  if (!legacy_null) {
    rx = *std::min_element(xv.begin(), xv.end());
    ry = *std::min_element(yv.begin(), yv.end());
  }
  double tx = 0.0, ty = 0.0, ax = 0.0, ay = 0.0;
  for (int i = 0; i < n; ++i) {
    tx += xv[i] - rx;
    ty += yv[i] - ry;
    ax = std::max(ax, xv[i] - rx);
    ay = std::max(ay, yv[i] - ry);
  }
  tx /= n;
  ty /= n;
  // The mean's rounding depends on row order, so a value equal to the mean can't fall either side.
  const double tie_x = legacy_null ? 0.0 : (n + 1) * kEps * ax;
  const double tie_y = legacy_null ? 0.0 : (n + 1) * kEps * ay;

  std::vector<double> dx(n), dy(n);
  int a = 0, b = 0;
  double conc = 0.0, s = 0.0, d = 0.0;
  for (int i = 0; i < n; ++i) {
    dx[i] = (xv[i] - rx) - tx;
    dy[i] = (yv[i] - ry) - ty;
    const bool lo_x = dx[i] <= tie_x, lo_y = dy[i] <= tie_y;
    a += lo_x;
    b += lo_y;
    // Ties go to the lower side so the degree-0 counts partition the cell.
    if (lo_x == lo_y) conc += 1.0;
    const double p = dx[i] * dy[i];
    s += p;
    d += std::abs(p);
  }

  // Compare raw values: a repeated 0.1 need not centre to exactly 0.
  bool x_const = true, y_const = true;
  for (int i = 1; i < n; ++i) {
    if (xv[i] != xv[0]) x_const = false;
    if (yv[i] != yv[0]) y_const = false;
  }
  if (x_const) {
    kind = y_const ? CELL_BOTH_CONST : CELL_X_CONST;
    return {0.0, 0.0};
  }

  // s is summed in partition order, so a zero covariance can round either way.
  const double s_tol = legacy_null ? 0.0 : 1e-12 * d;
  const double slope_sign = (s > s_tol) ? 1.0 : (s < -s_tol) ? -1.0 : 0.0;
  const double e0 = null_T0(n, a, b);

  // One attainable count is uninformative; compare counts, as 49 * (1/49) != 1.
  if (e0 < 1.0 - 1e-12) {
    if (conc == static_cast<double>(n)) return {1.0, 1.0};
    if (conc == 0.0) return {-1.0, 1.0};
  }

  const double md = static_cast<double>(n);
  const double t0 = std::abs(2.0 * conc / md - 1.0);
  const double t1 = d > 0.0 ? std::abs(s) / d : 0.0;
  const double e1 = null_T1(dx, dy, x_ties, y_ties, legacy_null);

  const double dep = 0.5 * (debias_term(t0, e0) + debias_term(t1, e1));
  if (slope_sign == 0.0) return {0.0, 0.0};
  return {dep * slope_sign, dep};
}

struct LevelBest {
  double dependence;
  double correlation;
};

enum PartitionKind { CONDITIONAL, JOINT };

// Chance that y is constant in a cell of m points under the permutation null.
double p_constant(int m, const std::vector<int>& counts, int n, std::vector<double>& cache) {
  if (cache[m] >= 0.0) return cache[m];
  const double lchoose_n = log_factorial(n) - log_factorial(m) - log_factorial(n - m);
  double p = 0.0;
  for (int c : counts) {
    if (c < m) break;
    p += std::exp(log_factorial(c) - log_factorial(m) - log_factorial(c - m) - lchoose_n);
  }
  cache[m] = std::min(1.0, p);
  return cache[m];
}

std::vector<int> value_counts(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  std::vector<int> counts;
  for (std::size_t i = 0; i < v.size();) {
    std::size_t j = i;
    while (j < v.size() && v[j] == v[i]) ++j;
    counts.push_back(static_cast<int>(j - i));
    i = j;
  }
  std::sort(counts.begin(), counts.end(), std::greater<int>());
  return counts;
}

// The first len characters of a and b agree (a shorter label counts whole).
inline bool same_prefix(const std::string& a, const std::string& b, std::size_t len) {
  const std::size_t la = std::min(a.size(), len), lb = std::min(b.size(), len);
  return la == lb && a.compare(0, la, b, 0, lb) == 0;
}

// Strongest level of a prefix-coded partition ("q", "q1", "q12", ...); ties go to the finer level.
LevelBest best_level(const std::vector<double>& u, const std::vector<double>& v,
                     const std::vector<std::string>& lab, double between_sign,
                     PartitionKind kind, bool u_ties, bool v_ties,
                     std::vector<double>* levels = nullptr) {
  const int n = static_cast<int>(u.size());
  std::vector<int> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  std::stable_sort(ord.begin(), ord.end(), [&](int i, int j) { return lab[i] < lab[j]; });
  std::size_t L = 1;
  for (const auto& s : lab) L = std::max(L, s.size());

  std::vector<int> v_counts;
  std::vector<double> p0_cache;
  if (kind == CONDITIONAL && v_ties) {
    v_counts = value_counts(v);
    p0_cache.assign(n + 1, -1.0);
  }

  LevelBest best{-1e300, 0.0};
  int prev_cells = -1;
  std::vector<double> uq, vq;
  std::vector<int> starts;
  for (std::size_t len = 1; len <= L; ++len) {
    starts.clear();
    starts.push_back(0);
    for (int k = 1; k < n; ++k) {
      if (!same_prefix(lab[ord[k]], lab[ord[k - 1]], len)) starts.push_back(k);
    }
    const int cells = static_cast<int>(starts.size());
    if (cells == prev_cells) continue;  // nested partition unchanged: same value
    prev_cells = cells;
    starts.push_back(n);

    double dep = 0.0, corr = 0.0;
    for (int c = 0; c < cells; ++c) {
      const int lo = starts[c], hi = starts[c + 1], m = hi - lo;
      uq.resize(m);
      vq.resize(m);
      for (int k = 0; k < m; ++k) { uq[k] = u[ord[lo + k]]; vq[k] = v[ord[lo + k]]; }
      CellKind ck;
      CellCopula cop = cell_dependence(uq, vq, u_ties, v_ties, ck);
      if (ck != CELL_VARYING) {
        double mag = 0.0;
        if (kind == CONDITIONAL) {
          const double p0 = v_ties ? p_constant(m, v_counts, n, p0_cache) : 0.0;
          const double score = (ck == CELL_BOTH_CONST) ? 1.0 : 0.0;
          if (p0 < 1.0 - 1e-12) mag = (score - p0) / (1.0 - p0);
        }
        // Conditional: is v determined here, against chance; joint: constant by construction.
        cop = {between_sign * mag, mag};  // no within-cell sign
      }
      if (!std::isfinite(cop.signed_value) || !std::isfinite(cop.magnitude)) cop = {0.0, 0.0};
      const double w = static_cast<double>(m) / static_cast<double>(n);
      dep += w * cop.magnitude;
      corr += w * cop.signed_value;
    }
    if (levels) levels->push_back(dep);
    if (dep >= best.dependence) best = {dep, corr};
  }
  // Dependence is a magnitude: floored at 0, with |Correlation| kept within it.
  best.dependence = std::max(0.0, best.dependence);
  best.correlation = std::min(best.dependence, std::max(-best.dependence, best.correlation));
  return best;
}

bool has_ties(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  return std::adjacent_find(v.begin(), v.end()) != v.end();
}

int count_unique(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  return static_cast<int>(std::unique(v.begin(), v.end()) - v.begin());
}

bool either_constant(const std::vector<double>& xv, const std::vector<double>& yv) {
  bool cx = true, cy = true;
  for (std::size_t i = 1; i < xv.size(); ++i) {
    if (xv[i] != xv[0]) cx = false;
    if (yv[i] != yv[0]) cy = false;
    if (!cx && !cy) break;
  }
  return xv.size() < 2 || cx || cy;
}

// NaN never equals itself, so the run scans over tied values would not advance.
void check_finite(const double* x, const double* y, std::size_t n) {
  for (std::size_t i = 0; i < n; ++i) {
    if (!std::isfinite(x[i])) throw std::invalid_argument("x must be finite.");
    if (!std::isfinite(y[i])) throw std::invalid_argument("y must be finite.");
  }
}

void check_labels(const std::vector<std::string>& q, std::size_t n, const char* msg) {
  if (q.size() != n) throw std::invalid_argument(msg);
}

double copula_degree0_unsigned(const std::vector<double>& xv, const std::vector<double>& yv) {
  const int n = static_cast<int>(xv.size());
  if (n < 2) return 0.0;
  double tx = 0.0, ty = 0.0;
  for (int i = 0; i < n; ++i) { tx += xv[i]; ty += yv[i]; }
  tx /= n;
  ty /= n;
  double conc = 0.0, disc = 0.0;
  for (int i = 0; i < n; ++i) {
    const bool lo_x = (xv[i] - tx) <= 0.0, lo_y = (yv[i] - ty) <= 0.0;
    if (lo_x == lo_y) conc += 1.0; else disc += 1.0;
  }
  const double d0_Co = conc / n, dpm_d0 = disc / n;
  const double disc_dep = std::min(1.0, std::max(0.0, std::abs(d0_Co - 0.5) / 0.5));
  const double nd_disc = std::abs(dpm_d0 - 0.5) / 0.5;
  return std::sqrt((disc_dep + nd_disc) / 2.0);
}

// Finest cells only, each as a signed square root: a fine-scale signal-to-noise measure.
double local_dep_dir(const std::vector<double>& u, const std::vector<double>& v,
                     const std::vector<std::string>& lab, bool u_ties, bool v_ties) {
  const int n = static_cast<int>(u.size());
  auto root = [](double inner) { return (inner >= 0.0 ? 1.0 : -1.0) * std::sqrt(std::abs(inner)); };
  auto cell_value = [&](const std::vector<double>& uq, const std::vector<double>& vq) -> double {
    CellKind ck;
    CellCopula c = cell_dependence(uq, vq, u_ties, v_ties, ck, true);
    if (ck == CELL_BOTH_CONST) return 1.0;
    if (ck == CELL_X_CONST) return 0.0;
    return root(c.magnitude);
  };
  double global = cell_value(u, v);
  if (!std::isfinite(global)) global = 0.0;

  std::vector<int> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  std::stable_sort(ord.begin(), ord.end(), [&](int i, int j) { return lab[i] < lab[j]; });
  double dep = 0.0;
  std::vector<double> uq, vq;
  for (int lo = 0; lo < n;) {
    int hi = lo + 1;
    while (hi < n && lab[ord[hi]] == lab[ord[lo]]) ++hi;
    const int m = hi - lo;
    uq.resize(m);
    vq.resize(m);
    for (int k = 0; k < m; ++k) { uq[k] = u[ord[lo + k]]; vq[k] = v[ord[lo + k]]; }
    double val = cell_value(uq, vq);
    if (!std::isfinite(val)) val = global;
    dep += val * static_cast<double>(m) / static_cast<double>(n);
    lo = hi;
  }
  return std::max(0.0, dep);
}

// NNS.dep's partitions: the mean split on `a` alone (XONLY) or on both.
std::vector<std::string> dep_partition(const std::vector<double>& a, const std::vector<double>& b,
                                       bool xonly) {
  const std::size_t n = a.size();
  const int order = std::max(static_cast<int>(std::ceil(std::log2(static_cast<double>(n)))), 1);
  const int obs = std::max(8, static_cast<int>(static_cast<double>(n) / 8.0));
  const std::optional<std::string> type =
      xonly ? std::optional<std::string>("XONLY") : std::nullopt;
  return partition(a.data(), b.data(), n, type, order, obs, false, "off", true).quadrant;
}

} // namespace

// ---------- Pairwise Dependence Kernel ----------

DepResult dep_pair(const double* x, const double* y, std::size_t n,
                   const std::vector<std::string>& quad_xy,
                   const std::vector<std::string>& quad_yx,
                   const std::vector<std::string>& quad_joint, bool asym) {
  check_finite(x, y, n);
  check_labels(quad_xy, n, "x, y and quad_xy must all have the same length");
  if (!asym && (quad_yx.size() != n || quad_joint.size() != n))
    throw std::invalid_argument("quad_yx and quad_joint must have the length of x when asym = FALSE");

  const std::vector<double> xv(x, x + n), yv(y, y + n);
  if (either_constant(xv, yv)) return {0.0, 0.0};

  // Centred on the first observation so a large common level doesn't round the sign away.
  const double x0 = xv[0], y0 = yv[0];
  double mx = 0.0, my = 0.0, cov = 0.0, cov_abs = 0.0;
  for (std::size_t i = 0; i < n; ++i) { mx += xv[i] - x0; my += yv[i] - y0; }
  mx /= static_cast<double>(n);
  my /= static_cast<double>(n);
  for (std::size_t i = 0; i < n; ++i) {
    const double p = ((xv[i] - x0) - mx) * ((yv[i] - y0) - my);
    cov += p;
    cov_abs += std::abs(p);
  }
  const double cov_tol = 1e-12 * cov_abs;
  const double between_sign = (cov > cov_tol) ? 1.0 : (cov < -cov_tol) ? -1.0 : 0.0;

  const bool x_ties = has_ties(xv), y_ties = has_ties(yv);

  LevelBest best = best_level(xv, yv, quad_xy, between_sign, CONDITIONAL, x_ties, y_ties);
  if (!asym) {
    // The joint quadrants catch a circle or y = x * noise; they win ties.
    LevelBest out = best_level(xv, yv, quad_joint, between_sign, JOINT, x_ties, y_ties);
    LevelBest on_y = best_level(yv, xv, quad_yx, between_sign, CONDITIONAL, y_ties, x_ties);
    if (best.dependence > out.dependence) out = best;
    if (on_y.dependence > out.dependence) out = on_y;
    best = out;
  }
  return {best.correlation, best.dependence};
}

std::vector<double> dep_levels(const double* x, const double* y, std::size_t n,
                               const std::vector<std::string>& quad_xy) {
  check_finite(x, y, n);
  check_labels(quad_xy, n, "x, y and quad_xy must all have the same length");
  const std::vector<double> xv(x, x + n), yv(y, y + n);
  if (either_constant(xv, yv)) return {0.0};
  std::vector<double> levels;
  best_level(xv, yv, quad_xy, 0.0, CONDITIONAL, has_ties(xv), has_ties(yv), &levels);
  return levels;
}

double dep_local(const double* x, const double* y, std::size_t n,
                 const std::vector<std::string>& quad_xy,
                 const std::vector<std::string>& quad_yx, bool asym) {
  check_finite(x, y, n);
  check_labels(quad_xy, n, "x, y and quad_xy must all have the same length");
  if (!asym) check_labels(quad_yx, n, "quad_yx must have the length of x when asym = FALSE");
  const std::vector<double> xv(x, x + n), yv(y, y + n);
  if (either_constant(xv, yv)) return 0.0;
  const bool x_ties = has_ties(xv), y_ties = has_ties(yv);
  double dep = local_dep_dir(xv, yv, quad_xy, x_ties, y_ties);
  if (!asym) dep = std::max(dep, local_dep_dir(yv, xv, quad_yx, y_ties, x_ties));
  const double rn = std::sqrt(static_cast<double>(n));
  if (count_unique(xv) < rn && count_unique(yv) < rn) {
    double disc = copula_degree0_unsigned(xv, yv);
    if (!std::isfinite(disc)) disc = dep;
    dep = 0.5 * (dep + disc);
  }
  return dep;
}

// ---------- NNS.dep / NNS.dep.matrix ----------

DepResult dep(const double* x, const double* y, std::size_t n, bool asym) {
  check_finite(x, y, n);
  if (n < 2) return {0.0, 0.0};
  const std::vector<double> xv(x, x + n), yv(y, y + n);
  const std::vector<std::string> qxy = dep_partition(xv, yv, true);
  if (asym) return dep_pair(x, y, n, qxy, {}, {}, true);
  return dep_pair(x, y, n, qxy, dep_partition(yv, xv, true), dep_partition(xv, yv, false), false);
}

DepMatrixResult dep_matrix(const double* X, std::size_t n, std::size_t p,
                           bool asym, int nthreads) {
  if (p < 2) throw std::invalid_argument("dep_matrix: X must have at least 2 columns");

  const std::size_t n_pairs = p * (p - 1) / 2;
  std::vector<std::size_t> pair_i, pair_j;
  pair_i.reserve(n_pairs);
  pair_j.reserve(n_pairs);
  for (std::size_t i = 0; i + 1 < p; ++i) {
    for (std::size_t j = i + 1; j < p; ++j) {
      pair_i.push_back(i);
      pair_j.push_back(j);
    }
  }

  std::vector<DepResult> upper(n_pairs), lower(n_pairs);
  parallel_for(0, n_pairs, [&](std::size_t begin, std::size_t end) {
    for (std::size_t k = begin; k < end; ++k) {
      const double* xi = X + pair_i[k] * n;
      const double* xj = X + pair_j[k] * n;
      upper[k] = dep(xi, xj, n, asym);
      lower[k] = asym ? dep(xj, xi, n, true) : upper[k];
    }
  }, nthreads, n * 64U);  // three partitions and several passes per pair

  DepMatrixResult result;
  result.p = p;
  result.correlation.assign(p * p, 0.0);
  result.dependence.assign(p * p, 0.0);
  for (std::size_t i = 0; i < p; ++i) {
    result.correlation[i * p + i] = 1.0;
    result.dependence[i * p + i] = 1.0;
  }
  for (std::size_t k = 0; k < n_pairs; ++k) {
    const std::size_t i = pair_i[k], j = pair_j[k];
    result.correlation[j * p + i] = upper[k].correlation;
    result.dependence[j * p + i] = upper[k].dependence;
    result.correlation[i * p + j] = lower[k].correlation;
    result.dependence[i * p + j] = lower[k].dependence;
  }
  return result;
}

} // namespace nns
