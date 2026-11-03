// src/central_tendencies.cpp
//
// Implementation reconstructed from original_src/central_tendencies.cpp; covers NNS_gravity_cpp, NNS_rescale_cpp, and NNS_mode_cpp. Decoupled from Rcpp.
//
// SPDX-License-Identifier: GPL-3.0-only
#include "nns/central_tendencies.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_map>

namespace nns {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// ---------- helpers ----------

inline double frac_part(double x) {
  return x - std::floor(x);
}

inline double mean_vec(const std::vector<double>& v) {
  if (v.empty()) return kNaN;
  long double s = 0.0L;
  for (double x : v) s += x;
  return static_cast<double>(s / v.size());
}

inline double nearest_int_half_up(double x) {
  double f = std::floor(x);
  return ((x - f) < 0.5) ? f : std::ceil(x);
}

// Given a sorted vector xs, reproduce the q1, q2, q3 *exactly* as in the R code.
void quartiles_like_R_code(const std::vector<double>& xs, double& q1, double& q2, double& q3) {
  const int l = static_cast<int>(xs.size());
  const double l25 = l * 0.25;
  const double l50 = l * 0.50;
  const double l75 = l * 0.75;
  
  if (l % 2 == 0) {
    int i25 = std::max(1, static_cast<int>(std::floor(l25))) - 1;
    int i50 = std::max(1, static_cast<int>(std::floor(l50))) - 1;
    int i75 = std::max(1, static_cast<int>(std::floor(l75))) - 1;
    q1 = xs[i25];
    q2 = xs[i50];
    q3 = xs[i75];
  } else {
    int f25 = static_cast<int>(std::floor(l25));
    int c25 = static_cast<int>(std::ceil(l25));
    f25 = std::min(std::max(1, f25), l);
    c25 = std::min(std::max(1, c25), l);
    double w25 = frac_part(l25);
    q1 = xs[f25 - 1] + w25 * (xs[c25 - 1] - xs[f25 - 1]);
    
    int f50 = static_cast<int>(std::floor(l50));
    int c50 = static_cast<int>(std::ceil(l50));
    f50 = std::min(std::max(1, f50), l);
    c50 = std::min(std::max(1, c50), l);
    q2 = 0.5 * (xs[f50 - 1] + xs[c50 - 1]);
    
    int f75 = static_cast<int>(std::floor(l75));
    int c75 = static_cast<int>(std::ceil(l75));
    f75 = std::min(std::max(1, f75), l);
    c75 = std::min(std::max(1, c75), l);
    double w75 = frac_part(l75);
    q3 = xs[f75 - 1] + w75 * (xs[c75 - 1] - xs[f75 - 1]);
  }
}

// R's default quantile (type 7) of sorted xs: position 1 + (n - 1) p.
inline double quantile_type7(const std::vector<double>& xs, double p) {
  const int n = static_cast<int>(xs.size());
  const double h = (n - 1) * p;
  const int lo = static_cast<int>(std::floor(h));
  const int hi = std::min(lo + 1, n - 1);
  return xs[lo] + (h - lo) * (xs[hi] - xs[lo]);
}

// Minimal replacement for NNS_bin used by mode/gravity. Bins are
// [origin + k*width, origin + (k+1)*width), k = 0 .. nbins-1, with
// nbins = floor((max - origin) / width + 1e-12) + 1, as in NNS 13.3.
//
// Gravity and mode only read counts within 6 bins of a non-empty bin, so runs
// of more than kMaxEmptyRun empty bins keep only their first and last
// kMaxEmptyRun / 2 bins: results are bit-identical and memory is O(n). A width
// so small that the bin count exceeds 2^52 is widened to (max - origin) / 2^52.
constexpr std::int64_t kMaxEmptyRun = 16;
constexpr double kMaxBinIndex = 4503599627370496.0;  // 2^52

void simple_bin_counts(const std::vector<double>& xs, double width, double origin,
                       std::vector<double>& bin_names, std::vector<int>& counts) {
  bin_names.clear();
  counts.clear();
  if (xs.empty()) return;

  const double xmax = xs.back();
  double last = std::floor((xmax - origin) / width + 1e-12);
  if (!(last < kMaxBinIndex)) {
    width = (xmax - origin) / kMaxBinIndex;
    last = std::floor((xmax - origin) / width + 1e-12);
    if (!(last < kMaxBinIndex)) last = kMaxBinIndex;
  }
  const std::int64_t nbins = (last >= 0.0) ? static_cast<std::int64_t>(last) + 1 : 1;
  const std::int64_t cap = static_cast<std::int64_t>(xs.size()) * (kMaxEmptyRun + 1) + kMaxEmptyRun;
  bin_names.reserve(static_cast<std::size_t>(std::min(nbins, cap)));
  counts.reserve(static_cast<std::size_t>(std::min(nbins, cap)));

  auto push_bin = [&](std::int64_t k, int c) {
    bin_names.push_back(origin + static_cast<double>(k) * width);
    counts.push_back(c);
  };
  auto push_empty_run = [&](std::int64_t first, std::int64_t last_bin) {
    const std::int64_t len = last_bin - first + 1;
    if (len <= 0) return;
    if (len <= kMaxEmptyRun) {
      for (std::int64_t k = first; k <= last_bin; ++k) push_bin(k, 0);
      return;
    }
    const std::int64_t half = kMaxEmptyRun / 2;
    for (std::int64_t k = first; k < first + half; ++k) push_bin(k, 0);
    for (std::int64_t k = last_bin - half + 1; k <= last_bin; ++k) push_bin(k, 0);
  };

  // xs is sorted, so bin indices are non-decreasing.
  const double top = static_cast<double>(nbins - 1);
  std::int64_t prev = -1;
  for (double v : xs) {
    const double q = std::floor((v - origin) / width);
    const std::int64_t k = (q < 0.0) ? 0 : (q >= top ? nbins - 1 : static_cast<std::int64_t>(q));
    if (k == prev) {
      ++counts.back();
      continue;
    }
    push_empty_run(prev + 1, k - 1);
    push_bin(k, 1);
    prev = k;
  }
  push_empty_run(prev + 1, nbins - 1);
}

// Triangular smoothing helper: 7-tap [1,2,3,4,3,2,1] with mirrored edges
void smooth_counts_tri7(const std::vector<int>& counts, std::vector<double>& smooth) {
  static const int w[7] = {1, 2, 3, 4, 3, 2, 1};
  static const int Wsum = 16;
  const int n = static_cast<int>(counts.size());
  smooth.assign(n, 0.0);
  if (n == 0) return;
  
  // Mirror at the edges. The clamp only binds for n <= 3 bins, whose smoothed
  // values are unused (the peak scan needs n >= 7).
  auto at = [&](int idx) -> int {
    if (idx < 0) idx = -idx;
    if (idx >= n) idx = 2 * n - 2 - idx;
    idx = std::min(std::max(idx, 0), n - 1);
    return counts[idx];
  };
  
  for (int i = 0; i < n; ++i) {
    int acc = 0;
    acc += w[0] * at(i - 3); acc += w[1] * at(i - 2); acc += w[2] * at(i - 1);
    acc += w[3] * at(i);
    acc += w[4] * at(i + 1); acc += w[5] * at(i + 2); acc += w[6] * at(i + 3);
    smooth[i] = static_cast<double>(acc) / static_cast<double>(Wsum);
  }
}

} // namespace

// ---------- NNS.gravity ----------

double gravity(const double* x_in, std::size_t n, bool discrete) {
  std::vector<double> x;
  x.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    if (std::isfinite(x_in[i])) x.push_back(x_in[i]);
  }
  
  const int l = static_cast<int>(x.size());
  if (l == 0) return kNaN;
  if (l <= 3) {
    std::vector<double> t = x;
    std::sort(t.begin(), t.end());
    double med = (l % 2) ? t[l / 2] : 0.5 * (t[l / 2 - 1] + t[l / 2]);
    if (discrete) return nearest_int_half_up(med);
    return med;
  }
  
  bool all_eq = true;
  for (int i = 1; i < l; ++i) {
    if (x[i] != x[0]) { all_eq = false; break; }
  }
  if (all_eq) return x[0];
  
  std::sort(x.begin(), x.end());
  double range = std::fabs(x.back() - x.front());
  if (range == 0.0) return x.front();
  
  // The floor-position median and quartiles put gravity below the centre of symmetric data.
  const double q1 = quantile_type7(x, 0.25);
  const double q2 = quantile_type7(x, 0.50);
  const double q3 = quantile_type7(x, 0.75);
  
  double width = (q3 - q1) * std::pow(static_cast<double>(l), -0.5);
  if (!(width > 0.0) || !std::isfinite(width)) width = range / 128.0;
  
  std::vector<double> z_names;
  std::vector<int> counts;
  simple_bin_counts(x, width, x.front(), z_names, counts);
  const int lz = static_cast<int>(counts.size());
  
  int maxc = 0;
  for (int c : counts) if (c > maxc) maxc = c;
  int ties = 0;
  for (int c : counts) if (c == maxc) ++ties;
  
  int lo = 0, hi = lz - 1;
  if (ties == 1) {
    int zc = 0;
    for (int i = 0; i < lz; ++i) {
      if (counts[i] == maxc) { zc = i; break; }
    }
    lo = std::max(0, zc - 1);
    hi = std::min(lz - 1, zc + 1);
  }
  
  // Mean of the observations in bins lo..hi (a contiguous slice of sorted x): bin left edges sit half a bin low.
  int a = 0;
  for (int i = 0; i < lo; ++i) a += counts[i];
  int e = a;
  for (int i = lo; i <= hi; ++i) e += counts[i];
  long double num = 0.0L;
  for (int t = a; t < e; ++t) num += static_cast<long double>(x[t]);
  double m = (e > a) ? static_cast<double>(num / static_cast<long double>(e - a)) : z_names[(lo + hi) / 2];
  
  double mu = mean_vec(x);
  double mid = 0.25 * (q2 + m + mu + 0.5 * (q1 + q3));
  
  double out = std::isfinite(mid) ? mid : q2;
  if (discrete) out = nearest_int_half_up(out);
  return out;
}

// ---------- NNS.rescale ----------

std::vector<double> rescale(const double* x_in, std::size_t n, double a, double b,
                            const std::string& method,
                            std::optional<double> T,
                            const std::string& type) {
  std::vector<double> out(n, kNaN);
  
  std::string method_lower = method;
  std::transform(method_lower.begin(), method_lower.end(), method_lower.begin(),
                 [](unsigned char c){ return std::tolower(c); });
                 
  std::string type_lower = type;
  std::transform(type_lower.begin(), type_lower.end(), type_lower.begin(),
                 [](unsigned char c){ return std::tolower(c); });
  
  if (method_lower == "minmax") {
    double xmin = std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();
    
    for (std::size_t i = 0; i < n; ++i) {
      if (std::isfinite(x_in[i])) {
        if (x_in[i] < xmin) xmin = x_in[i];
        if (x_in[i] > xmax) xmax = x_in[i];
      }
    }
    
    // Fallback if all values identical or no valid values
    if (!std::isfinite(xmin) || !std::isfinite(xmax) || xmax == xmin) {
      for (std::size_t i = 0; i < n; ++i) out[i] = (a + b) / 2.0;
      return out;
    }
    
    for (std::size_t i = 0; i < n; ++i) {
      out[i] = a + (b - a) * ((x_in[i] - xmin) / (xmax - xmin));
    }
    return out;
  }
  
  if (method_lower == "riskneutral") {
    if (!T.has_value()) {
      throw std::invalid_argument("T (time to maturity) must be provided for riskneutral method");
    }
    double T_val = T.value();
    
    if (!(a > 0.0)) {
      throw std::invalid_argument("S_0 (a) must be positive for riskneutral method");
    }
    
    double S0 = a;
    double r = b;
    
    long double s = 0.0L; 
    int cnt = 0;
    for (std::size_t i = 0; i < n; ++i) {
      if (std::isfinite(x_in[i])) { s += x_in[i]; ++cnt; }
    }
    double mx = (cnt > 0) ? static_cast<double>(s / cnt) : kNaN;
    
    if (!std::isfinite(mx) || mx <= 0.0) {
      throw std::invalid_argument("Mean(x) must be positive/finite for riskneutral scaling");
    }
    
    double target = (type_lower == "discounted") ? S0 : (S0 * std::exp(r * T_val));
    double theta = std::log(target / mx);
    
    for (std::size_t i = 0; i < n; ++i) {
      out[i] = x_in[i] * std::exp(theta);
    }
    return out;
  }
  
  throw std::invalid_argument("Invalid method: use 'minmax' or 'riskneutral'");
}

// ---------- NNS.mode ----------

std::vector<double> mode(const double* x_in, std::size_t n, bool discrete, bool multi) {
  std::vector<double> xnum; 
  xnum.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    if (std::isfinite(x_in[i])) xnum.push_back(x_in[i]);
  }
  
  const int l = static_cast<int>(xnum.size());
  if (l == 0) return {kNaN};
  
  // ====================== DISCRETE PATH ======================
  if (discrete) {
    // Key on the rounded value itself: NNS 13.3's int key overflowed beyond the
    // int range. "+ 0.0" folds -0.0 into the 0 key.
    std::unordered_map<double, int> freq;
    freq.reserve(l * 2u);
    for (double v : xnum) ++freq[nearest_int_half_up(v) + 0.0];
    
    int maxf = 0; 
    for (const auto& kv : freq) if (kv.second > maxf) maxf = kv.second;
    
    std::vector<double> modes_int;
    for (const auto& kv : freq) if (kv.second == maxf) modes_int.push_back(kv.first);
    std::sort(modes_int.begin(), modes_int.end());
    
    if (multi) {
      return modes_int;
    } else {
      // Tied modes: nearest the gravity, then the lower value; a mean of tied classes is no class.
      if (modes_int.size() == 1) return modes_int;
      const double g = gravity(xnum.data(), xnum.size(), false);
      double best_c = kNaN, best_d = std::numeric_limits<double>::infinity();
      for (double m : modes_int) {
        const double dist = std::fabs(m - g);
        if (dist < best_d) { best_d = dist; best_c = m; }
      }
      return {best_c};
    }
  }
  
  // ====================== CONTINUOUS PATH ======================
  if (l <= 3) {
    std::vector<double> tmp = xnum; 
    std::sort(tmp.begin(), tmp.end());
    double med = (l % 2 == 1) ? tmp[l / 2] : 0.5 * (tmp[l / 2 - 1] + tmp[l / 2]);
    return {med};
  }
  
  bool all_eq = true;
  for (int i = 1; i < l; ++i) {
    if (xnum[i] != xnum[0]) { all_eq = false; break; }
  }
  if (all_eq) return {xnum[0]};
  
  std::sort(xnum.begin(), xnum.end());
  double range = std::fabs(xnum.back() - xnum.front());
  if (range == 0.0) return {xnum.front()};
  
  double q1, q2, q3;
  quartiles_like_R_code(xnum, q1, q2, q3);
  double width = (q3 - q1) * std::pow(static_cast<double>(l), -0.5);
  if (!(width > 0.0) || !std::isfinite(width)) width = range / 128.0;

  // Every gap wider than a bin: smoothing would move the peaks off the values, so use their frequencies.
  {
    std::vector<double> vals;
    std::vector<int> cnt;
    for (int i = 0; i < l;) {
      int j = i + 1;
      while (j < l && xnum[j] == xnum[i]) ++j;
      vals.push_back(xnum[i]);
      cnt.push_back(j - i);
      i = j;
    }
    const int k = static_cast<int>(vals.size());
    bool coarse = k >= 2;
    for (int i = 1; coarse && i < k; ++i) {
      if (!(vals[i] - vals[i - 1] > width)) coarse = false;
    }
    if (coarse) {
      std::vector<int> modes;
      for (int i = 0; i < k;) {
        int j = i + 1;
        while (j < k && cnt[j] == cnt[i]) ++j;
        const bool left_ok = (i == 0) || cnt[i - 1] < cnt[i];
        const bool right_ok = (j == k) || cnt[j] < cnt[i];
        if (left_ok && right_ok) for (int t = i; t < j; ++t) modes.push_back(t);
        i = j;
      }
      if (multi) {
        std::vector<double> out;
        out.reserve(modes.size());
        for (int t : modes) out.push_back(vals[t]);
        return out;
      }
      // Equally frequent modes: nearest the gravity, then the lower value.
      int top = 0;
      for (int t : modes) top = std::max(top, cnt[t]);
      int n_top = 0, only = -1;
      for (int t : modes) {
        if (cnt[t] == top) { ++n_top; only = t; }
      }
      if (n_top == 1) return {vals[only]};
      const double g = gravity(xnum.data(), xnum.size(), false);
      double best_c = kNaN, best_d = std::numeric_limits<double>::infinity();
      for (int t : modes) {
        if (cnt[t] != top) continue;
        const double dist = std::fabs(vals[t] - g);
        if (dist < best_d) { best_d = dist; best_c = vals[t]; }
      }
      return {best_c};
    }
  }

  std::vector<double> z_names;
  std::vector<int> counts;
  simple_bin_counts(xnum, width, xnum.front(), z_names, counts);
  
  const int lz = static_cast<int>(counts.size());
  if (lz == 0) return {kNaN};

  // Centre on the observations in bins lo..hi, a slice of sorted xnum; bin left edges sit below the data.
  std::vector<int> start(lz + 1, 0);
  for (int b = 0; b < lz; ++b) start[b + 1] = start[b] + counts[b];
  auto obs_centre = [&](int lo, int hi) -> double {
    const int a = start[lo], e = start[hi + 1];
    if (e <= a) return z_names[lo];
    if (xnum[a] == xnum[e - 1]) return xnum[a];
    long double s = 0.0L;
    for (int t = a; t < e; ++t) s += static_cast<long double>(xnum[t]);
    return static_cast<double>(s / static_cast<long double>(e - a));
  };
  
  int maxc = 0; 
  for (int c : counts) if (c > maxc) maxc = c;
  
  std::vector<double> cs; 
  smooth_counts_tri7(counts, cs);
  
  const double MARGIN = 0.0;
  std::vector<int> peak_idx; 
  peak_idx.reserve(lz);
  
  // Near an end a peak must be the tallest bin: a density may peak at its minimum, a lone extreme value not.
  const double cs_max = *std::max_element(cs.begin(), cs.end());
  for (int i = 0; i < lz && lz >= 7; ++i) {
    const double ci = cs[i];
    if (ci <= 0.0) continue;
    if ((i < 3 || i > lz - 4) && ci < cs_max) continue;
    
    double Ls = -std::numeric_limits<double>::infinity();
    double Rs = -std::numeric_limits<double>::infinity();
    for (int d = 1; d <= 3; ++d) {
      if (i - d >= 0) Ls = std::max(Ls, cs[i - d]);
      if (i + d < lz) Rs = std::max(Rs, cs[i + d]);
    }
    if (!(ci > Ls + MARGIN && ci > Rs + MARGIN)) continue;
    
    const double left = (i > 0) ? cs[i - 1] : ((i + 1 < lz) ? cs[i + 1] : ci);
    const double right = (i + 1 < lz) ? cs[i + 1] : ((i > 0) ? cs[i - 1] : ci);
    if (!(left - 2.0 * ci + right < 0.0)) continue;
    
    peak_idx.push_back(i);
  }
  
  if (!peak_idx.empty()) {
    // Equal heights keep bin order, so the kept peaks don't depend on the sort implementation.
    std::stable_sort(peak_idx.begin(), peak_idx.end(), [&](int a, int b){ return cs[a] > cs[b]; });
    std::vector<int> kept;
    for (int idx : peak_idx) {
      bool too_close = false;
      for (int jdx : kept) if (std::abs(idx - jdx) <= 3) { too_close = true; break; }
      if (!too_close) kept.push_back(idx);
    }
    
    if (!kept.empty()) {
      std::vector<double> centers(kept.size());
      for (std::size_t t = 0; t < kept.size(); ++t) {
        const int zc = kept[t];
        centers[t] = obs_centre(std::max(0, zc - 3), std::min(lz - 1, zc + 3));
      }
      
      if (multi) {
        std::vector<double> out = centers;
        std::sort(out.begin(), out.end());
        return out;
      } else {
        // Tied peaks: nearest the gravity, then the lower value; an average can fall in a gap.
        double best_h = cs[kept[0]];
        for (std::size_t t = 1; t < kept.size(); ++t) best_h = std::max(best_h, cs[kept[t]]);
        const double g = gravity(xnum.data(), xnum.size(), false);
        double best_c = kNaN, best_d = std::numeric_limits<double>::infinity();
        for (std::size_t t = 0; t < kept.size(); ++t) {
          if (cs[kept[t]] != best_h) continue;
          const double dist = std::fabs(centers[t] - g);
          if (dist < best_d || (dist == best_d && centers[t] < best_c)) { best_d = dist; best_c = centers[t]; }
        }
        return {best_c};
      }
    }
  }
  
  int ties = 0; 
  for (int c : counts) if (c == maxc) ++ties;
  
  if (ties > 1) {
    if (multi) {
      std::vector<double> out;
      out.reserve(ties);
      for (int i = 0; i < lz; ++i) if (counts[i] == maxc) out.push_back(obs_centre(i, i));
      std::sort(out.begin(), out.end());
      return out;
    } else {
      // Tied highest bins: nearest the gravity, then the lower value.
      const double g = gravity(xnum.data(), xnum.size(), false);
      double best_c = kNaN, best_d = std::numeric_limits<double>::infinity();
      for (int i = 0; i < lz; ++i) {
        if (counts[i] != maxc) continue;
        const double c = obs_centre(i, i);
        const double dist = std::fabs(c - g);
        if (dist < best_d || (dist == best_d && c < best_c)) { best_d = dist; best_c = c; }
      }
      return {best_c};
    }
  }
  
  int zc = 0; 
  for (int i = 0; i < lz; ++i) if (counts[i] == maxc) { zc = i; break; }
  return {obs_centre(std::max(0, zc - 1), std::min(lz - 1, zc + 1))};
}

} // namespace nns