// src/seasonality.cpp
//
// Implementation extracted from NNS NNS_seas.cpp (synced to R 13.3, 79c60272).
// Decoupled from Rcpp.
//
// SPDX-License-Identifier: GPL-3.0-only
#include "nns/seasonality.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>

namespace nns {

namespace {

inline bool any_na_or_inf(const double* x, std::size_t n) {
  for (std::size_t i = 0; i < n; ++i) {
    if (!std::isfinite(x[i])) return true;
  }
  return false;
}

} // namespace

// ---------- Core Seasonality API ----------

SeasonalityResult seasonality(const double* x, std::size_t n,
                              const int* modulo, std::size_t mod_len,
                              bool mod_only, std::optional<int> max_period) {
  if (n == 0) throw std::invalid_argument("Variable must be numeric and non-empty");
  if (any_na_or_inf(x, n)) throw std::invalid_argument("You have some missing or infinite values, please address.");

  if (n < 5) {
    return {
      {0}, {0.0}, {0.0}, // all.periods (DataFrame cols)
      0,                 // best.period
      {0}                // periods (upstream returns c(0), not empty)
    };
  }

  const int ni = static_cast<int>(n);
  const int half_n = ni / 2;

  // Gini mean difference: unbiased at any component length, and one sort serves every lag.
  std::vector<int> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  std::stable_sort(ord.begin(), ord.end(), [&](int a, int b) { return x[a] < x[b]; });

  // CV needs a ratio scale: used only when the series stays on one side of 0.
  const bool use_cv = (x[ord[0]] > 0.0) || (x[ord[n - 1]] < 0.0);

  // Every component series of lag m keeps at least kMinComponentLength points.
  constexpr int kMinComponentLength = 5;
  int max_m = std::min(half_n, ni / kMinComponentLength);
  if (max_period) max_m = std::max(1, std::min(max_m, *max_period));

  const std::size_t buf = static_cast<std::size_t>(half_n) + 1;
  std::vector<int> rank(buf), len(buf);
  std::vector<double> sum(buf), acc(buf), first(buf), last(buf);
  // Length-weighted mean CV (or Gini) of the m component series x[j], x[j + m], ...
  auto score = [&](int m) {
    const double inv_m = 1.0 / static_cast<double>(m);
    for (int j = 0; j < m; ++j) { rank[j] = 0; len[j] = (ni - 1 - j) / m + 1; sum[j] = 0.0; acc[j] = 0.0; }
    for (int t = 0; t < ni; ++t) {
      const int i = ord[t];
      int j = i - m * static_cast<int>(static_cast<double>(i) * inv_m);
      if (j >= m) j -= m; else if (j < 0) j += m;
      const double v = x[i];
      if (rank[j] == 0) first[j] = v;
      last[j] = v;
      // The weights sum to 0, so measuring from the component minimum avoids cancelling a large level.
      acc[j] += (2.0 * rank[j] + 1.0 - len[j]) * (v - first[j]);
      sum[j] += v;
      ++rank[j];
    }
    double wsum = 0.0;
    for (int j = 0; j < m; ++j) {
      const double k = static_cast<double>(len[j]);
      const double g = (len[j] < 2 || first[j] == last[j]) ? 0.0 : std::max(0.0, 2.0 * acc[j] / (k * (k - 1.0)));
      wsum += k * (use_cv ? (g == 0.0 ? 0.0 : g / std::fabs(sum[j] / k)) : g);
    }
    return wsum / static_cast<double>(ni);
  };
  // The series itself is the single component of lag 1.
  const double var_cov = score(1);

  // Scored past the scan limit too, so ARMA.seas.weighting doesn't give it zero weight.
  auto added_cv = [&](int m) {
    if (m > max_m && m <= half_n) {
      const double c = score(m);
      if (c < var_cov) return c;
    }
    return var_cov;
  };

  std::vector<int> periods_vec;
  std::vector<double> cvmean_vec;
  for (int m = 1; m <= max_m; ++m) {
    const double cv_m = score(m);
    if (cv_m < var_cov) {
      periods_vec.push_back(m);
      cvmean_vec.push_back(cv_m);
    }
  }

  std::vector<int> Period;
  std::vector<double> CoefVar;
  std::vector<double> VarCoefVar;

  if (!periods_vec.empty()) {
    int m = static_cast<int>(periods_vec.size());
    Period = periods_vec;
    CoefVar = cvmean_vec;
    VarCoefVar.assign(m, var_cov);

    std::vector<int> ord(m);
    std::iota(ord.begin(), ord.end(), 0);
    // Equal CVs keep the shorter period first.
    std::stable_sort(ord.begin(), ord.end(), [&](int a, int b) { return CoefVar[a] < CoefVar[b]; });

    std::vector<int> sortedP(m);
    std::vector<double> sortedCV(m);
    for (int k = 0; k < m; ++k) {
      sortedP[k] = Period[ord[k]];
      sortedCV[k] = CoefVar[ord[k]];
    }
    Period = std::move(sortedP);
    CoefVar = std::move(sortedCV);
  } else {
    Period = {1};
    CoefVar = {var_cov};
    VarCoefVar = {var_cov};
  }

  // Modulo Handling
  if (modulo != nullptr && mod_len > 0) {
    std::set<int> per_set;
    for (std::size_t i = 0; i < Period.size(); ++i) {
      for (std::size_t j = 0; j < mod_len; ++j) {
        int m_val = modulo[j];
        if (m_val <= 0) continue;
        int minus = Period[i] - (Period[i] % m_val);
        int plus  = Period[i] + (m_val - (Period[i] % m_val));
        if (minus > 0) per_set.insert(minus);
        if (plus > 0) per_set.insert(plus);
      }
    }

    if (mod_only) {
      std::set<int> curr(Period.begin(), Period.end());
      std::vector<int> keptP;
      std::vector<double> keptCV;

      for (std::size_t i = 0; i < Period.size(); ++i) {
        if (per_set.count(Period[i])) {
          keptP.push_back(Period[i]);
          keptCV.push_back(CoefVar[i]);
        }
      }
      for (int s : per_set) {
        if (!curr.count(s)) {
          keptP.push_back(s);
          keptCV.push_back(added_cv(s));
        }
      }

      if (keptP.empty()) {
        Period = {1};
        CoefVar = {var_cov};
        VarCoefVar = {var_cov};
      } else {
        int m = static_cast<int>(keptP.size());
        Period = keptP;
        CoefVar = keptCV;
        VarCoefVar.assign(m, var_cov);

        std::vector<int> ord(m);
        std::iota(ord.begin(), ord.end(), 0);
        std::stable_sort(ord.begin(), ord.end(), [&](int a, int b) { return CoefVar[a] < CoefVar[b]; });

        std::vector<int> sortedP(m);
        std::vector<double> sortedCV(m);
        for (int k = 0; k < m; ++k) {
          sortedP[k] = Period[ord[k]];
          sortedCV[k] = CoefVar[ord[k]];
        }
        Period = std::move(sortedP);
        CoefVar = std::move(sortedCV);
      }
    } else {
      per_set.insert(1);
      std::set<int> curr(Period.begin(), Period.end());
      std::vector<int> add;
      for (int s : per_set) {
        if (!curr.count(s)) add.push_back(s);
      }

      if (!add.empty()) {
        for (int a : add) {
          Period.push_back(a);
          CoefVar.push_back(added_cv(a));
          VarCoefVar.push_back(var_cov);
        }

        int m = static_cast<int>(Period.size());
        std::vector<int> ord(m);
        std::iota(ord.begin(), ord.end(), 0);
        std::stable_sort(ord.begin(), ord.end(), [&](int a, int b) { return CoefVar[a] < CoefVar[b]; });

        std::vector<int> sortedP(m);
        std::vector<double> sortedCV(m);
        for (int k = 0; k < m; ++k) {
          sortedP[k] = Period[ord[k]];
          sortedCV[k] = CoefVar[ord[k]];
        }
        Period = std::move(sortedP);
        CoefVar = std::move(sortedCV);
      }
    }
  }

  // Strict cap: Period < n/2
  {
    std::vector<int> P;
    std::vector<double> CV;
    std::vector<double> VCV;
    for (std::size_t i = 0; i < Period.size(); ++i) {
      if (Period[i] < static_cast<int>(n) / 2) {
        P.push_back(Period[i]);
        CV.push_back(CoefVar[i]);
        VCV.push_back(VarCoefVar[i]);
      }
    }

    if (!P.empty()) {
      int m = static_cast<int>(P.size());
      Period = P;
      CoefVar = CV;
      VarCoefVar = VCV;

      std::vector<int> ord(m);
      std::iota(ord.begin(), ord.end(), 0);
      std::stable_sort(ord.begin(), ord.end(), [&](int a, int b) { return CoefVar[a] < CoefVar[b]; });

      std::vector<int> sortedP(m);
      std::vector<double> sortedCV(m);
      for (int k = 0; k < m; ++k) {
        sortedP[k] = Period[ord[k]];
        sortedCV[k] = CoefVar[ord[k]];
      }
      Period = std::move(sortedP);
      CoefVar = std::move(sortedCV);
    } else {
      Period = {1};
      CoefVar = {var_cov};
      VarCoefVar = {var_cov};
    }
  }

  SeasonalityResult res;
  res.all_periods = Period;
  res.all_coef_var = CoefVar;
  res.all_var_coef_var = VarCoefVar;
  res.best_period = Period.empty() ? 0 : Period[0];
  res.periods = Period;

  return res;
}

} // namespace nns