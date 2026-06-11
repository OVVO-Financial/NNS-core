// src/numerics.cpp
//
// SPDX-License-Identifier: GPL-3.0-only
#include "nns/numerics.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

namespace nns {

namespace {
  inline double at(const double* M, std::size_t rows, std::size_t r, std::size_t c) {
    return M[c * rows + r];
  }

  // Gaussian CDF inverse approximation for CLT
  double qnorm_approx(double p) {
    // A standard rational approximation for the inverse normal CDF
    if (p <= 0.0) return -8.0;
    if (p >= 1.0) return 8.0;
    double c0 = 2.515517, c1 = 0.802853, c2 = 0.010328;
    double d1 = 1.432788, d2 = 0.189269, d3 = 0.001308;
    double t = std::sqrt(-2.0 * std::log(std::min(p, 1.0 - p)));
    double val = t - ((c2 * t + c1) * t + c0) / (((d3 * t + d2) * t + d1) * t + 1.0);
    return (p < 0.5) ? -val : val;
  }
}

// ---------- Basic Utilities ----------

double vec_sd(const double* x, std::size_t n) {
  if (n <= 1) return std::numeric_limits<double>::quiet_NaN();
  double mu = 0.0;
  for (std::size_t i = 0; i < n; ++i) mu += x[i];
  mu /= static_cast<double>(n);
  double ss = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    double d = x[i] - mu;
    ss += d * d;
  }
  return std::sqrt(ss / static_cast<double>(n - 1));
}

std::vector<double> col_sd(const double* X, std::size_t n, std::size_t p) {
  std::vector<double> sds(p, std::numeric_limits<double>::quiet_NaN());
  if (n <= 1) return sds;
  
  for (std::size_t j = 0; j < p; ++j) {
    double mu = 0.0;
    for (std::size_t i = 0; i < n; ++i) mu += at(X, n, i, j);
    mu /= static_cast<double>(n);
    double ss = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      double d = at(X, n, i, j) - mu;
      ss += d * d;
    }
    sds[j] = std::sqrt(ss / static_cast<double>(n - 1));
  }
  return sds;
}

bool is_discrete(const double* x, std::size_t n) {
  for (std::size_t i = 0; i < n; ++i) {
    if (std::isfinite(x[i]) && x[i] != std::trunc(x[i])) return false;
  }
  return true;
}

// ---------- Vector Generation ----------

TimeSeriesVectors generate_vectors(const double* x, std::size_t n, const int* lags, std::size_t num_lags) {
  TimeSeriesVectors res;
  res.series.resize(num_lags);
  res.index.resize(num_lags);
  
  for (std::size_t t = 0; t < num_lags; ++t) {
    int lag = lags[t];
    if (lag <= 0) continue;
    
    int start = (n % lag) + 1;
    int m = ((n - start) / lag) + 1;
    
    std::vector<double> s(m);
    std::vector<int> idx(m);
    
    int pos = start;
    for (int i = 0; i < m; ++i, pos += lag) {
      s[i] = x[pos - 1];
      idx[i] = i + 1;
    }
    res.series[t] = std::move(s);
    res.index[t] = std::move(idx);
  }
  return res;
}

ForecastVectors generate_lin_vectors(const double* x, std::size_t n, int l, int h) {
  int max_fcast = std::min(h, l);
  ForecastVectors res;
  res.series.resize(max_fcast);
  res.index.resize(max_fcast);
  
  for (int i = 1; i <= max_fcast; ++i) {
    int start = ((n + i - 1) % l) + 1;
    int m = ((n - start) / l) + 1;
    std::vector<double> s(m);
    std::vector<int> idx(m);
    int pos = start;
    for (int k = 0; k < m; ++k, pos += l) {
      s[k] = x[pos - 1];
      idx[k] = k + 1;
    }
    res.series[i - 1] = std::move(s);
    res.index[i - 1] = std::move(idx);
  }
  
  res.forecast_index.resize(max_fcast);
  for (int i = 0; i < h; ++i) {
    res.forecast_index[i % max_fcast].push_back(i + 1);
  }
  
  res.forecast_values.resize(l);
  for (int i = 1; i <= h; ++i) {
    int ci = ((((i - 1) % l)) % std::max(1, max_fcast)) + 1;
    int last_val = res.index[ci - 1].size();
    double fval = static_cast<double>(last_val) + std::ceil(static_cast<double>(i) / static_cast<double>(l));
    res.forecast_values[(i - 1) % l].push_back(fval);
  }
  return res;
}

// ---------- ARMA Weighting ----------

ARMAWeights arma_seas_weighting(const double* periods, const double* covar, const double* varcovar, std::size_t m) {
  if (m == 0) return {{1.0}, {1.0}};
  
  std::vector<double> obs_weight(m);
  for (std::size_t i = 0; i < m; ++i) obs_weight[i] = 1.0 / std::sqrt(periods[i]);
  
  std::vector<double> lag_weight(m, 1.0);
  if (covar != nullptr && varcovar != nullptr) {
    for (std::size_t i = 0; i < m; ++i) lag_weight[i] = varcovar[i] - covar[i];
  }
  
  std::vector<double> wprod(m);
  double denom = 0.0;
  for (std::size_t i = 0; i < m; ++i) {
    wprod[i] = lag_weight[i] * obs_weight[i];
    denom += wprod[i];
  }
  
  ARMAWeights res;
  res.lags.assign(periods, periods + m);
  if (denom == 0.0) {
    res.weights.assign(m, 0.0);
  } else {
    res.weights.resize(m);
    for (std::size_t i = 0; i < m; ++i) res.weights[i] = wprod[i] / denom;
  }
  return res;
}

// ---------- MEBOOT Core ----------

std::vector<double> meboot_part(const double* xx, std::size_t m, std::size_t n, 
                                const double* z, std::size_t z_len,
                                double xmin, double xmax, 
                                const double* desintxb, bool reachbnd, int seed) {
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  
  std::vector<double> p(n);
  for (std::size_t i = 0; i < n; ++i) p[i] = dist(gen);
  
  std::vector<double> q(n, std::numeric_limits<double>::quiet_NaN());
  if (m == 1) {
    std::fill(q.begin(), q.end(), xx[0]);
  } else if (m > 1) {
    for (std::size_t i = 0; i < n; ++i) {
      double pi = p[i];
      if (pi <= 0.0) { q[i] = xx[0]; continue; }
      if (pi >= 1.0) { q[i] = xx[m - 1]; continue; }
      double h = 1.0 + (m - 1.0) * pi;
      int j = static_cast<int>(std::floor(h));
      if (j < 1) j = 1; else if (j > static_cast<int>(m) - 1) j = m - 1;
      q[i] = (1.0 - (h - j)) * xx[j - 1] + (h - j) * xx[j];
    }
  }
  
  double invn = 1.0 / static_cast<double>(n);
  double edge = static_cast<double>(n - 1) / static_cast<double>(n);
  
  for (std::size_t i = 0; i < n; ++i) {
    if (p[i] <= invn) {
      double val = xmin + (p[i] - 0.0) * (z[0] - xmin) / (invn - 0.0);
      if (!reachbnd) val = val + desintxb[0] - 0.5 * (z[0] + xmin);
      q[i] = val;
    } else if (p[i] >= edge) {
      double val = z[z_len - 2] + (p[i] - edge) * (xmax - z[z_len - 2]) / (1.0 - edge);
      if (!reachbnd) val = val + desintxb[m - 1] - 0.5 * (z[z_len - 2] + xmax);
      q[i] = val;
    }
  }
  return q;
}

void meboot_expand_sd(double* ensemble, std::size_t n, std::size_t J, 
                      const double* orig_sd, std::size_t orig_p, double fiv, int seed) {
  std::vector<double> ens_sd = col_sd(ensemble, n, J);
  std::vector<double> sdf;
  sdf.reserve(orig_p + J);
  for (std::size_t i = 0; i < orig_p; ++i) sdf.push_back(orig_sd[i]);
  for (std::size_t j = 0; j < J; ++j) sdf.push_back(ens_sd[j]);
  
  std::vector<double> sdfa(sdf.size()), sdfd(sdf.size());
  for (std::size_t i = 0; i < sdf.size(); ++i) {
    sdfa[i] = sdf[i] / sdf[0];
    sdfd[i] = sdf[0] / sdf[i];
  }
  
  std::mt19937 gen(seed);
  double mx = 1.0 + (fiv / 100.0);
  std::uniform_real_distribution<double> dist(1.0, mx);
  
  for (std::size_t i = 0; i < sdfa.size(); ++i) {
    if (sdfa[i] < 1.0) sdfa[i] = dist(gen);
  }
  
  for (std::size_t j = 0; j < J; ++j) {
    double a = sdfd[j + 1] * sdfa[j + 1];
    if (std::floor(a) > 0.0) {
      for (std::size_t i = 0; i < n; ++i) {
        ensemble[j * n + i] *= a;
      }
    }
  }
}

void force_clt(double* ensemble, std::size_t n, std::size_t J, 
               double orig_gm, const double* orig_sd, std::size_t orig_p) {
  std::vector<double> xbar(J);
  for (std::size_t j = 0; j < J; ++j) {
    double mu = 0.0;
    for (std::size_t i = 0; i < n; ++i) mu += at(ensemble, n, i, j);
    xbar[j] = mu / static_cast<double>(n);
  }
  
  std::vector<int> oo(J);
  std::iota(oo.begin(), oo.end(), 0);
  std::sort(oo.begin(), oo.end(), [&](int a, int b){ return xbar[a] < xbar[b]; });
  
  std::vector<double> sortxbar = xbar;
  std::sort(sortxbar.begin(), sortxbar.end());
  
  std::vector<double> smean(orig_p);
  for (std::size_t i = 0; i < orig_p; ++i) smean[i] = orig_sd[i] / std::sqrt(static_cast<double>(J));
  double smean_scalar = smean.empty() ? 0.0 : smean[0];
  
  std::vector<double> newbar(J);
  for (std::size_t j = 0; j < J; ++j) {
    double sm = (orig_p == 1) ? smean_scalar : smean[j % orig_p];
    newbar[j] = orig_gm + qnorm_approx(static_cast<double>(j + 1) / static_cast<double>(J + 1)) * sm;
  }
  
  double mu_nb = 0.0, ss_nb = 0.0;
  for (std::size_t j = 0; j < J; ++j) mu_nb += newbar[j];
  mu_nb /= static_cast<double>(J);
  for (std::size_t j = 0; j < J; ++j) ss_nb += (newbar[j] - mu_nb) * (newbar[j] - mu_nb);
  double sd_nb = std::sqrt(ss_nb / static_cast<double>(J - 1));
  
  std::vector<double> out_ensemble(n * J);
  for (std::size_t i = 0; i < J; ++i) {
    int col = oo[i];
    double sm = (orig_p == 1) ? smean_scalar : smean[i % orig_p];
    double add = (((newbar[i] - mu_nb) / sd_nb) * sm + orig_gm) - sortxbar[i];
    for (std::size_t r = 0; r < n; ++r) {
      out_ensemble[col * n + r] = at(ensemble, n, r, col) + add;
    }
  }
  std::copy(out_ensemble.begin(), out_ensemble.end(), ensemble);
}

// ---------- Class Resampling ----------

SampleResult down_sample(const double* X, const int* y, std::size_t n, std::size_t p, int seed) {
  if (n == 0) return {std::vector<double>(), std::vector<int>(), 0, p};
  
  std::map<int, std::vector<std::size_t>> per_class;
  for (std::size_t i = 0; i < n; ++i) per_class[y[i]].push_back(i);
  
  std::size_t min_class = n;
  for (const auto& kv : per_class) {
    if (kv.second.size() < min_class && !kv.second.empty()) min_class = kv.second.size();
  }
  
  if (min_class == n || min_class == 0) throw std::invalid_argument("down_sample: no valid class distribution.");
  
  std::vector<std::size_t> rows_out;
  rows_out.reserve(per_class.size() * min_class);
  
  std::mt19937 gen(static_cast<unsigned int>(seed));
  for (const auto& kv : per_class) {
    std::vector<std::size_t> indices = kv.second;
    std::shuffle(indices.begin(), indices.end(), gen);
    for (std::size_t i = 0; i < min_class; ++i) rows_out.push_back(indices[i]);
  }
  
  std::size_t new_n = rows_out.size();
  SampleResult res{std::vector<double>(new_n * p), std::vector<int>(new_n), new_n, p};
  
  for (std::size_t i = 0; i < new_n; ++i) {
    std::size_t orig_row = rows_out[i];
    res.y[i] = y[orig_row];
    for (std::size_t j = 0; j < p; ++j) res.x[j * new_n + i] = at(X, n, orig_row, j);
  }
  return res;
}

SampleResult up_sample(const double* X, const int* y, std::size_t n, std::size_t p, int seed) {
  if (n == 0) return {std::vector<double>(), std::vector<int>(), 0, p};
  
  std::map<int, std::vector<std::size_t>> per_class;
  for (std::size_t i = 0; i < n; ++i) per_class[y[i]].push_back(i);
  
  std::size_t max_class = 0;
  for (const auto& kv : per_class) {
    if (kv.second.size() > max_class) max_class = kv.second.size();
  }
  
  if (max_class == 0) throw std::invalid_argument("up_sample: no valid class distribution.");
  
  std::vector<std::size_t> rows_out;
  rows_out.reserve(per_class.size() * max_class);
  std::mt19937 gen(static_cast<unsigned int>(seed));
  
  for (const auto& kv : per_class) {
    const auto& indices = kv.second;
    std::size_t sz = indices.size();
    for (std::size_t i = 0; i < sz; ++i) rows_out.push_back(indices[i]);
    
    std::size_t needed = max_class - sz;
    if (needed > 0) {
      std::uniform_int_distribution<std::size_t> dist(0, sz - 1);
      for (std::size_t i = 0; i < needed; ++i) rows_out.push_back(indices[dist(gen)]);
    }
  }
  
  std::size_t new_n = rows_out.size();
  SampleResult res{std::vector<double>(new_n * p), std::vector<int>(new_n), new_n, p};
  
  for (std::size_t i = 0; i < new_n; ++i) {
    std::size_t orig_row = rows_out[i];
    res.y[i] = y[orig_row];
    for (std::size_t j = 0; j < p; ++j) res.x[j * new_n + i] = at(X, n, orig_row, j);
  }
  return res;
}

} // namespace nns