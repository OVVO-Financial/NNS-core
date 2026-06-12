// src/stochastic_dominance.cpp
//
// Implementation reconstructed from original_src/SD.cpp and original_src/stoch_sup.cpp; covers FSD/SSD/TSD and stochastic superiority. Decoupled from Rcpp.
//
// SPDX-License-Identifier: GPL-3.0-only
#include "nns/stochastic_dominance.hpp"
#include "nns/parallel.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace nns {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

inline double at(const double* M, std::size_t rows, std::size_t r, std::size_t c) {
  return M[c * rows + r];
}

// =====================================================================
// Per-column precompute: sorted values, prefix sums, basic stats
// =====================================================================
struct ColPre {
  std::vector<double> vals;      // sorted ascending, length m
  std::vector<double> P1;        // prefix sum of vals; length m+1, P1[0]=0
  std::vector<double> P2;        // prefix sum of vals^2; length m+1
  double S1{0.0}, S2{0.0};
  double mn{std::numeric_limits<double>::infinity()}, mean{kNaN};
  std::size_t m{0};
};

ColPre precompute_col(const double* X, std::size_t n, std::size_t p, std::size_t j) {
  ColPre c; 
  c.m = n;
  c.vals.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    double v = at(X, n, i, j);
    if (!std::isfinite(v)) throw std::invalid_argument("NAs or non-finite values not allowed in SD");
    c.vals[i] = v;
  }
  std::sort(c.vals.begin(), c.vals.end());
  
  c.P1.assign(n + 1, 0.0);
  c.P2.assign(n + 1, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    double v = c.vals[i];
    c.P1[i + 1] = c.P1[i] + v;
    c.P2[i + 1] = c.P2[i] + v * v;
  }
  c.S1 = c.P1[n];
  c.S2 = c.P2[n];
  if (n > 0) {
    c.mn = c.vals[0];
    c.mean = c.S1 / static_cast<double>(n);
  }
  return c;
}

ColPre precompute_vec(const double* x, std::size_t n) {
  ColPre c; 
  c.m = n;
  c.vals.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    if (!std::isfinite(x[i])) throw std::invalid_argument("NAs or non-finite values not allowed in SD");
    c.vals[i] = x[i];
  }
  std::sort(c.vals.begin(), c.vals.end());
  
  c.P1.assign(n + 1, 0.0);
  c.P2.assign(n + 1, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    double v = c.vals[i];
    c.P1[i + 1] = c.P1[i] + v;
    c.P2[i + 1] = c.P2[i] + v * v;
  }
  c.S1 = c.P1[n];
  c.S2 = c.P2[n];
  if (n > 0) {
    c.mn = c.vals[0];
    c.mean = c.S1 / static_cast<double>(n);
  }
  return c;
}

// =====================================================================
// Core Pairwise Check
// degree = 1 (FSD), 2 (SSD), 3 (TSD)
// Returns 1 if X dominates Y, else 0
// =====================================================================
int sd_dom_pair(const ColPre& X, const ColPre& Y, int degree, bool discrete = false) {
  if (X.m == 0 || Y.m == 0) return 0;
  
  if (degree == 1) {
    if (X.mn < Y.mn) return 0;
    std::size_t nX = X.m, nY = Y.m;
    for (std::size_t i = 0; i < nY; ++i) {
      double t = Y.vals[i];
      auto itX = std::upper_bound(X.vals.begin(), X.vals.end(), t);
      auto itY = std::upper_bound(Y.vals.begin(), Y.vals.end(), t);
      int countX = static_cast<int>(std::distance(X.vals.begin(), itX));
      int countY = static_cast<int>(std::distance(Y.vals.begin(), itY));
      double cdfX = static_cast<double>(countX) / nX;
      double cdfY = static_cast<double>(countY) / nY;
      
      if (discrete) {
        if (cdfX > cdfY) return 0;
      } else {
        if (cdfX > cdfY && std::abs(cdfX - cdfY) > 1e-7) return 0;
      }
    }
    return 1;
  }
  
  if (degree == 2) {
    if (X.mn < Y.mn) return 0;
    if (X.mean < Y.mean) return 0;
    
    std::size_t nX = X.m, nY = Y.m;
    for (std::size_t i = 0; i < nY; ++i) {
      double t = Y.vals[i];
      auto itX = std::upper_bound(X.vals.begin(), X.vals.end(), t);
      int countX = static_cast<int>(std::distance(X.vals.begin(), itX));
      
      double areaX = (static_cast<double>(countX) * t - X.P1[countX]) / nX;
      double areaY = (static_cast<double>(i + 1) * t - Y.P1[i + 1]) / nY;
      
      if (areaX > areaY && std::abs(areaX - areaY) > 1e-7) return 0;
    }
    return 1;
  }
  
  if (degree == 3) {
    if (X.mn < Y.mn) return 0;
    if (X.mean < Y.mean) return 0;
    
    std::size_t nX = X.m, nY = Y.m;
    for (std::size_t i = 0; i < nY; ++i) {
      double t = Y.vals[i];
      auto itX = std::upper_bound(X.vals.begin(), X.vals.end(), t);
      int countX = static_cast<int>(std::distance(X.vals.begin(), itX));
      
      double term1X = t * t * countX;
      double term2X = 2.0 * t * X.P1[countX];
      double term3X = X.P2[countX];
      double volX   = (term1X - term2X + term3X) / (2.0 * nX);
      
      double term1Y = t * t * (i + 1);
      double term2Y = 2.0 * t * Y.P1[i + 1];
      double term3Y = Y.P2[i + 1];
      double volY   = (term1Y - term2Y + term3Y) / (2.0 * nY);
      
      if (volX > volY && std::abs(volX - volY) > 1e-7) return 0;
    }
    return 1;
  }
  
  return 0;
}

} // namespace

// ---------- Univariate Wrappers ----------

int fsd_uni(const double* x, const double* y, std::size_t n, bool discrete) {
  ColPre X = precompute_vec(x, n);
  ColPre Y = precompute_vec(y, n);
  return sd_dom_pair(X, Y, 1, discrete);
}

int ssd_uni(const double* x, const double* y, std::size_t n) {
  ColPre X = precompute_vec(x, n);
  ColPre Y = precompute_vec(y, n);
  return sd_dom_pair(X, Y, 2, false);
}

int tsd_uni(const double* x, const double* y, std::size_t n) {
  ColPre X = precompute_vec(x, n);
  ColPre Y = precompute_vec(y, n);
  return sd_dom_pair(X, Y, 3, false);
}

// ---------- Multivariate Wrappers ----------

std::vector<int> fsd(const double* X, std::size_t n, std::size_t p, bool discrete, int nthreads) {
  std::vector<ColPre> cols(p);
  for (std::size_t j = 0; j < p; ++j) cols[j] = precompute_col(X, n, p, j);
  
  std::vector<int> dom(p * p, 0);
  parallel_for(0, p, [&](std::size_t begin, std::size_t end) {
    for (std::size_t j = begin; j < end; ++j) {
      for (std::size_t k = 0; k < p; ++k) {
        if (j != k) {
          dom[j * p + k] = sd_dom_pair(cols[j], cols[k], 1, discrete);
        }
      }
    }
  }, nthreads);
  
  std::vector<int> kept;
  for (std::size_t j = 0; j < p; ++j) {
    bool dominated = false;
    for (std::size_t k = 0; k < p; ++k) {
      if (j != k && dom[k * p + j] == 1) { dominated = true; break; }
    }
    if (!dominated) kept.push_back(static_cast<int>(j));
  }
  return kept;
}

std::vector<int> ssd(const double* X, std::size_t n, std::size_t p, int nthreads) {
  std::vector<ColPre> cols(p);
  for (std::size_t j = 0; j < p; ++j) cols[j] = precompute_col(X, n, p, j);
  
  std::vector<int> dom(p * p, 0);
  parallel_for(0, p, [&](std::size_t begin, std::size_t end) {
    for (std::size_t j = begin; j < end; ++j) {
      for (std::size_t k = 0; k < p; ++k) {
        if (j != k) {
          dom[j * p + k] = sd_dom_pair(cols[j], cols[k], 2, false);
        }
      }
    }
  }, nthreads);
  
  std::vector<int> kept;
  for (std::size_t j = 0; j < p; ++j) {
    bool dominated = false;
    for (std::size_t k = 0; k < p; ++k) {
      if (j != k && dom[k * p + j] == 1) { dominated = true; break; }
    }
    if (!dominated) kept.push_back(static_cast<int>(j));
  }
  return kept;
}

std::vector<int> tsd(const double* X, std::size_t n, std::size_t p, int nthreads) {
  std::vector<ColPre> cols(p);
  for (std::size_t j = 0; j < p; ++j) cols[j] = precompute_col(X, n, p, j);
  
  std::vector<int> dom(p * p, 0);
  parallel_for(0, p, [&](std::size_t begin, std::size_t end) {
    for (std::size_t j = begin; j < end; ++j) {
      for (std::size_t k = 0; k < p; ++k) {
        if (j != k) {
          dom[j * p + k] = sd_dom_pair(cols[j], cols[k], 3, false);
        }
      }
    }
  }, nthreads);
  
  std::vector<int> kept;
  for (std::size_t j = 0; j < p; ++j) {
    bool dominated = false;
    for (std::size_t k = 0; k < p; ++k) {
      if (j != k && dom[k * p + j] == 1) { dominated = true; break; }
    }
    if (!dominated) kept.push_back(static_cast<int>(j));
  }
  return kept;
}

// ---------- Stochastic Superiority ----------

StochSupResult stochastic_superiority(const double* x, std::size_t n_x, 
                                      const double* y, std::size_t n_y) {
  if (n_x == 0 || n_y == 0) {
    throw std::invalid_argument("stochastic_superiority: x and y must both have positive length.");
  }
  
  // Clone and sort the arrays natively
  std::vector<double> xs(x, x + n_x);
  std::vector<double> ys(y, y + n_y);
  
  std::sort(xs.begin(), xs.end());
  std::sort(ys.begin(), ys.end());
  
  long double less_count = 0.0L;
  long double tie_count  = 0.0L;
  
  std::size_t left  = 0;   // number of elements in y strictly less than x[i]
  std::size_t right = 0;   // number of elements in y less than or equal to x[i]
  
  for (std::size_t i = 0; i < n_x; ++i) {
    const double xi = xs[i];
    
    // Advance 'left' pointer for strictly less
    while (left < n_y && ys[left] < xi) {
      ++left;
    }
    // Advance 'right' pointer for less than or equal
    while (right < n_y && ys[right] <= xi) {
      ++right;
    }
    
    less_count += left;
    tie_count  += (right - left);
  }
  
  const long double denom = static_cast<long double>(n_x) * static_cast<long double>(n_y);
  
  const double p_gt = static_cast<double>(less_count / denom);
  const double p_tie = static_cast<double>(tie_count / denom);
  const double p_star = p_gt + 0.5 * p_tie;
  
  return {p_gt, p_tie, p_star};
}

} // namespace nns