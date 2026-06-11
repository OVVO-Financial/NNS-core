// include/nns/fast_lm.hpp
//
// SPDX-License-Identifier: GPL-3.0-only
#ifndef NNS_FAST_LM_HPP
#define NNS_FAST_LM_HPP

#include <cstddef>
#include <vector>

namespace nns {

struct FastLmResult {
  std::vector<double> coef;          // length 2: [intercept, slope]
  std::vector<double> fitted_values; // length n
  std::vector<double> residuals;     // length n
  int df_residual;
};

struct FastLmMultResult {
  std::vector<double> coefficients;  // length p + 1 (intercept + slopes)
  std::vector<double> fitted_values; // length n
  std::vector<double> residuals;     // length n
  double r_squared;
};

/// Fast univariate ordinary least squares
///
/// @param x Pointer to the predictor array.
/// @param y Pointer to the response array.
/// @param n Length of the arrays.
/// @return FastLmResult struct containing coefficients, fitted values, and residuals.
FastLmResult fast_lm(const double* x, const double* y, std::size_t n);

/// Fast multivariate ordinary least squares
///
/// @param X Pointer to the column-major predictor matrix (n x p).
/// @param y Pointer to the response array.
/// @param n Number of rows in X.
/// @param p Number of columns in X.
/// @return FastLmMultResult struct containing coefficients, fitted values, residuals, and R^2.
FastLmMultResult fast_lm_mult(const double* X, const double* y, std::size_t n, std::size_t p);

} // namespace nns

#endif // NNS_FAST_LM_HPP