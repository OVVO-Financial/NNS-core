// include/nns/stochastic_dominance.hpp
//
// SPDX-License-Identifier: GPL-3.0-only
#ifndef NNS_STOCHASTIC_DOMINANCE_HPP
#define NNS_STOCHASTIC_DOMINANCE_HPP

#include <cstddef>
#include <vector>

namespace nns {

// --- Univariate Dominance Tests ---
// Return 1 if x dominates y, otherwise 0.

/// First-degree Stochastic Dominance (Univariate)
/// @param x Pointer to the first array.
/// @param y Pointer to the second array.
/// @param n Length of the arrays.
/// @param discrete Treat data as discrete (true) or continuous (false).
int fsd_uni(const double* x, const double* y, std::size_t n, bool discrete);

/// Second-degree Stochastic Dominance (Univariate)
/// @param x Pointer to the first array.
/// @param y Pointer to the second array.
/// @param n Length of the arrays.
int ssd_uni(const double* x, const double* y, std::size_t n);

/// Third-degree Stochastic Dominance (Univariate)
/// @param x Pointer to the first array.
/// @param y Pointer to the second array.
/// @param n Length of the arrays.
int tsd_uni(const double* x, const double* y, std::size_t n);

// --- Multivariate Dominance Filters ---
// Returns a vector of 0-based column indices representing variables 
// that are NOT dominated by any other variable in the matrix.

/// First-degree Stochastic Dominance (Multivariate)
/// @param X Pointer to the column-major matrix data.
/// @param n Number of rows in X.
/// @param p Number of columns in X.
/// @param discrete Treat data as discrete (true) or continuous (false).
/// @param nthreads Number of parallel threads to use (-1 for hardware max).
std::vector<int> fsd(const double* X, std::size_t n, std::size_t p, bool discrete, int nthreads = -1);

/// Second-degree Stochastic Dominance (Multivariate)
/// @param X Pointer to the column-major matrix data.
/// @param n Number of rows in X.
/// @param p Number of columns in X.
/// @param nthreads Number of parallel threads to use (-1 for hardware max).
std::vector<int> ssd(const double* X, std::size_t n, std::size_t p, int nthreads = -1);

/// Third-degree Stochastic Dominance (Multivariate)
/// @param X Pointer to the column-major matrix data.
/// @param n Number of rows in X.
/// @param p Number of columns in X.
/// @param nthreads Number of parallel threads to use (-1 for hardware max).
std::vector<int> tsd(const double* X, std::size_t n, std::size_t p, int nthreads = -1);

// --- Stochastic Superiority ---

struct StochSupResult {
  double p_x_gt_y; // Probability that X > Y
  double p_x_eq_y; // Probability that X == Y
  double p_x_lt_y; // Probability that X < Y
};

/// Compute the stochastic superiority of array X over array Y.
///
/// @param x Pointer to the first numeric array (X).
/// @param n_x Length of array X.
/// @param y Pointer to the second numeric array (Y).
/// @param n_y Length of array Y.
/// @return StochSupResult containing the exact probabilities.
StochSupResult stochastic_superiority(const double* x, std::size_t n_x, 
                                      const double* y, std::size_t n_y);

} // namespace nns

#endif // NNS_STOCHASTIC_DOMINANCE_HPP