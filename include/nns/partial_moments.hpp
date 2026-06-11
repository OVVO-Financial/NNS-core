// include/nns/partial_moments.hpp
//
// SPDX-License-Identifier: GPL-3.0-only
#ifndef NNS_PARTIAL_MOMENTS_HPP
#define NNS_PARTIAL_MOMENTS_HPP

#include <cstddef>
#include <vector>

namespace nns {

// --- Output Data Structures ---

struct PMMatrixResult {
  std::vector<double> cupm; // Co-Upper Partial Moment matrix (column-major)
  std::vector<double> dupm; // Divergent Upper Partial Moment matrix (column-major)
  std::vector<double> dlpm; // Divergent Lower Partial Moment matrix (column-major)
  std::vector<double> clpm; // Co-Lower Partial Moment matrix (column-major)
  std::vector<double> cov;  // Covariance matrix (column-major)
  std::size_t p;            // Number of columns/variables
};

// --- Univariate Partial Moments ---

double lpm(double degree, double target, const double* x, std::size_t n);
double upm(double degree, double target, const double* x, std::size_t n);

// --- Bivariate Co-Partial Moments ---

double co_lpm(double degree, double target_x, double target_y, 
              const double* x, const double* y, std::size_t n);

double co_upm(double degree, double target_x, double target_y, 
              const double* x, const double* y, std::size_t n);

double d_lpm(double degree_lpm, double degree_upm, double target_x, double target_y, 
             const double* x, const double* y, std::size_t n);

double d_upm(double degree_lpm, double degree_upm, double target_x, double target_y, 
             const double* x, const double* y, std::size_t n);

// --- N-Dimensional Partial Moments ---

double clpm_nD(const double* X, std::size_t n, std::size_t p, const double* target, double degree, bool norm);
double cupm_nD(const double* X, std::size_t n, std::size_t p, const double* target, double degree, bool norm);
double dlpm_nD(const double* X, std::size_t n, std::size_t p, const double* target, double degree, bool norm);
double dupm_nD(const double* X, std::size_t n, std::size_t p, const double* target, double degree, bool norm);

// --- Full Matrix Computation ---

/// Compute the full suite of Partial Moment matrices.
///
/// @param degree_lpm Degree for the lower partial moments.
/// @param degree_upm Degree for the upper partial moments.
/// @param target Pointer to the target array (length p).
/// @param X Pointer to the column-major data matrix.
/// @param n Number of rows in X.
/// @param p Number of columns in X.
/// @param pop_adj Apply population adjustment.
/// @param norm Normalize the resulting matrices.
/// @param nthreads Number of threads to use (-1 for hardware max).
/// @return PMMatrixResult containing the 5 resulting flat column-major matrices.
PMMatrixResult pm_matrix(double degree_lpm, double degree_upm, const double* target, 
                         const double* X, std::size_t n, std::size_t p, 
                         bool pop_adj, bool norm, int nthreads = -1);

} // namespace nns

#endif // NNS_PARTIAL_MOMENTS_HPP
