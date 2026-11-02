// include/nns/dependence.hpp
//
// SPDX-License-Identifier: GPL-3.0-only
#ifndef NNS_DEPENDENCE_HPP
#define NNS_DEPENDENCE_HPP

#include <cstddef>
#include <string>
#include <vector>

namespace nns {

// --- Output Data Structures ---

struct DepResult {
  double correlation;
  double dependence;
};

struct DepMatrixResult {
  std::vector<double> correlation; // p x p column-major matrix
  std::vector<double> dependence;  // p x p column-major matrix
  std::size_t p;
};

// --- Core API ---

/// Bivariate dependence from pre-computed partition labels (R NNS_dep_pair_cpp).
///
/// Labels are the prefix-coded quadrants of nns::partition ("q", "q1",
/// "q12", ...); every level of each partition is read from the prefixes.
///
/// @param x Pointer to the first array.
/// @param y Pointer to the second array.
/// @param n Length of the arrays.
/// @param quad_xy Labels of the x-only partition of (x, y): y as a function of x.
/// @param quad_yx Labels of the x-only partition of (y, x); unused when asym.
/// @param quad_joint Labels of the joint partition of (x, y); unused when asym.
/// @param asym Only the partition on x (the dependence of y on x).
/// @return DepResult containing the scalar correlation and dependence.
DepResult dep_pair(const double* x, const double* y, std::size_t n,
                   const std::vector<std::string>& quad_xy,
                   const std::vector<std::string>& quad_yx,
                   const std::vector<std::string>& quad_joint, bool asym = false);

/// Unfloored dependence of y on x at each distinct level of the partition on
/// x (R NNS_dep_levels_cpp); {0} when x or y is constant.
std::vector<double> dep_levels(const double* x, const double* y, std::size_t n,
                               const std::vector<std::string>& quad_xy);

/// Finest-cell dependence (R NNS_dep_local_cpp), the measure NNS.reg's default
/// order and dy.d_'s step rule read. quad_yx is unused when asym.
double dep_local(const double* x, const double* y, std::size_t n,
                 const std::vector<std::string>& quad_xy,
                 const std::vector<std::string>& quad_yx, bool asym = true);

/// NNS.dep(x, y, asym): dep_pair on the partitions R's native path builds.
/// Throws std::invalid_argument on non-finite input.
DepResult dep(const double* x, const double* y, std::size_t n, bool asym = false);

/// The pairwise dependence matrix of R NNS.dep(X, asym). Entry [i, j] is
/// dep(column i, column j); symmetric unless asym.
///
/// @param X Pointer to the column-major data matrix.
/// @param n Number of rows in X.
/// @param p Number of columns in X.
/// @param asym Calculate asymmetric dependence (true/false).
/// @param nthreads Number of parallel threads to use (-1 for hardware max).
/// @return DepMatrixResult containing the column-major correlation and dependence matrices.
DepMatrixResult dep_matrix(const double* X, std::size_t n, std::size_t p,
                           bool asym = false, int nthreads = -1);

} // namespace nns

#endif // NNS_DEPENDENCE_HPP
