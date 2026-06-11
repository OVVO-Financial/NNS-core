// include/nns/partition.hpp
//
// SPDX-License-Identifier: GPL-3.0-only
#ifndef NNS_PARTITION_HPP
#define NNS_PARTITION_HPP

#include <cstddef>
#include <string>
#include <vector>

namespace nns {

// --- Output Data Structures ---

struct RegressionPoint {
  std::string quadrant;
  double x;
  double y;
};

struct SegmentH {
  double x0;
  double x1;
  double y;
};

struct SegmentV {
  double x;
  double y0;
  double y1;
};

struct PartitionNode {
  int id;
  double x;
  double y;
  std::string quadrant;
};

struct PartitionResult {
  std::vector<RegressionPoint> regression_points;
  std::vector<SegmentH> segment_h;
  std::vector<SegmentV> segment_v;
  std::vector<PartitionNode> nodes;
};

// --- Core API ---

/// Iteratively partition an (x, y) space based on central gravity.
///
/// @param x Pointer to the predictor array.
/// @param y Pointer to the response array.
/// @param n Number of observations.
/// @param Voronoi Flag to enable fine-grained Voronoi splitting (1 = yes, 0 = no).
/// @param min_obs Minimum number of observations required to split a node.
/// @param type Aggregation rule for y ("mean", "median", "sum", "last").
/// @return A PartitionResult struct containing quadrants, regression points, and boundaries.
PartitionResult partition(const double* x, const double* y, std::size_t n, 
                          int Voronoi, int min_obs, const std::string& type = "mean");

} // namespace nns

#endif // NNS_PARTITION_HPP