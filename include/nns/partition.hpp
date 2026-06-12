// include/nns/partition.hpp
//
// SPDX-License-Identifier: GPL-3.0-only
#ifndef NNS_PARTITION_HPP
#define NNS_PARTITION_HPP

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace nns {

struct PartitionRow {
  double x;
  double y;
  std::string quadrant;
  std::string prior_quadrant;
};

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

struct PartitionResult {
  int order = 0;
  bool quadrants_only = false;
  std::vector<std::string> quadrant;
  std::vector<PartitionRow> dt;
  std::vector<RegressionPoint> regression_points;
  std::vector<SegmentH> segments_h;
  std::vector<SegmentV> segments_v;
  std::vector<double> vlines;
};

/// Pure C++ port of original NNS_part_cpp.  x and y are observation vectors of
/// length n.  A present `type` optional enables the upstream x-only path; its
/// string contents are intentionally ignored.  Labels and output field names map
/// to original R payload names (`quadrant`, `prior.quadrant`, `segments_h`,
/// `segments_v`).
PartitionResult partition(const double* x,
                          const double* y,
                          std::size_t n,
                          const std::optional<std::string>& type,
                          const std::optional<int>& order_in,
                          int obs_req,
                          bool min_obs_stop,
                          const std::string& noise_reduction,
                          bool quadrants_only = false);

} // namespace nns

#endif // NNS_PARTITION_HPP
