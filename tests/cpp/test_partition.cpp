// tests/cpp/test_partition.cpp
//
// NNS.part treats +-Inf as values (R 13.2, complete_case_value = !ISNAN).
// Expected labels: R 13.3
//   x <- c(1:20, Inf, -Inf); y <- c(20:1, 5, 7)
//   NNS:::NNS_part_cpp(x, y, NULL, 2L, 8L, FALSE, "off")$dt$quadrant
//
// SPDX-License-Identifier: GPL-3.0-only
#include "nns/partition.hpp"

#include <limits>
#include <string>
#include <vector>

#include "check.hpp"

using nns_test::check;

int main() {
  std::vector<double> x, y;
  for (int i = 1; i <= 20; ++i) {
    x.push_back(i);
    y.push_back(21 - i);
  }
  x.push_back(std::numeric_limits<double>::infinity());
  y.push_back(5);
  x.push_back(-std::numeric_limits<double>::infinity());
  y.push_back(7);

  const nns::PartitionResult r =
      nns::partition(x.data(), y.data(), x.size(), std::nullopt, 2, 8, false, "off");
  const std::vector<std::string> want = {
      "q22", "q22", "q22", "q22", "q21", "q23", "q23", "q23", "q23", "q1", "q1",
      "q32", "q32", "q32", "q32", "q31", "q33", "q33", "q33", "q33", "q31", "q4"};
  check(r.dt.size() == want.size(), "row count");
  for (std::size_t i = 0; i < want.size() && i < r.dt.size(); ++i)
    check(r.dt[i].quadrant == want[i], "quadrant row " + std::to_string(i + 1) + ": got " +
                                           r.dt[i].quadrant + ", want " + want[i]);
  return nns_test::check_summary("test_partition");
}
