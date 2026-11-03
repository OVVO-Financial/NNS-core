// tests/cpp/test_partition.cpp
//
// NNS.part treats +-Inf as values (R 13.2, complete_case_value = !ISNAN), and
// its "mode" reducer is the continuous mode.
// Expected labels: R 499f82b1
//   x <- c(1:20, Inf, -Inf); y <- c(20:1, 5, 7)
//   NNS:::NNS_part_cpp(x, y, NULL, 2L, 8L, FALSE, "off")$dt$quadrant
//
// SPDX-License-Identifier: GPL-3.0-only
#include "nns/partition.hpp"

#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "check.hpp"

using nns_test::check;
using nns_test::check_close;

namespace {

// noise_reduction = "mode" reduces each cell with the continuous mode; the
// discrete one rounded data in [0, 1] to integers. Expected: R 499f82b1
//   x <- (0:199) / 199; NNS:::NNS_part_cpp(x, sin(3 * x), NULL, NULL, 8L, FALSE, "mode")
void mode_reducer_case() {
  std::vector<double> x, y;
  for (int i = 0; i < 200; ++i) {
    x.push_back(i / 199.0);
    y.push_back(std::sin(3.0 * x.back()));
  }
  const nns::PartitionResult r =
      nns::partition(x.data(), y.data(), x.size(), std::nullopt, std::nullopt, 8, false, "mode");
  const std::vector<double> want_x = {
      0.97738693467336679, 0.92964824120603018, 0.90452261306532666, 0.86180904522613067,
      0.80653266331658291, 0.73618090452261309, 0.66080402010050254, 0.59045226130653261,
      0.53517587939698497, 0.44221105527638194, 0.40954773869346733, 0.3592964824120603,
      0.31407035175879394, 0.26633165829145727, 0.23366834170854273, 0.17839195979899497,
      0.14572864321608039, 0.097989949748743727, 0.055276381909547742, 0.027638190954773871};
  const std::vector<double> want_y = {
      0.2226196579024336, 0.31691197620799005, 0.38742999225541136, 0.50207658265289423,
      0.63757768365444356, 0.80341674315749834, 0.91031428944862003, 0.98864842480038284,
      0.99969933928628785, 0.99079815904141599, 0.9420076846615717, 0.88093605337862879,
      0.80883720099370371, 0.71660123129169018, 0.64496772445162398, 0.50997784708592475,
      0.42335967337573188, 0.28974581417958145, 0.1650576542399845, 0.082817248874378677};
  check(r.order == 7, "mode reducer order " + std::to_string(r.order));
  check(r.regression_points.size() == want_x.size(),
        "mode reducer point count " + std::to_string(r.regression_points.size()));
  for (std::size_t i = 0; i < want_x.size() && i < r.regression_points.size(); ++i) {
    check_close("mode reducer x " + std::to_string(i), r.regression_points[i].x, want_x[i], 1e-13);
    check_close("mode reducer y " + std::to_string(i), r.regression_points[i].y, want_y[i], 1e-13);
  }
}

}  // namespace

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
      "q22", "q22", "q22", "q22", "q22", "q23", "q23", "q23", "q23", "q23", "q1",
      "q32", "q32", "q32", "q32", "q34", "q33", "q33", "q33", "q33", "q33", "q4"};
  check(r.dt.size() == want.size(), "row count");
  for (std::size_t i = 0; i < want.size() && i < r.dt.size(); ++i)
    check(r.dt[i].quadrant == want[i], "quadrant row " + std::to_string(i + 1) + ": got " +
                                           r.dt[i].quadrant + ", want " + want[i]);
  mode_reducer_case();
  return nns_test::check_summary("test_partition");
}
