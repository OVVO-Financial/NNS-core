// tests/cpp/test_dependence.cpp
//
// dep_pair on single-cell inputs where the copula's degree-0 counts decide
// the result. Expected values: NNS:::NNS_dep_pair_cpp from R 13.3 with the
// DEP-01, DEP-08 and CORE-09 patches.
//
// SPDX-License-Identifier: GPL-3.0-only
#include "nns/dependence.hpp"

#include <cstdint>
#include <string>
#include <vector>

#include "check.hpp"

using nns_test::check_close;

namespace {

void check_dep(const std::string& name, const std::vector<double>& x,
               const std::vector<double>& y, const std::vector<uint64_t>& cells,
               double correlation, double dependence) {
  for (bool asym : {false, true}) {
    const nns::DepResult r = nns::dep_pair(x.data(), y.data(), x.size(), cells.data(),
                                           cells.data(), asym);
    const std::string tag = name + (asym ? " asym" : "");
    check_close(tag + " correlation", r.correlation, correlation, 1e-12);
    check_close(tag + " dependence", r.dependence, dependence, 1e-12);
  }
}

}  // namespace

int main() {
  // Every point discordant (DEP-01): R 13.3 returned Correlation +1.
  std::vector<double> x(100), neg(100);
  for (int i = 0; i < 100; ++i) {
    x[i] = i + 1;
    neg[i] = -(i + 1);
  }
  check_dep("x vs -x", x, neg, std::vector<uint64_t>(100, 1), -1.0, 1.0);

  // Binary self-dependence in cells of 51 and 49 points (DEP-08): 49 * (1/49)
  // is not 1, so R 13.3 returned 0.51 / 0.755.
  std::vector<double> b(100, 0.0);
  std::vector<uint64_t> bcells(100, 1);
  for (int i = 51; i < 100; ++i) {
    b[i] = 1.0;
    bcells[i] = 2;
  }
  check_dep("binary self", b, b, bcells, 1.0, 1.0);

  // Points on the cell mean (CORE-09): R 13.3 counted them as co-lower and
  // discordant, giving -0.68718427093627676 / 0.67692546880147175.
  const double pattern[9] = {0, 1, -1, 0, -1, 1, 1, 0, -1};
  std::vector<double> tx(99), ty(99);
  for (int i = 0; i < 99; ++i) {
    tx[i] = (i % 3) - 1.0;
    ty[i] = pattern[i % 9];
  }
  check_dep("ties at the mean", tx, ty, std::vector<uint64_t>(99, 1), -0.55277079839256671,
            0.44305206586294998);

  return nns_test::check_summary("test_dependence");
}
