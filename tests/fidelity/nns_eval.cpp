// tests/fidelity/nns_eval.cpp
//
// Tiny CLI that evaluates NNS-core partial moment functions so that their
// outputs can be compared against live R (Rcpp) truth. It reads a function name
// and flag arguments, computes a single scalar with the portable C++ core, and
// prints it with full precision. The fidelity orchestrator (scripts/run_fidelity.py)
// drives this binary once per truth case and compares the printed value to the
// value produced by the installed R NNS package.
//
// Usage examples:
//   nns_eval lpm --degree 1 --target 3 --x 1,2,3,4,5
//   nns_eval co_lpm --degree 1 --x 1,2,3 --y 4,5,6 --target_x 2 --target_y 5
//   nns_eval d_lpm --degree_lpm 1 --degree_upm 1 --x 1,2,3 --y 4,5,6 \
//            --target_x 2 --target_y 5
//
// SPDX-License-Identifier: GPL-3.0-only
#include "nns/partial_moments.hpp"

#include <cstdio>
#include <cstdlib>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::vector<double> parse_list(const std::string& s) {
  std::vector<double> out;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, ',')) {
    if (!item.empty()) {
      out.push_back(std::strtod(item.c_str(), nullptr));
    }
  }
  return out;
}

[[noreturn]] void die(const std::string& msg) {
  std::fprintf(stderr, "nns_eval: %s\n", msg.c_str());
  std::exit(2);
}

double need_scalar(const std::map<std::string, std::string>& flags,
                   const std::string& key) {
  auto it = flags.find(key);
  if (it == flags.end()) {
    die("missing required flag --" + key);
  }
  return std::strtod(it->second.c_str(), nullptr);
}

std::vector<double> need_list(const std::map<std::string, std::string>& flags,
                              const std::string& key) {
  auto it = flags.find(key);
  if (it == flags.end()) {
    die("missing required flag --" + key);
  }
  return parse_list(it->second);
}

} // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    die("usage: nns_eval <func> [--flag value ...]");
  }
  const std::string func = argv[1];

  std::map<std::string, std::string> flags;
  for (int i = 2; i + 1 < argc; i += 2) {
    std::string key = argv[i];
    if (key.rfind("--", 0) != 0) {
      die("expected --flag, got: " + key);
    }
    flags[key.substr(2)] = argv[i + 1];
  }

  double result = 0.0;

  if (func == "lpm" || func == "upm") {
    const double degree = need_scalar(flags, "degree");
    const double target = need_scalar(flags, "target");
    const std::vector<double> x = need_list(flags, "x");
    result = (func == "lpm")
                 ? nns::lpm(degree, target, x.data(), x.size())
                 : nns::upm(degree, target, x.data(), x.size());
  } else if (func == "co_lpm" || func == "co_upm") {
    const double degree = need_scalar(flags, "degree");
    const std::vector<double> x = need_list(flags, "x");
    const std::vector<double> y = need_list(flags, "y");
    const double tx = need_scalar(flags, "target_x");
    const double ty = need_scalar(flags, "target_y");
    // R wrappers default degree_y to the single supplied degree.
    result = (func == "co_lpm")
                 ? nns::co_lpm(degree, degree, x.data(), y.data(), x.size(),
                               y.size(), tx, ty)
                 : nns::co_upm(degree, degree, x.data(), y.data(), x.size(),
                               y.size(), tx, ty);
  } else if (func == "d_lpm" || func == "d_upm") {
    const double dl = need_scalar(flags, "degree_lpm");
    const double du = need_scalar(flags, "degree_upm");
    const std::vector<double> x = need_list(flags, "x");
    const std::vector<double> y = need_list(flags, "y");
    const double tx = need_scalar(flags, "target_x");
    const double ty = need_scalar(flags, "target_y");
    result = (func == "d_lpm")
                 ? nns::d_lpm(dl, du, x.data(), y.data(), x.size(), y.size(),
                              tx, ty)
                 : nns::d_upm(dl, du, x.data(), y.data(), x.size(), y.size(),
                              tx, ty);
  } else {
    die("unknown function: " + func);
  }

  std::printf("%.17g\n", result);
  return 0;
}
