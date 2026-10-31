// tests/cpp/check.hpp
//
// Minimal dependency-free check helpers shared by the C++ tests. Each test
// binary calls the helpers and returns check_summary() from main(), so ctest
// sees a non-zero exit code on any failure.
//
// SPDX-License-Identifier: GPL-3.0-only
#ifndef NNS_TESTS_CHECK_HPP
#define NNS_TESTS_CHECK_HPP

#include <cmath>
#include <cstdio>
#include <exception>
#include <string>

namespace nns_test {

inline int& failures() {
  static int n = 0;
  return n;
}

inline int& checks() {
  static int n = 0;
  return n;
}

inline void fail(const std::string& label, const std::string& detail) {
  ++failures();
  std::fprintf(stderr, "FAIL %s: %s\n", label.c_str(), detail.c_str());
}

inline void check(bool ok, const std::string& label) {
  ++checks();
  if (!ok) fail(label, "condition is false");
}

// |got - want| <= abs_tol + rel_tol * |want|; NaN only matches NaN.
inline void check_close(const std::string& label, double got, double want,
                        double rel_tol, double abs_tol = 0.0) {
  ++checks();
  if (std::isnan(got) || std::isnan(want)) {
    if (!(std::isnan(got) && std::isnan(want))) {
      char buf[128];
      std::snprintf(buf, sizeof buf, "got %.17g, want %.17g", got, want);
      fail(label, buf);
    }
    return;
  }
  if (got == want) return;
  const double diff = std::fabs(got - want);
  if (!(diff <= abs_tol + rel_tol * std::fabs(want))) {
    char buf[160];
    std::snprintf(buf, sizeof buf, "got %.17g, want %.17g (diff %.3g)", got, want, diff);
    fail(label, buf);
  }
}

template <typename Fn>
void check_throws(const std::string& label, Fn&& fn) {
  ++checks();
  try {
    fn();
  } catch (const std::exception&) {
    return;
  }
  fail(label, "expected an exception");
}

inline int check_summary(const char* name) {
  if (failures() != 0) {
    std::fprintf(stderr, "%s: %d of %d check(s) failed\n", name, failures(), checks());
    return 1;
  }
  std::fprintf(stderr, "%s: all %d checks passed\n", name, checks());
  return 0;
}

}  // namespace nns_test

#endif  // NNS_TESTS_CHECK_HPP
