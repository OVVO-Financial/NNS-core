// tests/fidelity/nns_eval.cpp
//
// Tiny CLI that evaluates NNS-core functions so that their outputs can be
// compared against live R (Rcpp) truth. It reads a function name and named
// arguments, computes the result with the portable C++ core, and prints it as
// a comma-separated list of full-precision numbers (a scalar is a list of
// length one). The fidelity orchestrator (scripts/run_fidelity.py) drives this
// binary once per truth case and compares the values with those produced by
// the installed R NNS package.
//
// Arguments are given either on the command line (`--key value`) or, with
// `--stdin`, one per line on standard input as `key<TAB>value` (no length
// limit; run_fidelity.py uses this form). Vectors are comma-separated;
// matrices are column-major vectors plus a `<name>_nrow` scalar.
//
// Usage examples:
//   nns_eval lpm --degree 1 --target 3 --x 1,2,3,4,5
//   nns_eval upm_v --degree 3 --target 0.1,0.2,... --x ...
//   nns_eval co_lpm --degree 1 --x 1,2,3 --y 4,5,6 --target_x 2 --target_y 5
//   nns_eval d_lpm --degree_lpm 1 --degree_upm 1 --x 1,2,3 --y 4,5,6 \
//            --target_x 2 --target_y 5
//
// Functions: lpm upm lpm_v upm_v lpm_ratio_v upm_ratio_v co_lpm co_upm d_lpm
// d_upm co_lpm_v co_upm_v d_lpm_v d_upm_v clpm_nd cupm_nd dpm_nd
// clpm_nd_batch pm_matrix gravity mode partition dep_matrix dep_pair distance_path
// distance_bulk fast_lm duplicate_column_map mreg_setup mreg_predict
// mreg_predict_path xstar_path.
//
// SPDX-License-Identifier: GPL-3.0-only
#include "nns/nns.hpp"

#include <cstdio>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

namespace {

using Flags = std::map<std::string, std::string>;

std::vector<double> parse_list(const std::string& s) {
  std::vector<double> out;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, ',')) {
    if (!item.empty()) out.push_back(std::strtod(item.c_str(), nullptr));
  }
  return out;
}

[[noreturn]] void die(const std::string& msg) {
  std::fprintf(stderr, "nns_eval: %s\n", msg.c_str());
  std::exit(2);
}

bool has(const Flags& f, const std::string& k) { return f.find(k) != f.end(); }

const std::string& need(const Flags& f, const std::string& key) {
  auto it = f.find(key);
  if (it == f.end()) die("missing required flag --" + key);
  return it->second;
}

double scalar(const Flags& f, const std::string& key) {
  return std::strtod(need(f, key).c_str(), nullptr);
}

double scalar_or(const Flags& f, const std::string& key, double dflt) {
  return has(f, key) ? scalar(f, key) : dflt;
}

std::vector<double> list(const Flags& f, const std::string& key) {
  return parse_list(need(f, key));
}

struct Matrix {
  std::vector<double> v;
  std::size_t nrow = 0, ncol = 0;
  const double* data() const { return v.data(); }
};

Matrix matrix(const Flags& f, const std::string& key) {
  Matrix m;
  m.v = list(f, key);
  m.nrow = static_cast<std::size_t>(scalar(f, key + "_nrow"));
  if (m.nrow == 0 || m.v.size() % m.nrow != 0) die("bad matrix shape for --" + key);
  m.ncol = m.v.size() / m.nrow;
  return m;
}

bool flag(const Flags& f, const std::string& key, bool dflt = false) {
  return has(f, key) ? scalar(f, key) != 0.0 : dflt;
}

// "q" -> 0, "q1243" -> 1243 (labels only contain the digits 1-4).
double quadrant_code(const std::string& q) {
  return q.size() > 1 ? std::strtod(q.c_str() + 1, nullptr) : 0.0;
}

void append(std::vector<double>& out, const std::vector<double>& v) {
  out.insert(out.end(), v.begin(), v.end());
}

std::vector<double> evaluate(const std::string& func, const Flags& f) {
  std::vector<double> out;

  if (func == "lpm" || func == "upm") {
    const std::vector<double> x = list(f, "x");
    const double d = scalar(f, "degree"), t = scalar(f, "target");
    out.push_back(func == "lpm" ? nns::lpm(d, t, x.data(), x.size())
                                : nns::upm(d, t, x.data(), x.size()));
  } else if (func == "lpm_v" || func == "upm_v" || func == "lpm_ratio_v" ||
             func == "upm_ratio_v") {
    const std::vector<double> x = list(f, "x"), t = list(f, "target");
    const double d = scalar(f, "degree");
    out.resize(t.size());
    if (func == "lpm_v") nns::lpm_v(d, t.data(), t.size(), x.data(), x.size(), out.data());
    if (func == "upm_v") nns::upm_v(d, t.data(), t.size(), x.data(), x.size(), out.data());
    if (func == "lpm_ratio_v")
      nns::lpm_ratio_v(d, t.data(), t.size(), x.data(), x.size(), out.data());
    if (func == "upm_ratio_v")
      nns::upm_ratio_v(d, t.data(), t.size(), x.data(), x.size(), out.data());
  } else if (func == "co_lpm" || func == "co_upm") {
    const std::vector<double> x = list(f, "x"), y = list(f, "y");
    const double dx = scalar(f, "degree");
    // R wrappers default degree_y to the single supplied degree.
    const double dy = scalar_or(f, "degree_y", dx);
    const double tx = scalar(f, "target_x"), ty = scalar(f, "target_y");
    out.push_back(func == "co_lpm"
                      ? nns::co_lpm(dx, dy, x.data(), y.data(), x.size(), y.size(), tx, ty)
                      : nns::co_upm(dx, dy, x.data(), y.data(), x.size(), y.size(), tx, ty));
  } else if (func == "d_lpm" || func == "d_upm") {
    const std::vector<double> x = list(f, "x"), y = list(f, "y");
    const double dl = scalar(f, "degree_lpm"), du = scalar(f, "degree_upm");
    const double tx = scalar(f, "target_x"), ty = scalar(f, "target_y");
    out.push_back(func == "d_lpm"
                      ? nns::d_lpm(dl, du, x.data(), y.data(), x.size(), y.size(), tx, ty)
                      : nns::d_upm(dl, du, x.data(), y.data(), x.size(), y.size(), tx, ty));
  } else if (func == "co_lpm_v" || func == "co_upm_v" || func == "d_lpm_v" ||
             func == "d_upm_v") {
    const std::vector<double> x = list(f, "x"), y = list(f, "y");
    const std::vector<double> tx = list(f, "target_x"), ty = list(f, "target_y");
    const double d1 = scalar(f, "degree_1"), d2 = scalar(f, "degree_2");
    out.resize(std::max(tx.size(), ty.size()));
    auto call = [&](auto fn) {
      fn(d1, d2, x.data(), y.data(), x.size(), y.size(), tx.data(), tx.size(), ty.data(),
         ty.size(), out.data(), -1);
    };
    if (func == "co_lpm_v") call(nns::co_lpm_v);
    if (func == "co_upm_v") call(nns::co_upm_v);
    if (func == "d_lpm_v") call(nns::d_lpm_v);
    if (func == "d_upm_v") call(nns::d_upm_v);
  } else if (func == "clpm_nd" || func == "cupm_nd" || func == "dpm_nd") {
    const Matrix X = matrix(f, "data");
    const std::vector<double> t = list(f, "target");
    const double d = scalar(f, "degree");
    const bool norm = flag(f, "norm", true);
    if (func == "clpm_nd") out.push_back(nns::clpm_nd(X.data(), X.nrow, X.ncol, t.data(), d, norm));
    if (func == "cupm_nd") out.push_back(nns::cupm_nd(X.data(), X.nrow, X.ncol, t.data(), d, norm));
    if (func == "dpm_nd") out.push_back(nns::dpm_nd(X.data(), X.nrow, X.ncol, t.data(), d, norm));
  } else if (func == "clpm_nd_batch") {
    const Matrix X = matrix(f, "data"), T = matrix(f, "targets");
    out.resize(T.nrow);
    nns::clpm_nd_batch(X.data(), X.nrow, X.ncol, T.data(), T.nrow, scalar(f, "degree"),
                       flag(f, "norm", true), out.data());
  } else if (func == "pm_matrix") {
    const Matrix X = matrix(f, "variable");
    const std::vector<double> t = list(f, "target");
    const nns::PMMatrixResult r =
        nns::pm_matrix(scalar(f, "degree_lpm"), scalar(f, "degree_upm"), t.data(), X.data(),
                       X.nrow, X.ncol, flag(f, "pop_adj"), flag(f, "norm"));
    append(out, r.cupm);
    append(out, r.dupm);
    append(out, r.dlpm);
    append(out, r.clpm);
    append(out, r.cov);
  } else if (func == "gravity") {
    const std::vector<double> x = list(f, "x");
    out.push_back(nns::gravity(x.data(), x.size(), flag(f, "discrete")));
  } else if (func == "mode") {
    const std::vector<double> x = list(f, "x");
    out = nns::mode(x.data(), x.size(), flag(f, "discrete"), flag(f, "multi"));
  } else if (func == "partition") {
    // Output: order, quadrant codes (one per row), then the regression points
    // as x values followed by y values.
    const std::vector<double> x = list(f, "x"), y = list(f, "y");
    const std::optional<std::string> type =
        flag(f, "xonly") ? std::optional<std::string>("XONLY") : std::nullopt;
    const std::optional<int> order =
        has(f, "order") ? std::optional<int>(static_cast<int>(scalar(f, "order"))) : std::nullopt;
    const nns::PartitionResult r =
        nns::partition(x.data(), y.data(), x.size(), type, order,
                       static_cast<int>(scalar_or(f, "obs_req", 8)), flag(f, "min_obs_stop"),
                       has(f, "noise_reduction") ? need(f, "noise_reduction") : "off");
    out.push_back(r.order);
    for (const auto& row : r.dt) out.push_back(quadrant_code(row.quadrant));
    for (const auto& rp : r.regression_points) out.push_back(rp.x);
    for (const auto& rp : r.regression_points) out.push_back(rp.y);
  } else if (func == "dep_matrix") {
    const Matrix X = matrix(f, "X");
    const nns::DepMatrixResult r = nns::dep_matrix(X.data(), X.nrow, X.ncol, flag(f, "asym"));
    append(out, r.correlation);
    append(out, r.dependence);
  } else if (func == "dep_pair") {
    // quad_xy / quad_yx: partition cell ids, compared only for equality.
    const std::vector<double> x = list(f, "x"), y = list(f, "y");
    const std::vector<double> qa = list(f, "quad_xy"), qb = list(f, "quad_yx");
    if (y.size() != x.size() || qa.size() != x.size() || qb.size() != x.size())
      die("dep_pair: x, y, quad_xy and quad_yx must have the same length");
    const std::vector<uint64_t> a(qa.begin(), qa.end()), b(qb.begin(), qb.end());
    const nns::DepResult r =
        nns::dep_pair(x.data(), y.data(), x.size(), a.data(), b.data(), flag(f, "asym"));
    out.push_back(r.correlation);
    out.push_back(r.dependence);
  } else if (func == "distance_path" || func == "distance_bulk") {
    const Matrix R = matrix(f, "RPM"), T = matrix(f, "Xtest");
    const std::vector<double> yhat = list(f, "yhat");
    const int k = static_cast<int>(scalar(f, "k"));
    out = func == "distance_path"
              ? nns::distance_path(R.data(), R.nrow, R.ncol, yhat.data(), T.data(), T.nrow, k,
                                   flag(f, "is_class"))
              : nns::distance_bulk(R.data(), R.nrow, R.ncol, yhat.data(), T.data(), T.nrow, k,
                                   flag(f, "is_class"));
  } else if (func == "fast_lm") {
    const std::vector<double> x = list(f, "x"), y = list(f, "y");
    const nns::FastLmResult r = nns::fast_lm(x.data(), y.data(), x.size());
    append(out, r.coef);
    append(out, r.fitted_values);
  } else if (func == "duplicate_column_map") {
    const Matrix X = matrix(f, "X");
    for (std::size_t r : nns::duplicate_column_map(X.data(), X.nrow, X.ncol))
      out.push_back(static_cast<double>(r + 1));  // 1-based like R
  } else if (func == "mreg_setup") {
    // Boundaries: concatenated per column, with their lengths in
    // boundary_lengths. Output: the RPM matrix (column-major).
    const Matrix X = matrix(f, "X");
    const std::vector<double> y = list(f, "y"), b = list(f, "boundaries"),
                              bl = list(f, "boundary_lengths");
    std::vector<std::vector<double>> bounds;
    std::size_t pos = 0;
    for (double len : bl) {
      const std::size_t l = static_cast<std::size_t>(len);
      bounds.emplace_back(b.begin() + pos, b.begin() + pos + l);
      pos += l;
    }
    const nns::MregSetupResult r =
        nns::mreg_setup(X.data(), X.nrow, X.ncol, y.data(), bounds,
                        static_cast<int>(scalar(f, "reducer")), flag(f, "is_class"));
    out = r.rpm;
  } else if (func == "mreg_predict" || func == "mreg_predict_path") {
    const Matrix R = matrix(f, "rpm_x"), T = matrix(f, "Xtest");
    const std::vector<double> yhat = list(f, "yhat"), mins = list(f, "mins"),
                              maxs = list(f, "maxs");
    const int k = static_cast<int>(scalar(f, "k"));
    const int dist = static_cast<int>(scalar(f, "dist_code"));
    const int threads = static_cast<int>(scalar_or(f, "n_threads", 1));
    out = func == "mreg_predict"
              ? nns::mreg_predict(R.data(), R.nrow, R.ncol, yhat.data(), T.data(), T.nrow, k,
                                  dist, mins.data(), maxs.data(), flag(f, "is_class"), threads)
              : nns::mreg_predict_path(R.data(), R.nrow, R.ncol, yhat.data(), T.data(),
                                       T.nrow, k, dist, mins.data(), maxs.data(),
                                       flag(f, "is_class"), threads);
  } else if (func == "xstar_path") {
    // column_order is 1-based as in R. Output: train, test, denominator.
    const Matrix A = matrix(f, "train"), B = matrix(f, "test");
    const std::vector<double> coef = list(f, "coef"), ord = list(f, "column_order");
    std::vector<std::size_t> order;
    for (double v : ord) order.push_back(static_cast<std::size_t>(v) - 1U);
    const nns::XstarPathResult r = nns::xstar_path(A.data(), A.nrow, B.data(), B.nrow,
                                                   A.ncol, coef.data(), order.data());
    append(out, r.train);
    append(out, r.test);
    append(out, r.denominator);
  } else {
    die("unknown function: " + func);
  }
  return out;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) die("usage: nns_eval <func> [--stdin | --flag value ...]");
  const std::string func = argv[1];

  Flags flags;
  bool from_stdin = false;
  for (int i = 2; i < argc;) {
    std::string key = argv[i];
    if (key == "--stdin") {
      from_stdin = true;
      ++i;
      continue;
    }
    if (key.rfind("--", 0) != 0 || i + 1 >= argc) die("expected --flag value, got: " + key);
    flags[key.substr(2)] = argv[i + 1];
    i += 2;
  }
  if (from_stdin) {
    std::string line;
    while (std::getline(std::cin, line)) {
      if (line.empty()) continue;
      const std::size_t tab = line.find('\t');
      if (tab == std::string::npos) die("stdin lines must be key<TAB>value");
      flags[line.substr(0, tab)] = line.substr(tab + 1);
    }
  }

  std::vector<double> result;
  try {
    result = evaluate(func, flags);
  } catch (const std::exception& e) {
    die(std::string("error: ") + e.what());
  }
  for (std::size_t i = 0; i < result.size(); ++i)
    std::printf(i ? ",%.17g" : "%.17g", result[i]);
  std::printf("\n");
  return 0;
}
