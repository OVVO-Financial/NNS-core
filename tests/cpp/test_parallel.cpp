// tests/cpp/test_parallel.cpp
//
// nns::parallel_for honours the global thread setting and the serial work
// threshold, processes every index exactly once, and kernel results do not
// depend on the thread count.
//
// SPDX-License-Identifier: GPL-3.0-only
#include "nns/nns.hpp"

#include <atomic>
#include <mutex>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "check.hpp"

using nns_test::check;

namespace {

std::size_t distinct_threads(std::size_t n, int n_threads, std::size_t work) {
  std::mutex mu;
  std::set<std::thread::id> ids;
  std::vector<int> hits(n, 0);
  nns::parallel_for(
      0, n,
      [&](std::size_t b, std::size_t e) {
        for (std::size_t i = b; i < e; ++i) ++hits[i];
        std::lock_guard<std::mutex> lock(mu);
        ids.insert(std::this_thread::get_id());
      },
      n_threads, work);
  bool once = true;
  for (int h : hits) once = once && h == 1;
  check(once, "every index processed exactly once");
  return ids.size();
}

}  // namespace

int main() {
  check(nns::get_num_threads() >= 1, "automatic thread count >= 1");
  nns::set_num_threads(3);
  check(nns::get_num_threads() == 3, "set_num_threads(3)");
  check(distinct_threads(1000, -1, 1000) == 3, "default uses the global setting");
  check(distinct_threads(1000, 2, 1000) == 2, "explicit n_threads wins");
  check(distinct_threads(2, 8, 1000000) == 2, "never more threads than items");
  check(distinct_threads(1000, 4, 10) == 1, "small work runs serially");
  check(distinct_threads(1000, 4, 0) == 4, "unknown work keeps the requested threads");
  nns::set_num_threads(1);
  check(distinct_threads(100000, -1, 1000) == 1, "set_num_threads(1) disables threading");
  nns::set_num_threads(0);
  check(nns::get_num_threads() >= 1, "set_num_threads(0) restores the automatic count");
  distinct_threads(0, 4, 1);  // empty range is a no-op

  // A throwing chunk reaches the caller instead of std::terminate, whether a
  // spawned worker (index 0) or the calling thread (index 999) runs it.
  for (std::size_t bad : {std::size_t{0}, std::size_t{999}}) {
    std::atomic<int> chunks{0};
    bool caught = false;
    try {
      nns::parallel_for(
          0, 1000,
          [&](std::size_t b, std::size_t e) {
            ++chunks;
            if (b <= bad && bad < e) throw std::runtime_error("chunk failed");
          },
          4, 0);
    } catch (const std::runtime_error& ex) {
      caught = std::string(ex.what()) == "chunk failed";
    }
    check(caught, "worker exception rethrown on the caller, index " + std::to_string(bad));
    check(chunks == 4, "every chunk ran before the rethrow");
  }
  check(distinct_threads(1000, 4, 0) == 4, "parallel_for works after an exception");

  // Results are identical for any thread count.
  std::mt19937_64 g(3);
  std::normal_distribution<double> nd;
  const std::size_t n = 3000, d = 6;
  std::vector<double> X(n * d), t(500), tgt(d, 0.0);
  for (double& v : X) v = nd(g);
  for (double& v : t) v = nd(g);
  std::vector<double> ref(t.size()), out(t.size());
  nns::lpm_v(2.5, t.data(), t.size(), X.data(), n, ref.data(), 1);
  const nns::PMMatrixResult pm1 = nns::pm_matrix(1, 2, tgt.data(), X.data(), n, d, true, false, 1);
  const double c1 = nns::clpm_nd(X.data(), n, d, tgt.data(), 1.5, true, 1);
  for (int th : {2, 4, 7, -1}) {
    nns::lpm_v(2.5, t.data(), t.size(), X.data(), n, out.data(), th);
    check(out == ref, "lpm_v independent of thread count");
    const nns::PMMatrixResult pm = nns::pm_matrix(1, 2, tgt.data(), X.data(), n, d, true, false, th);
    check(pm.cov == pm1.cov && pm.dupm == pm1.dupm, "pm_matrix independent of thread count");
    check(nns::clpm_nd(X.data(), n, d, tgt.data(), 1.5, true, th) == c1,
          "clpm_nd independent of thread count");
  }
  return nns_test::check_summary("test_parallel");
}
