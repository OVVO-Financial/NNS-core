// src/partition.cpp
//
// Implementation extracted from NNS 13.0. Decoupled from Rcpp.
//
// SPDX-License-Identifier: GPL-3.0-only
#include "nns/partition.hpp"
#include "nns/central_tendencies.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace nns {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// ---------- Aggregation Helpers ----------

inline double mean_no_na(const std::vector<double>& v) {
  long double s = 0.0L; 
  std::size_t m = 0;
  for (double xi : v) {
    if (std::isfinite(xi)) { s += xi; ++m; }
  }
  return m ? static_cast<double>(s / m) : kNaN;
}

inline double median_no_na(const std::vector<double>& v) {
  std::vector<double> a; 
  a.reserve(v.size());
  for (double xi : v) {
    if (std::isfinite(xi)) a.push_back(xi);
  }
  if (a.empty()) return kNaN;
  std::size_t n = a.size();
  std::nth_element(a.begin(), a.begin() + n / 2, a.end());
  double hi = a[n / 2];
  if (n % 2 == 1) return hi;
  auto lm = std::max_element(a.begin(), a.begin() + n / 2);
  return (*lm + hi) * 0.5;
}

inline double sum_no_na(const std::vector<double>& v) {
  long double s = 0.0L; 
  std::size_t m = 0;
  for (double xi : v) {
    if (std::isfinite(xi)) { s += xi; ++m; }
  }
  return m ? static_cast<double>(s) : kNaN;
}

struct Agg {
  std::string type;
  Agg(std::string t) : type(t) {}
  
  double for_x(const std::vector<double>& v) {
    if (type == "median") return median_no_na(v);
    return mean_no_na(v);
  }
  
  double for_y(const std::vector<double>& v) {
    if (type == "sum") return sum_no_na(v);
    if (type == "median") return median_no_na(v);
    if (type == "last") {
      for (auto it = v.rbegin(); it != v.rend(); ++it) {
        if (std::isfinite(*it)) return *it;
      }
      return kNaN;
    }
    return mean_no_na(v);
  }
};

// ---------- Point Tracking Struct ----------

struct Point {
  int id;
  double x;
  double y;
  std::string quad;
};

} // namespace

// ---------- NNS Partition ----------

PartitionResult partition(const double* x, const double* y, std::size_t n, 
                          int Voronoi, int min_obs, const std::string& type) {
  if (n == 0) return PartitionResult{};
  
  std::vector<Point> df;
  df.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    df.push_back({static_cast<int>(i + 1), x[i], y[i], "1"});
  }
  
  std::vector<double> H_x0, H_x1, H_y;
  std::vector<double> V_x, V_y0, V_y1;
  
  bool any_split = true;
  while (any_split) {
    any_split = false;
    std::unordered_map<std::string, std::vector<int>> grouped;
    for (std::size_t i = 0; i < df.size(); ++i) grouped[df[i].quad].push_back(i);
    
    std::vector<Point> new_df = df;
    
    for (const auto& kv : grouped) {
      const std::string& q = kv.first;
      const std::vector<int>& idx = kv.second;
      
      if (static_cast<int>(idx.size()) <= min_obs) continue;
      
      std::vector<double> cx, cy;
      cx.reserve(idx.size());
      cy.reserve(idx.size());
      
      double x_min = std::numeric_limits<double>::infinity(), x_max = -std::numeric_limits<double>::infinity();
      double y_min = std::numeric_limits<double>::infinity(), y_max = -std::numeric_limits<double>::infinity();
      
      for (int i : idx) {
        double xi = df[i].x;
        double yi = df[i].y;
        cx.push_back(xi);
        cy.push_back(yi);
        if (std::isfinite(xi)) {
          if (xi < x_min) x_min = xi;
          if (xi > x_max) x_max = xi;
        }
        if (std::isfinite(yi)) {
          if (yi < y_min) y_min = yi;
          if (yi > y_max) y_max = yi;
        }
      }
      
      if (x_min == x_max && y_min == y_max) continue;
      
      double x_center = nns::gravity(cx.data(), cx.size(), false);
      double y_center = nns::gravity(cy.data(), cy.size(), false);
      
      if (!std::isfinite(x_center) || !std::isfinite(y_center)) continue;
      
      int tl = 0, tr = 0, bl = 0, br = 0;
      for (std::size_t k = 0; k < idx.size(); ++k) {
        double xi = cx[k], yi = cy[k];
        if (xi <= x_center && yi >= y_center) tl++;
        else if (xi > x_center && yi >= y_center) tr++;
        else if (xi <= x_center && yi < y_center) bl++;
        else if (xi > x_center && yi < y_center) br++;
      }
      
      bool will_split = false;
      if (Voronoi == 1) {
        int non_zero = (tl > 0) + (tr > 0) + (bl > 0) + (br > 0);
        will_split = (non_zero > 1);
      } else {
        int non_zero = (tl > 0) + (tr > 0) + (bl > 0) + (br > 0);
        bool sum_gt = ((tl + tr + bl + br) > min_obs);
        will_split = (non_zero > 1 && sum_gt);
      }
      
      if (will_split) {
        any_split = true;
        for (std::size_t k = 0; k < idx.size(); ++k) {
          int i = idx[k];
          double xi = df[i].x;
          double yi = df[i].y;
          std::string ext;
          if (xi <= x_center && yi >= y_center) ext = "1";
          else if (xi > x_center && yi >= y_center) ext = "2";
          else if (xi <= x_center && yi < y_center) ext = "3";
          else ext = "4";
          
          new_df[i].quad = q + ext;
        }
        
        H_x0.push_back(x_min); H_x1.push_back(x_max); H_y.push_back(y_center);
        V_x.push_back(x_center); V_y0.push_back(y_min); V_y1.push_back(y_max);
      }
    }
    df = new_df;
  }
  
  std::map<std::string, std::vector<int>> by_prior;
  for (std::size_t i = 0; i < df.size(); ++i) {
    by_prior[df[i].quad].push_back(i);
  }
  
  Agg agg(type);
  PartitionResult res;
  res.regression_points.reserve(by_prior.size());
  
  for (const auto& kv : by_prior) {
    const auto& idx = kv.second;
    std::vector<double> xv(idx.size()), yv(idx.size());
    for (std::size_t k = 0; k < idx.size(); ++k) {
      xv[k] = df[idx[k]].x;
      yv[k] = df[idx[k]].y;
    }
    res.regression_points.push_back({
      kv.first,
      agg.for_x(xv),
      agg.for_y(yv)
    });
  }
  
  res.segment_h.reserve(H_x0.size());
  for (std::size_t i = 0; i < H_x0.size(); ++i) {
    res.segment_h.push_back({H_x0[i], H_x1[i], H_y[i]});
  }
  
  res.segment_v.reserve(V_x.size());
  for (std::size_t i = 0; i < V_x.size(); ++i) {
    res.segment_v.push_back({V_x[i], V_y0[i], V_y1[i]});
  }
  
  res.nodes.reserve(df.size());
  for (const auto& row : df) {
    res.nodes.push_back({row.id, row.x, row.y, row.quad});
  }
  
  return res;
}

} // namespace nns