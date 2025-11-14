#ifndef MULTIVIEW_STATE_H
#define MULTIVIEW_STATE_H

#include <vector>

// Per-view state
struct ViewState {
  int K;
  std::vector<int> n_vk;
  std::vector<int> l_vk;
  std::vector<double> sum_y;
  std::vector<double> sum_y2;
  std::vector<std::vector<int>> customers_at_dish;
  double alpha_v;
  double sigma_v;
  double tau_v;
};

// Global state variables (dichiarazioni)
extern int n, d;
extern std::vector<std::vector<double>> y;

extern int T;
extern std::vector<int> table_of;
extern std::vector<int> n_t;
extern std::vector<std::vector<int>> customers_at_table;
extern std::vector<std::vector<int>> dish_of;

extern std::vector<ViewState> views;

// Saved states
extern std::vector<std::vector<int>> saved_table_of;
extern std::vector<std::vector<std::vector<int>>> saved_dish_of;
extern std::vector<double> saved_loglik;

#endif