//multiview_state.cpp
#include "multiview_state.h"

int n = 0, d = 0;
std::vector<std::vector<double>> y;

double alpha_global = 1.0;
double sigma_global = 0.5;

int T = 0;
std::vector<int> table_of;
std::vector<int> n_t;
std::vector<std::vector<int>> customers_at_table;
std::vector<std::vector<int>> dish_of;

std::vector<ViewState> views;

std::vector<std::vector<int>> saved_table_of;
std::vector<std::vector<std::vector<int>>> saved_dish_of;
std::vector<double> saved_loglik;