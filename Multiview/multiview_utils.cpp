//multiview_utils.cpp
#include <Rcpp.h>
#include <algorithm>        // std::remove
#include <unordered_map>    // std::unordered_map
#include <cmath>            // std::log, std::exp
using namespace Rcpp;

#include "multiview_utils.h"
#include "multiview_state.h"


double compute_f_vk(int v, int k, int i);
double compute_f_vk_new(int v, int i);

void compute_table_probs_with_cache(
    int i,
    std::vector<double> &prob_existing,
    double &prob_new
) {
  std::vector<std::unordered_map<int, double>> cache_fvk(d);
  
  for (int t = 0; t < T; ++t) {
    if (n_t[t] == 0) {
      prob_existing[t] = 0.0;
      continue; 
    }
    
    double log_prob_t = 0.0;
    
    for (int v = 0; v < d; ++v) {
      int k = dish_of[v][t];  
      
      auto &cache_v = cache_fvk[v];
      
      double f_vk;
      auto it = cache_v.find(k);
      if (it == cache_v.end()) {
        f_vk = compute_f_vk(v, k, i);
        cache_v[k] = f_vk;
      } else {
        f_vk = it->second;
      }
      
      log_prob_t += std::log(f_vk);
    }
    
    double mass_t = n_t[t] - sigma_global; 
    if (mass_t <= 0.0) {
      prob_existing[t] = 0.0;
    } else {
      prob_existing[t] = mass_t * std::exp(log_prob_t);
    }
  }
  
  double log_prob_new = 0.0;
  for (int v = 0; v < d; ++v) {
    double f_vk_new = compute_f_vk_new(v, i);
    log_prob_new += std::log(f_vk_new);
  }
  
  int K_nonempty = 0;
  for (int t = 0; t < T; ++t) {
    if (n_t[t] > 0) ++K_nonempty;
  }
  
  double mass_new = alpha_global + sigma_global * K_nonempty;
  
  if (mass_new <= 0.0) {
    prob_new = 0.0;
  } else {
    prob_new = mass_new * std::exp(log_prob_new);
  }
}

void remove_customer(int i) {
  int t = table_of[i];  
  
  if (t < 0 || t >= T) {
    Rcpp::stop("remove_customer: invalid table index for customer");
  }
  
  auto &list_t = customers_at_table[t];
  auto it = std::find(list_t.begin(), list_t.end(), i);
  if (it == list_t.end()) {
    Rcpp::stop("remove_customer: customer not found in customers_at_table[t]");
  }
  std::swap(*it, list_t.back());
  list_t.pop_back();
  n_t[t]--;
  
  for (int v = 0; v < d; ++v) {
    int k = dish_of[v][t];              // dish della view v al tavolo t
    if (k < 0 || k >= (int)views[v].n_vk.size()) {
      Rcpp::stop("remove_customer: invalid dish index");
    }
    
    ViewState &V = views[v];
    
    V.n_vk[k]--;
    V.sum_y[k]  -= y[v][i];
    V.sum_y2[k] -= y[v][i] * y[v][i];
    
    auto &list_k = V.customers_at_dish[k];
    auto it2 = std::find(list_k.begin(), list_k.end(), i);
    if (it2 == list_k.end()) {
      Rcpp::stop("remove_customer: customer not found in customers_at_dish[k]");
    }
    std::swap(*it2, list_k.back());
    list_k.pop_back();

    if (V.n_vk[k] == 0 && V.l_vk[k] > 0) {
      V.l_vk[k]--;
    }
  }
  
  table_of[i] = -1;
  
  if (n_t[t] == 0) {
    int last = T - 1;
    if (t != last) {
      for (int v = 0; v < d; ++v) {
        int k_old = dish_of[v][t];
        if (k_old >= 0 && k_old < (int)views[v].l_vk.size()
              && views[v].l_vk[k_old] > 0) {
          views[v].l_vk[k_old]--;
        }
      }
      
      customers_at_table[t] = std::move(customers_at_table[last]);
      n_t[t] = n_t[last];
      for (int v = 0; v < d; ++v) {
        dish_of[v][t] = dish_of[v][last];
      }
      
      for (int j : customers_at_table[t]) {
        table_of[j] = t;
      }
      
      for (int v = 0; v < d; ++v) {
        int k_new = dish_of[v][t];
        if (k_new >= 0 && k_new < (int)views[v].l_vk.size()) {
          views[v].l_vk[k_new]++;
        }
      }
    }
    
    customers_at_table.pop_back();
    n_t.pop_back();
    for (int v = 0; v < d; ++v) {
      dish_of[v].pop_back();  
    }
    T--;
  }
}

void add_customer_to_existing_table(int i, int t) {
  if (t < 0 || t >= T) {
    Rcpp::stop("add_customer_to_existing_table: invalid table index");
  }
  
  table_of[i] = t;                     
  customers_at_table[t].push_back(i);
  n_t[t]++;
  
  for (int v = 0; v < d; ++v) {
    int k = dish_of[v][t];
    if (k < 0 || k >= (int)views[v].n_vk.size()) {
      Rcpp::stop("add_customer_to_existing_table: invalid dish index");
    }
    
    ViewState &V = views[v];
    V.n_vk[k]++;
    V.sum_y[k]  += y[v][i];
    V.sum_y2[k] += y[v][i] * y[v][i];
    V.customers_at_dish[k].push_back(i);
    table_of[i] = t;
  }
}

int create_empty_table() {
  int t_new = T;           
  T++;
  
  n_t.push_back(0);
  customers_at_table.emplace_back();
  for (int v = 0; v < d; ++v)
    dish_of[v].push_back(-1);
  
  return t_new;
}

void add_customer_to_new_table(int i, int t_new) {
  if (t_new < 0 || t_new >= T) {
    Rcpp::stop("add_customer_to_new_table: invalid new table index");
  }
  table_of[i] = t_new;
  customers_at_table[t_new].push_back(i);
  n_t[t_new] = 1;
}

int sample_dish_for_new_table(int v, int i) {
  if (v < 0 || v >= d) {
    Rcpp::stop("sample_dish_for_new_table: invalid view index");
  }
  ViewState &V = views[v];
  
  std::vector<double> weights;
  std::vector<int> candidate_dishes;
  
  // 1. EXISTING dishes
  for (int k = 0; k < V.K; ++k) {
    if (V.l_vk[k] > 0) {
      double f_vk = compute_f_vk(v, k, i);
      double weight = (V.l_vk[k] - V.sigma_v) * f_vk;
      if (weight < 0.0) weight = 0.0;
      weights.push_back(weight);
      candidate_dishes.push_back(k);
    }
  }
  
  // 2. NEW dish
  double f_vk_new = compute_f_vk_new(v, i);
  double weight_new = (V.alpha_v + V.sigma_v * V.K) * f_vk_new;
  if (weight_new < 0.0) weight_new = 0.0;
  weights.push_back(weight_new);
  
  double total_weight = 0.0;
  for (double w : weights) total_weight += w;
  if (total_weight <= 0.0) {
    Rcpp::stop("sample_dish_for_new_table: total_weight <= 0");
  }
  
  double u = R::runif(0.0, total_weight);
  double cumulative = 0.0;
  int selected_dish = -1;
  
  for (size_t j = 0; j < candidate_dishes.size(); ++j) {
    cumulative += weights[j];
    if (u < cumulative) {
      selected_dish = candidate_dishes[j];
      break;
    }
  }
  
  if (selected_dish == -1) {
    // nuovo dish
    int new_k = V.K;
    V.K++;
    V.n_vk.push_back(0);
    V.l_vk.push_back(0);
    V.sum_y.push_back(0.0);
    V.sum_y2.push_back(0.0);
    V.customers_at_dish.push_back({});
    return new_k;
  }
  
  return selected_dish;
}

void assign_dishes_new_table(int i, int t_new) {
  if (t_new < 0 || t_new >= T) {
    Rcpp::stop("assign_dishes_new_table: invalid table index");
  }
  
  for (int v = 0; v < d; ++v) {
    int k = sample_dish_for_new_table(v, i);
    
    if (t_new >= (int)dish_of[v].size()) {
      Rcpp::stop("assign_dishes_new_table: dish_of[v] too small");
    }
    
    dish_of[v][t_new] = k;
    
    ViewState &V = views[v];
    if (k < 0 || k >= (int)V.n_vk.size()) {
      Rcpp::stop("assign_dishes_new_table: invalid k after sampling");
    }
    
    V.l_vk[k]++;
    V.n_vk[k]++;
    V.sum_y[k]  += y[v][i];
    V.sum_y2[k] += y[v][i] * y[v][i];
    V.customers_at_dish[k].push_back(i);
  }
}

void save_state() {
  saved_table_of.push_back(table_of);
  saved_dish_of.push_back(dish_of);
}

double uniform01() {
  return R::runif(0.0, 1.0);
}

double rnorm_scalar(double mean, double sd) {
  return R::rnorm(mean, sd);
}


double compute_f_vk(int v, int k, int i) {
  const ViewState &V = views[v];
  
  double y_vi = y[v][i];      
  double tau  = V.tau_v;      
  
  if (tau <= 0.0) {
    tau = 1.0; 
  }
  
  int n_k = 0;
  double mean_k = 0.0;
  
  if (k >= 0 && k < (int)V.n_vk.size()) {
    n_k = V.n_vk[k];
    if (n_k > 0) {
      mean_k = V.sum_y[k] / static_cast<double>(n_k);
    }
  }
  
  double diff   = y_vi - mean_k;
  double var    = tau;
  double logden = -0.5 * std::log(2.0 * M_PI * var)
    -0.5 * (diff * diff) / var;
  
  return std::exp(logden);
}

double compute_f_vk_new(int v, int i) {
  const ViewState &V = views[v];
  
  double y_vi = y[v][i];
  double tau  = V.tau_v;
  
  if (tau <= 0.0) {
    tau = 1.0;
  }
  
  double mean0 = 0.0;   
  double diff  = y_vi - mean0;
  double var   = tau;
  double logden = -0.5 * std::log(2.0 * M_PI * var)
    -0.5 * (diff * diff) / var;
  
  return std::exp(logden);
}