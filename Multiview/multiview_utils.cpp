// multiview_utils.cpp
#include <Rcpp.h>
#include <algorithm>        // std::remove
#include <unordered_map>    // std::unordered_map
#include <cmath>            // std::log, std::exp
#include <vector>
using namespace Rcpp;

#include "multiview_utils.h"
#include "multiview_state.h"

std::vector<double> global_view_variance;

void ensure_global_variances_calculated() {
  if ((int)global_view_variance.size() != d) {
    global_view_variance.assign(d, 10.0);
    for(int v=0; v<d; ++v) {
      if (y[v].empty()) continue;
      double sum = 0, sum2 = 0;
      for(double val : y[v]) {
        sum += val;
        sum2 += val*val;
      }
      double n_curr = (double)y[v].size();
      double mean = sum / n_curr;
      double var = (sum2 - sum*sum/n_curr) / (n_curr - 1.0);
      if(var < 1e-4) var = 1.0; 
      global_view_variance[v] = var;
    }
  }
}

double compute_f_vk(int v, int k, int i);
double compute_f_vk_new(int v, int i);

double compute_marginal_likelihood_new_table(int v, int i) {
  const ViewState &V = views[v];
  
  double total_tables_v = 0.0;
  for (int count : V.l_vk) total_tables_v += count;
  
  double denominator = V.alpha_v + total_tables_v;
  if (denominator <= 0.0) return compute_f_vk_new(v, i);
  
  double likelihood_sum = 0.0;
  
  for (int k = 0; k < V.K; ++k) {
    if (V.l_vk[k] > 0) {
      double weight = (V.l_vk[k] - V.sigma_v);
      
      likelihood_sum += weight * compute_f_vk(v, k, i);
    }
  }
  
  int K_active = 0;
  for(int c : V.l_vk) if(c > 0) K_active++;
  
  double weight_new = (V.alpha_v + K_active * V.sigma_v);
  
  likelihood_sum += weight_new * compute_f_vk_new(v, i);
  
  return likelihood_sum / denominator;
}


void compute_table_probs_with_cache(
    int i,
    std::vector<double> &prob_existing,
    double &prob_new,
    std::vector<std::unordered_map<int, double>> &cache_fvk 
) {
  for(auto &map : cache_fvk) {
    map.clear();
  }
  
  if(global_view_variance.empty()) ensure_global_variances_calculated();
  
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
      if (it != cache_v.end()) {
        f_vk = it->second;
      } else {

        f_vk = compute_f_vk(v, k, i);
        cache_v[k] = f_vk;
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
  
  double log_prob_new_table_data = 0.0;
  for (int v = 0; v < d; ++v) {
    double marg_lik = compute_marginal_likelihood_new_table(v, i);
    log_prob_new_table_data += std::log(marg_lik);
  }
  
  int T_nonempty = 0;
  for (int t = 0; t < T; ++t) {
    if (n_t[t] > 0) ++T_nonempty;
  }
  
  double mass_new = alpha_global + sigma_global * T_nonempty;
  
  if (mass_new <= 0.0) {
    prob_new = 0.0;
  } else {
    prob_new = mass_new * std::exp(log_prob_new_table_data);
  }
}


void remove_customer(int i) {
  int t = table_of[i];  
  
  if (t < 0 || t >= T) Rcpp::stop("remove_customer: invalid table");
  
  auto &list_t = customers_at_table[t];
  auto it = std::find(list_t.begin(), list_t.end(), i);
  if (it == list_t.end()) Rcpp::stop("customer not at table");
  
  std::swap(*it, list_t.back());
  list_t.pop_back();
  n_t[t]--;
  
  for (int v = 0; v < d; ++v) {
    int k = dish_of[v][t];              
    ViewState &V = views[v];
    
    V.n_vk[k]--;
    V.sum_y[k]  -= y[v][i];
    V.sum_y2[k] -= y[v][i] * y[v][i];
    
    auto &list_k = V.customers_at_dish[k];
    auto it2 = std::find(list_k.begin(), list_k.end(), i);
    std::swap(*it2, list_k.back());
    list_k.pop_back();
    
  }
  
  table_of[i] = -1;
  
  if (n_t[t] == 0) {

    for (int v = 0; v < d; ++v) {
      int k = dish_of[v][t];
      if(k >= 0 && views[v].l_vk[k] > 0) views[v].l_vk[k]--;
    }
    
    int last = T - 1;
    if (t != last) {
      customers_at_table[t] = std::move(customers_at_table[last]);
      n_t[t] = n_t[last];
      for (int v = 0; v < d; ++v) {
        dish_of[v][t] = dish_of[v][last];
      }
      for (int j : customers_at_table[t]) {
        table_of[j] = t;
      }
    }
    
    customers_at_table.pop_back();
    n_t.pop_back();
    for (int v = 0; v < d; ++v) dish_of[v].pop_back();  
    T--;
  }
}

void add_customer_to_existing_table(int i, int t) {
  table_of[i] = t;                    
  customers_at_table[t].push_back(i);
  n_t[t]++;
  
  for (int v = 0; v < d; ++v) {
    int k = dish_of[v][t];
    ViewState &V = views[v];
    V.n_vk[k]++;
    V.sum_y[k]  += y[v][i];
    V.sum_y2[k] += y[v][i] * y[v][i];
    V.customers_at_dish[k].push_back(i);
  }
}

int create_empty_table() {
  int t_new = T;          
  T++;
  n_t.push_back(0);
  customers_at_table.emplace_back();
  for (int v = 0; v < d; ++v) dish_of[v].push_back(-1);
  return t_new;
}

void add_customer_to_new_table(int i, int t_new) {
  table_of[i] = t_new;
  customers_at_table[t_new].push_back(i);
  n_t[t_new] = 1;
}



int sample_dish_for_new_table(int v, int i) {
  ViewState &V = views[v];
  std::vector<double> weights;
  std::vector<int> candidate_dishes;
  
  for (int k = 0; k < V.K; ++k) {
    if (V.l_vk[k] > 0) {
      double f_vk = compute_f_vk(v, k, i);
      double w = (V.l_vk[k] - V.sigma_v) * f_vk;
      if(w<0) w=0;
      weights.push_back(w);
      candidate_dishes.push_back(k);
    }
  }
  
  int K_active = 0; 
  for(int c : V.l_vk) if(c>0) K_active++;
  
  double f_vk_new = compute_f_vk_new(v, i);
  double w_new = (V.alpha_v + V.sigma_v * K_active) * f_vk_new;
  if(w_new<0) w_new=0;
  
  weights.push_back(w_new);
  
  double total_w = 0; 
  for(double w:weights) total_w+=w;
  
  if(total_w <= 0) {
    return (V.K); 
  }
  
  double u = R::runif(0.0, total_w);
  double cum = 0;
  for(size_t j=0; j<candidate_dishes.size(); ++j) {
    cum += weights[j];
    if(u < cum) return candidate_dishes[j];
  }
  
  int new_k = V.K;
  V.K++;
  V.n_vk.push_back(0);
  V.l_vk.push_back(0);
  V.sum_y.push_back(0.0);
  V.sum_y2.push_back(0.0);
  V.customers_at_dish.push_back({});
  return new_k;
}

void assign_dishes_new_table(int i, int t_new) {
  for (int v = 0; v < d; ++v) {
    int k = sample_dish_for_new_table(v, i);
    dish_of[v][t_new] = k;
    ViewState &V = views[v];
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

double uniform01() { return R::runif(0.0, 1.0); }
double rnorm_scalar(double mean, double sd) { return R::rnorm(mean, sd); }



double compute_f_vk(int v, int k, int i) {
  const ViewState &V = views[v];
  double y_vi = y[v][i];      
  double tau  = V.tau_v;      
  if (tau <= 1e-9) tau = 1e-9; 

  double mean_k = 0.0;
  if (V.n_vk[k] > 0) {
    mean_k = V.sum_y[k] / static_cast<double>(V.n_vk[k]);
  }
  
  double diff = y_vi - mean_k;
  double logden = -0.5 * std::log(2.0 * M_PI * tau) -0.5 * (diff * diff) / tau;
  return std::exp(logden);
}


double compute_f_vk_new(int v, int i) {
  if(global_view_variance.empty()) ensure_global_variances_calculated();
  
  const ViewState &V = views[v];
  double y_vi = y[v][i];
  double tau  = V.tau_v;
  
  double pred_var = tau + global_view_variance[v];
  
  double diff = y_vi - 0.0; 
  
  double logden = -0.5 * std::log(2.0 * M_PI * pred_var) -0.5 * (diff * diff) / pred_var;
  return std::exp(logden);
}