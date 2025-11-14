#include <Rcpp.h>
using namespace Rcpp;

double compute_f_vk(int v, int k, int i);
double compute_f_vk_new(int v, int i);

void compute_table_probs_with_cache(
    int i,
    std::vector<double> &prob_existing,
    double &prob_new
) {
  std::vector<std::unordered_map<int, double>> cache_fvk(d);
  
  for (int t = 0; t < T; ++t) {
    double log_prob_t = 0.0;
    for (int v = 0; v < d; ++v) {
      int k = dish_of[v][t];
      ViewState &V = views[v];
      
      double f_vk;
      auto it = cache_fvk[v].find(k);
      if (it == cache_fvk[v].end()) {
        f_vk = compute_f_vk(v, k, i);
        cache_fvk[v][k] = f_vk;
      } else {
        f_vk = it->second;
      }
      log_prob_t += std::log(f_vk);
    }
    double mass_t = n_t[t];
    prob_existing[t] = mass_t * std::exp(log_prob_t);
  }
  
  double log_prob_new = 0.0;
  for (int v = 0; v < d; ++v) {
    double f_vk_new = compute_f_vk_new(v, i);
    log_prob_new += std::log(f_vk_new);
  }
  double mass_new = 1.0;
  prob_new = mass_new * std::exp(log_prob_new);
}


void remove_customer(int i) { ... }

void add_customer_to_existing_table(int i, int t) { ... }

void create_empty_table(int t_new) { ... }

void add_customer_to_new_table(int i, int t_new) { ... }

int sample_dish_for_new_table(int v, int i);

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

