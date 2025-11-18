#include <Rcpp.h>
using namespace Rcpp;
#include "multiview_clustering.h"

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
    double mass_t = n_t[t]-sigma_global;
    prob_existing[t] = mass_t * std::exp(log_prob_t);
  }
  
  double log_prob_new = 0.0;
  for (int v = 0; v < d; ++v) {
    double f_vk_new = compute_f_vk_new(v, i);
    log_prob_new += std::log(f_vk_new);
  }
  double mass_new = alpha_global + sigma_global*T;
  prob_new = mass_new * std::exp(log_prob_new);
}


void remove_customer(int i) { 
  int t = table_of[i];
  auto &list_t = customers_at_table[t]}

void add_customer_to_existing_table(int i, int t) {
  customers_at_table[t].push_back(i);
  n_t++;
  
  for(int v = 0; v < d; v++){
    int k = dish_of[v][t];
    ViewState &V = views[v];
    V.n_vk[k]++;
    V.sum_y[k] += y[v][i];
    V.sum_y2[k] += y[v][i]*y[v][i];
    V.customer_at_dish[k].push_back(i);
  }
}

void create_empty_table(int t_new) { 
  T++;
  n_t.push_back(0);
  customers_at_table.emplace_back();
  for(int v = 0; v < d; v++)
    dish_of[v].push_back(-1);
}

void add_customer_to_new_table(int i, int t_new) {
  table_of[i]=t_new;
  customers_at_table[t_new].push_back(i);
  n_t[t_new]=1;
}

int sample_dish_for_new_table(int v, int i) {
    // Reference to the current view to simplify notation
    ViewState &V = views[v];

    std::vector<double> weights;
    std::vector<int> candidate_dishes; // Maps weight index -> actual dish index k

    // ---------------------------------------------------------
    // 1. Compute weights for EXISTING DISHES (k = old)
    // Formula: (l_vk - sigma) * Likelihood(y_vi | theta_k)
    // ---------------------------------------------------------
    for (int k = 0; k < V.K; ++k) {
        // We only consider dishes currently being served by at least one table (l_vk > 0).
        // Note: l_vk counts TABLES, n_vk counts CUSTOMERS.
        // Since we are assigning a *table* to a dish, the CRF dictates we use l_vk.
        if (V.l_vk[k] > 0) {
            double f_vk = compute_f_vk(v, k, i); // Likelihood helper

            // Pitman-Yor weight for an existing dish
            double weight = (V.l_vk[k] - V.sigma_v) * f_vk;

            weights.push_back(weight);
            candidate_dishes.push_back(k);
        }
    }

    // ---------------------------------------------------------
    // 2. Compute weight for a NEW DISH (k = new)
    // Formula: (alpha + sigma * K) * Marginal_Likelihood(y_vi)
    // ---------------------------------------------------------
    double f_vk_new = compute_f_vk_new(v, i); // Marginal likelihood helper

    // Pitman-Yor weight for a new dish
    // V.K is the current number of active dishes in this view
    double weight_new = (V.alpha_v + V.sigma_v * V.K) * f_vk_new;

    weights.push_back(weight_new);
    // We do not add anything to candidate_dishes for the new option;
    // we will identify it if the sampling falls into the last weight bucket.

    // ---------------------------------------------------------
    // 3. Sampling (Weighted Random Selection)
    // ---------------------------------------------------------
    double total_weight = 0.0;
    for (double w : weights) total_weight += w;

    std::uniform_real_distribution<double> dist(0.0, total_weight);
    double u = dist(rng);

    double current_sum = 0.0;
    int selected_dish = -1;

    // Iterate over existing dish options
    for (size_t j = 0; j < candidate_dishes.size(); ++j) {
        current_sum += weights[j];
        if (u < current_sum) {
            selected_dish = candidate_dishes[j];
            break;
        }
    }

    // ---------------------------------------------------------
    // 4. Result Handling & State Expansion
    // ---------------------------------------------------------

    // If selected_dish remains -1, it means u fell into the last weight range (New Dish)
    if (selected_dish == -1) {
        // --- CREATE NEW DISH ---
        int new_k = V.K; // The new index is the current count K

        // Update View structure (expand vectors to accommodate the new index)
        V.K++;
        V.n_vk.push_back(0);         // Initialize with 0 customers
        V.l_vk.push_back(0);         // Initialize with 0 tables (will be incremented by the caller)
        V.sum_y.push_back(0.0);      // Initialize empty sufficient statistics
        V.sum_y2.push_back(0.0);
        V.customers_at_dish.push_back({}); // Initialize empty customer list

        return new_k;
    }

    // If an existing dish was selected, return its index
    return selected_dish;
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
  saved_loglik.push_back(compute_log_likelihood());
}
