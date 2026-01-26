// multiview_merge_split.cpp
// Merge-Split Metropolis-Hastings moves for improved MCMC mixing

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include "multiview_merge_split.h"
#include "multiview_state.h"
#include "multiview_utils.h"

using namespace Rcpp;

namespace {

// Compute log marginal likelihood for a set of customers in view v
double compute_log_marginal_for_customers(int v,
                                          const std::vector<int> &customers) {
  if (customers.empty())
    return 0.0;

  const ViewState &V = views[v];
  double tau = V.tau_v;
  int n_k = customers.size();

  double sum_y = 0.0;
  double sum_y2 = 0.0;
  for (int i : customers) {
    double val = y[v][i];
    sum_y += val;
    sum_y2 += val * val;
  }

  // Log marginal likelihood under Normal-InverseGamma prior
  // p(y_1, ..., y_n) integrated over mean with variance tau
  double log_det =
      -0.5 * n_k * std::log(2.0 * M_PI * tau) - 0.5 * std::log(tau + n_k);
  double term1 = -0.5 * sum_y2 / tau;
  double term2 = 0.5 * (sum_y * sum_y) / (tau * (tau + n_k));

  return log_det + term1 + term2;
}

// Compute log CRF prior ratio for table configurations
double compute_log_crp_ratio_merge(int n_t1, int n_t2) {
  // Merge: going from 2 tables to 1 table
  // Prior ratio = P(1 table with n1+n2) / P(2 tables with n1, n2)
  // Under CRP, this involves the discount parameter

  int n_merged = n_t1 + n_t2;

  // Log Pochhammer symbol (rising factorial) for cluster sizes
  double log_prior_merged = 0.0;
  for (int m = 1; m < n_merged; ++m) {
    log_prior_merged += std::log(m - sigma_global);
  }

  double log_prior_separate = 0.0;
  for (int m = 1; m < n_t1; ++m) {
    log_prior_separate += std::log(m - sigma_global);
  }
  for (int m = 1; m < n_t2; ++m) {
    log_prior_separate += std::log(m - sigma_global);
  }

  // Table creation term: one less table after merge
  double log_table_term = -std::log(alpha_global + (T - 1) * sigma_global);

  return log_prior_merged - log_prior_separate + log_table_term;
}

double compute_log_crp_ratio_split(int n_t1, int n_t2) {
  // Split: going from 1 table to 2 tables (inverse of merge)
  return -compute_log_crp_ratio_merge(n_t1, n_t2);
}

// Propose merge move
bool propose_merge() {
  if (T < 2)
    return false; // Need at least 2 tables to merge

  // Sample two distinct tables uniformly
  int t1 = static_cast<int>(std::floor(R::runif(0.0, T)));
  int t2 = static_cast<int>(std::floor(R::runif(0.0, T - 1)));
  if (t2 >= t1)
    t2++;

  if (t1 < 0 || t1 >= T || t2 < 0 || t2 >= T)
    return false;
  if (n_t[t1] == 0 || n_t[t2] == 0)
    return false;

  // Get customers at each table
  std::vector<int> cust1 = customers_at_table[t1];
  std::vector<int> cust2 = customers_at_table[t2];
  std::vector<int> cust_merged;
  cust_merged.reserve(cust1.size() + cust2.size());
  cust_merged.insert(cust_merged.end(), cust1.begin(), cust1.end());
  cust_merged.insert(cust_merged.end(), cust2.begin(), cust2.end());

  // Compute log likelihood ratio for each view
  double log_lik_ratio = 0.0;
  for (int v = 0; v < d; ++v) {
    double log_lik_merged = compute_log_marginal_for_customers(v, cust_merged);
    double log_lik_t1 = compute_log_marginal_for_customers(v, cust1);
    double log_lik_t2 = compute_log_marginal_for_customers(v, cust2);

    log_lik_ratio += log_lik_merged - (log_lik_t1 + log_lik_t2);
  }

  // Compute prior ratio
  double log_prior_ratio = compute_log_crp_ratio_merge(n_t[t1], n_t[t2]);

  // Proposal ratio:
  // Forward: choose 2 tables from T tables = T*(T-1)/2 ways
  // Reverse: choose 1 table from (T-1) tables, then split customers
  // Reverse split has 2^(n-1) - 1 possible splits (excluding empty)
  int n_merged = cust_merged.size();
  double log_proposal_ratio = std::log(T * (T - 1) / 2.0) - std::log(T - 1) -
                              (n_merged - 1) * std::log(2.0);

  double log_accept = log_lik_ratio + log_prior_ratio + log_proposal_ratio;

  if (std::log(R::runif(0.0, 1.0)) < log_accept) {
    // Accept: merge t2 into t1
    // Move all customers from t2 to t1
    for (int i : cust2) {
      // Update dish stats for all views
      for (int v = 0; v < d; ++v) {
        int k_old = dish_of[v][t2];
        int k_new = dish_of[v][t1];
        ViewState &V = views[v];

        // Remove from old dish
        V.n_vk[k_old]--;
        V.sum_y[k_old] -= y[v][i];
        V.sum_y2[k_old] -= y[v][i] * y[v][i];
        auto &list_old = V.customers_at_dish[k_old];
        auto it = std::find(list_old.begin(), list_old.end(), i);
        if (it != list_old.end()) {
          std::swap(*it, list_old.back());
          list_old.pop_back();
        }

        // Add to new dish
        V.n_vk[k_new]++;
        V.sum_y[k_new] += y[v][i];
        V.sum_y2[k_new] += y[v][i] * y[v][i];
        V.customers_at_dish[k_new].push_back(i);
      }

      table_of[i] = t1;
      customers_at_table[t1].push_back(i);
      n_t[t1]++;
    }

    // Decrement dish counts for t2
    for (int v = 0; v < d; ++v) {
      int k = dish_of[v][t2];
      if (k >= 0 && views[v].l_vk[k] > 0) {
        views[v].l_vk[k]--;
      }
    }

    // Remove t2 by swapping with last table
    n_t[t2] = 0;
    customers_at_table[t2].clear();

    int last = T - 1;
    if (t2 != last) {
      customers_at_table[t2] = std::move(customers_at_table[last]);
      n_t[t2] = n_t[last];
      for (int v = 0; v < d; ++v) {
        dish_of[v][t2] = dish_of[v][last];
      }
      for (int j : customers_at_table[t2]) {
        table_of[j] = t2;
      }
    }

    customers_at_table.pop_back();
    n_t.pop_back();
    for (int v = 0; v < d; ++v)
      dish_of[v].pop_back();
    T--;

    return true;
  }

  return false;
}

// Propose split move
bool propose_split() {
  if (T == 0)
    return false;

  // Sample a table with at least 2 customers
  std::vector<int> eligible_tables;
  for (int t = 0; t < T; ++t) {
    if (n_t[t] >= 2) {
      eligible_tables.push_back(t);
    }
  }

  if (eligible_tables.empty())
    return false;

  int idx = static_cast<int>(std::floor(R::runif(0.0, eligible_tables.size())));
  if (idx >= (int)eligible_tables.size())
    idx = eligible_tables.size() - 1;
  int t_split = eligible_tables[idx];

  std::vector<int> customers = customers_at_table[t_split];
  int n_orig = customers.size();

  // Randomly assign customers to two groups
  std::vector<int> group1, group2;
  for (int i : customers) {
    if (R::runif(0.0, 1.0) < 0.5) {
      group1.push_back(i);
    } else {
      group2.push_back(i);
    }
  }

  // Ensure neither group is empty
  if (group1.empty()) {
    group1.push_back(group2.back());
    group2.pop_back();
  }
  if (group2.empty()) {
    group2.push_back(group1.back());
    group1.pop_back();
  }

  // Compute log likelihood ratio
  double log_lik_ratio = 0.0;
  for (int v = 0; v < d; ++v) {
    double log_lik_orig = compute_log_marginal_for_customers(v, customers);
    double log_lik_g1 = compute_log_marginal_for_customers(v, group1);
    double log_lik_g2 = compute_log_marginal_for_customers(v, group2);

    log_lik_ratio += (log_lik_g1 + log_lik_g2) - log_lik_orig;
  }

  // Compute prior ratio
  double log_prior_ratio =
      compute_log_crp_ratio_split(group1.size(), group2.size());

  // Proposal ratio (inverse of merge)
  double log_proposal_ratio = std::log(eligible_tables.size()) +
                              (n_orig - 1) * std::log(2.0) -
                              std::log((T + 1) * T / 2.0);

  double log_accept = log_lik_ratio + log_prior_ratio + log_proposal_ratio;

  if (std::log(R::runif(0.0, 1.0)) < log_accept) {
    // Accept: split into two tables
    // Keep group1 at t_split, create new table for group2

    // Create new table
    int t_new = T;
    T++;
    n_t.push_back(0);
    customers_at_table.emplace_back();
    for (int v = 0; v < d; ++v) {
      // New table gets same dish as original (simplification)
      dish_of[v].push_back(dish_of[v][t_split]);
      views[v].l_vk[dish_of[v][t_split]]++;
    }

    // Move group2 customers to new table
    for (int i : group2) {
      // Update table assignment
      table_of[i] = t_new;

      // Remove from old table's customer list
      auto &list_old = customers_at_table[t_split];
      auto it = std::find(list_old.begin(), list_old.end(), i);
      if (it != list_old.end()) {
        std::swap(*it, list_old.back());
        list_old.pop_back();
      }
      n_t[t_split]--;

      // Add to new table
      customers_at_table[t_new].push_back(i);
      n_t[t_new]++;
    }

    return true;
  }

  return false;
}

} // anonymous namespace

void merge_split_step() {
  // Alternate between merge and split proposals
  // Do a few attempts per iteration for better mixing

  for (int attempt = 0; attempt < 3; ++attempt) {
    if (R::runif(0.0, 1.0) < 0.5) {
      propose_merge();
    } else {
      propose_split();
    }
  }
}
