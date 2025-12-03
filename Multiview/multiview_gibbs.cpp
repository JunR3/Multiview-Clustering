// multiview_gibbs.cpp
#include <Rcpp.h>
#include <numeric> 
#include <cmath>
#include "multiview_state.h"
#include "multiview_gibbs.h"
#include "multiview_hyper.h"
#include "multiview_utils.h"

using namespace Rcpp;

static void initialize_state_from_data() {
  
  // Do we really want fixed initializations ?
  int K_init_tables = 4; 
  int K_init_dishes = 2; 
  
  
  T = K_init_tables;
  table_of.assign(n, 0);
  n_t.assign(T, 0);
  
  customers_at_table.clear();
  customers_at_table.resize(T);
  
  for (int i = 0; i < n; ++i) {
    int t = static_cast<int>(std::floor(R::runif(0.0, (double)T)));
    if (t < 0) t = 0;
    if (t >= T) t = T - 1;
    
    table_of[i] = t;
    customers_at_table[t].push_back(i);
    n_t[t]++;
  }
  
  dish_of.clear();
  dish_of.resize(d);
  for (int v = 0; v < d; ++v) {
    dish_of[v].assign(T, 0);
  }
  
  views.clear();
  views.resize(d);
  
  for (int v = 0; v < d; ++v) {
    ViewState &V = views[v];
    
    V.K = K_init_dishes; 
    V.n_vk.assign(V.K, 0);
    V.l_vk.assign(V.K, 0);
    V.sum_y.assign(V.K, 0.0);
    V.sum_y2.assign(V.K, 0.0);
    V.customers_at_dish.clear();
    V.customers_at_dish.resize(V.K);
    
    for (int t = 0; t < T; ++t) {
      int k = static_cast<int>(std::floor(R::runif(0.0, (double)V.K)));
      if (k < 0) k = 0;
      if (k >= V.K) k = V.K - 1;
      
      dish_of[v][t] = k;
      V.l_vk[k]++; 
    }
    
    for (int i = 0; i < n; ++i) {
      int t = table_of[i];
      int k = dish_of[v][t];
      
      double val = y[v][i];
      V.n_vk[k]++;
      V.sum_y[k] += val;
      V.sum_y2[k] += val * val;
      V.customers_at_dish[k].push_back(i);
    }
    
    V.alpha_v = 1.0;
    V.sigma_v = 0.5;
    
    double s1 = 0.0;
    for (int i = 0; i < n; ++i) {
      s1 += y[v][i];
    }
    double mean = s1 / std::max(1, n);
    double var = 0.0;
    if (n > 1) {
      for (int i = 0; i < n; ++i) {
        double diff = y[v][i] - mean;
        var += diff * diff;
      }
      var /= (n - 1);
    } else {
      var = 1.0;
    }
    if (var <= 0.0) var = 1.0;
    V.tau_v = var*0.25 * 0.01; // modified so the values of tau_v are closer to their final values
  }
  
  alpha_global = 1;
  sigma_global = 0.6;
  
  saved_table_of.clear();
  saved_dish_of.clear();
  saved_loglik.clear();
}

// [[Rcpp::export]]
Rcpp::List run_gibbs_cpp(const Rcpp::List& data_views,
                         int M, int burn_in, int thin) {
  
  d = data_views.size();
  n = Rcpp::as<Rcpp::NumericVector>(data_views[0]).size();
  
  y.clear();
  y.resize(d);
  for (int v = 0; v < d; ++v)
    y[v] = Rcpp::as<std::vector<double>>(data_views[v]);
  
  initialize_state_from_data();
  
  gibbs_sampler(M, burn_in, thin);
  
  return Rcpp::List::create(
    Rcpp::Named("table_of") = saved_table_of,
    Rcpp::Named("dish_of") = saved_dish_of,
    Rcpp::Named("loglik")  = saved_loglik,
    Rcpp::Named("alpha_v") = saved_alpha_v,
    Rcpp::Named("sigma_v") = saved_sigma_v,
    Rcpp::Named("tau_v")   = saved_tau_v,
    Rcpp::Named("alpha_global") = saved_alpha_global,
    Rcpp::Named("sigma_global") = saved_sigma_global
  );
}
// multiview_gibbs.cpp
// ... (include e init invariati) ...

void gibbs_sampler(int M, int burn_in, int thin) {
  
  saved_table_of.clear();
  saved_dish_of.clear();
  saved_loglik.clear();
  saved_alpha_v.assign(d, {});
  saved_sigma_v.assign(d, {});
  saved_tau_v.assign(d, {});
  saved_alpha_global.clear();
  saved_sigma_global.clear();
  
  
  std::vector<std::unordered_map<int, double>> workspace_cache(d);
  
  for(auto &m : workspace_cache) m.reserve(64);
  
  for (int iter = 0; iter < M; ++iter) {
    
    if ( (iter + 1) % 100 == 0 ) {
      Rcpp::Rcout << "Iteration " << (iter + 1)
                  << " / " << M << std::endl;
    }
    
    for (int i = 0; i < n; ++i) {
      
      remove_customer(i);
      
      
      std::vector<double> prob_existing(T, 0.0);
      double prob_new = 0.0;
      
      
      compute_table_probs_with_cache(i, prob_existing, prob_new, workspace_cache);
      
      
      double sum_p = prob_new;
      for (int t = 0; t < T; ++t) sum_p += prob_existing[t];
      
      if (sum_p <= 0.0) {
        
        add_customer_to_existing_table(i, 0);
        continue;
      }
      
      for (int t = 0; t < T; ++t) prob_existing[t] /= sum_p;
      prob_new /= sum_p;
      
      double u = uniform01();
      double cum = 0.0;
      int t_star = -1;
      
      for (int t = 0; t < T; ++t) {
        cum += prob_existing[t];
        if (u < cum) {
          t_star = t;
          break;
        }
      }
      
      if (t_star == -1) {
        int t_new = create_empty_table();
        add_customer_to_new_table(i, t_new);
        assign_dishes_new_table(i, t_new);
      } else {
        add_customer_to_existing_table(i, t_star);
      }
    }
    
    update_hyperparameters();
    
    // Save state after burn-in and thinning
    /*if (iter >= burn_in && ((iter - burn_in) % thin == 0)) {
      save_state();
    }*/
   if (iter % thin == 0) {
    save_state(); // Save one in every "thin" iteration for now
   }
    
  }
}
