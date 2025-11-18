#include <Rcpp.h>
#include "multiview_state.h"
#include "multiview_gibbs.h"
#include "multiview_hyper.h"
#include "multiview_utils.h"

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List run_gibbs_cpp(const Rcpp::List& data_views,
                         int M, int burn_in, int thin) {
  
  d = data_views.size();
  n = Rcpp::as<Rcpp::NumericVector>(data_views[0]).size();
  
  y.resize(d);
  for (int v = 0; v < d; ++v)
    y[v] = Rcpp::as<std::vector<double>>(data_views[v]);
  
  gibbs_sampler(M, burn_in, thin);
  
  return Rcpp::List::create(
    Rcpp::Named("table_of") = saved_table_of,
    Rcpp::Named("dish_of")  = saved_dish_of,
    Rcpp::Named("loglik")   = saved_loglik
  );
}

void gibbs_sampler(int M, int burn_in, int thin) {
  
  saved_table_of.clear();
  saved_dish_of.clear();
  saved_loglik.clear();
  
  for (int iter = 0; iter < M; ++iter) {
    for(int i = 0; i < n; ++i){
      // Step 1: remove customer i
      remove_customer(i);
      
      // Step 2: compute table probabilities (existing/new)
      std::vector<double> prob_existing(T);
      double prob_new = 0.0;
      compute_table_probs_with_cache(i, prob_existing, prob_new);
      
      // normalise
      double sum_p = prob_new;
      for(int t = 0; t <T; T++) sum_p += prob_existing[t];
      for(int t = 0; t <T; T++) prob_existing[t]/=sum_p;
      prob_new/= sum_p;
      
      double u = uniform01();
      double cum = 0.0;
      int t_star = -1;
      for(int t = 0; t <T; ++t){
        cum += prob_existing[t];
        if(u < cum){t_star = t; break;}
      }
      if(t_star == -1) t_star = T; // new table
      
      if(t_star < T){
        add_customer_to_existing_table(i, t_star);
      } else {
        int t_new = T;
        create_empty_table();
        add_customer_to_new_table(i, t_new);
        assign_dishes_new_table(i, t_new);
      }
    }
    
    // Step 5: update hyperparameters (alpha_v, sigma_v, tau_v) via MH
    update_hyperparameters_MH();
    
    // Step 6: store state if needed
    if (iter >= burn_in && ((iter - burn_in) % thin == 0)) {
      save_state();
    }
  }
}


double compute_log_likelihood() {
  double ll = 0.0;
  for (int i = 0; i < n; ++i) {
    for (int v = 0; v < d; ++v) {
      int t = table_of[i];
      int k = dish_of[v][t];
      double y_vi = y[v][i];
      // append likelihood contribution
    }
  }
  return ll;
}