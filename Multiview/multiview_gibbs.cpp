#include <Rcpp.h>
#include "multiview_state.h"
#include "multiview_gibbs.h"

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List run_gibbs_cpp(const Rcpp::List& data_views,
                         int M, int burn_in, int thin) {
  
  d = data_views.size();
  n = Rcpp::as<int>(data_views[0].attr("n"));
  
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
    
    // qui andranno:
    // remove_customer(i)
    // compute_table_probs_with_cache
    // sample table
    // add_customer
    
    // update_hyperparameters_MH()
    
    if (iter >= burn_in && ((iter - burn_in) % thin == 0)) {
      save_state(saved_table_of.size());
    }
  }
}

void save_state(int s) {
  saved_table_of.push_back(table_of);
  saved_dish_of.push_back(dish_of);
  saved_loglik.push_back(compute_log_likelihood());
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