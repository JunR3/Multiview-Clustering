// multiview_hyper.cpp
#include "multiview_hyper.h"
#include <cmath>

// Hyper–prior for tau_v ~ Inv-Gamma(a_tau, b_tau)
static const double a_tau = 2.0;
static const double b_tau = 1.0;

// Random-walk proposal for log(tau)
double propose_tau(double tau_old) {
  const double step_size = 0.1;         // to change if necessary
  double log_tau_old  = std::log(tau_old);
  double eps          = rnorm(0.0, step_size);   // N(0, step_size^2)
  double log_tau_prop = log_tau_old + eps;
  return std::exp(log_tau_prop);        // guarantees taht tau_prop > 0
}

// Log-posterior of tau_v 
double log_posterior_given_tau(int v, double tau_candidate) {
  const ViewState &V = views[v];
  
  if (tau_candidate <= 0.0) {
    return -INFINITY;
  }
  
  // 1. Log-likelihood
  double loglik = 0.0;

  for (int k = 0; k < V.K; ++k) {
    int    n_k    = V.n_vk[k];
    double sum_y2 = V.sum_y2[k];
    
    if (n_k == 0) continue;
    
    double term = -0.5 * n_k * std::log(2.0 * M_PI * tau_candidate)
      - 0.5 * (sum_y2 / tau_candidate);
    
    loglik += term;
  }
  
  // 2. Log-prior: tau ~ Inv-Gamma(a_tau, b_tau)

  // densità: p(tau) ∝ tau^{-(a_tau+1)} exp(-b_tau / tau)
  // log p(tau) = a_tau*log(b_tau) - lgamma(a_tau)
  //             - (a_tau + 1)*log(tau) - b_tau/tau
  
  double logprior =
    a_tau * std::log(b_tau)
    - std::lgamma(a_tau)
    - (a_tau + 1.0) * std::log(tau_candidate)
    - b_tau / tau_candidate;
    
    return loglik + logprior;
}

// MH update for all tau_v (one per view)
void update_tau_v_MH() {
  for (int v = 0; v < d; ++v) {
    ViewState &V = views[v];
    
    double tau_old       = V.tau_v;
    double log_post_old  = log_posterior_given_tau(v, tau_old);
    
    double tau_prop = propose_tau(tau_old);
    if (tau_prop <= 0.0) {
      continue;  
    }
    
    double log_post_new  = log_posterior_given_tau(v, tau_prop);
    double log_acc_ratio = log_post_new - log_post_old;
    
    double u = uniform01();   // U(0,1)
    if (std::log(u) < log_acc_ratio) {
      V.tau_v = tau_prop;    // accept the proposal
    }
    // Alternatively, leave tau_v = tau_old
  }
}