#include "multiview_hyper.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace {

constexpr double kEps = 1e-6;
  
  // ---------------------------------------------------------------------------
  // Empirical variance (used for initialising tau_v)
  // ---------------------------------------------------------------------------
  double empirical_variance(const std::vector<double> &vals) {
    if (vals.empty()) return 1.0;
    
    double mean = 0.0;
    for (double v : vals) mean += v;
    mean /= static_cast<double>(vals.size());
    
    double accum = 0.0;
    for (double v : vals) {
      double diff = v - mean;
      accum += diff * diff;
    }
    
    double denom = static_cast<double>(std::max<size_t>(1, vals.size() - 1));
    double var = accum / denom;
    return (var > kEps) ? var : 1.0;
  }
  
  // ---------------------------------------------------------------------------
  // View–specific log-posteriors: log p(α_v | partition_v), log p(σ_v | partition_v)
  // ---------------------------------------------------------------------------
  double log_posterior_alpha_view(int v, double alpha_candidate) {
    if (alpha_candidate <= 0.0) return -INFINITY;
    const double sigma = views[v].sigma_v;
    
    double loglik  = log_EPPF(v, alpha_candidate, sigma);
    double logprior = log_prior_alpha(alpha_candidate);
    return loglik + logprior;
  }
  
  double log_posterior_sigma_view(int v, double sigma_candidate) {
    if (sigma_candidate <= 0.0 || sigma_candidate >= 1.0) return -INFINITY;
    const double alpha = views[v].alpha_v;
    
    double loglik  = log_EPPF(v, alpha, sigma_candidate);
    double logprior = log_prior_sigma(sigma_candidate);
    return loglik + logprior;
  }
  
  // ---------------------------------------------------------------------------
  // Global EPPF for the franchise (global α, σ)
  // ---------------------------------------------------------------------------
  double log_global_EPPF(double alpha, double sigma) {
    if (!(sigma > 0.0 && sigma < 1.0)) return -INFINITY;
    if (alpha <= -sigma) return -INFINITY;
    if (T <= 0 || n_t.empty()) return 0.0;
    
    double logp = 0.0;
    
    // Product over table creation terms
    for (int j = 0; j < T; ++j) {
      double term = alpha + j * sigma;
      if (term <= 0.0) return -INFINITY;
      logp += std::log(term);
    }
    
    // Total customers inside tables
    int total_customers = 0;
    for (int count : n_t) total_customers += count;
    
    // Normalising term
    for (int i = 1; i < total_customers; ++i) {
      double term = alpha + i;
      if (term <= 0.0) return -INFINITY;
      logp -= std::log(term);
    }
    
    // Discount contributions
    for (int count : n_t) {
      for (int m = 1; m < count; ++m) {
        double term = static_cast<double>(m) - sigma;
        if (term <= 0.0) return -INFINITY;
        logp += std::log(term);
      }
    }
    
    return logp;
  }
  
  // ---------------------------------------------------------------------------
  // Global α, σ log-posteriors
  // ---------------------------------------------------------------------------
  double log_posterior_global_alpha(double alpha_candidate) {
    double loglik  = log_global_EPPF(alpha_candidate, sigma_global);
    double logprior = log_prior_alpha(alpha_candidate);
    return loglik + logprior;
  }
  
  double log_posterior_global_sigma(double sigma_candidate) {
    double loglik  = log_global_EPPF(alpha_global, sigma_candidate);
    double logprior = log_prior_sigma(sigma_candidate);
    return loglik + logprior;
  }
  
  // ---------------------------------------------------------------------------
  // Proposal for α_v and α_global  (log-normal random walk)
  // ---------------------------------------------------------------------------
  double propose_alpha(double alpha_old) {
    const double step = 0.1;
    
    double log_alpha = std::log(std::max(alpha_old, kEps));
    log_alpha += rnorm(0.0, step);
    
    double candidate = std::exp(log_alpha);
    return (candidate > kEps) ? candidate : kEps;
  }
  
  // ---------------------------------------------------------------------------
  // Reflection mechanism for σ proposals
  // ---------------------------------------------------------------------------
  double reflect_into_unit_interval(double value) {
    double prop = value;
    
    while (prop <= kEps || prop >= 1.0 - kEps) {
      if (prop <= kEps)
        prop = 2.0 * kEps - prop;
      
      if (prop >= 1.0 - kEps)
        prop = 2.0 * (1.0 - kEps) - prop;
    }
    
    return std::clamp(prop, kEps, 1.0 - kEps);
  }
  
  // ---------------------------------------------------------------------------
  // Proposal for σ_v and σ_global (random walk + reflection)
  // ---------------------------------------------------------------------------
  double propose_sigma(double sigma_old) {
    const double step = 0.05;
    double proposal = sigma_old + rnorm(0.0, step);
    return reflect_into_unit_interval(proposal);
  }
  
} // end namespace



// ============================================================================
// GLOBAL CONSTANTS
// ============================================================================
static const double a_tau = 2.0;
static const double b_tau = 1.0;



// ============================================================================
// INITIALISATION OF HYPERPARAMETERS
// ============================================================================
void initialize_hyperparameters() {
  // --- initialise global α, σ
  if (!(alpha_global > 0.0))
    alpha_global = 1.0;
  
  if (!(sigma_global > 0.0 && sigma_global < 1.0))
    sigma_global = 0.5;
  
  // --- ensure views exist
  if (static_cast<int>(views.size()) != d)
    views.assign(d, ViewState());
  
  // --- initialise each view
  for (int v = 0; v < d; ++v) {
    ViewState &state = views[v];
    
    if (state.alpha_v <= 0.0)
      state.alpha_v = 1.0;
    
    if (!(state.sigma_v > 0.0 && state.sigma_v < 1.0))
      state.sigma_v = 0.5;
    
    double tau_seed = 1.0;
    if (v < static_cast<int>(y.size()))
      tau_seed = empirical_variance(y[v]);
    
    state.tau_v = tau_seed;
  }
}



// ============================================================================
// TAU_v UPDATES (per-view Metropolis–Hastings)
// ============================================================================
double propose_tau(double tau_old) {
  const double step_size = 0.1;
  
  double log_tau_old = std::log(tau_old);
  double eps         = rnorm(0.0, step_size);
  double log_tau_prop = log_tau_old + eps;
  
  return std::exp(log_tau_prop);
}

double log_posterior_given_tau(int v, double tau_candidate) {
  const ViewState &V = views[v];
  
  if (tau_candidate <= 0.0)
    return -INFINITY;
  
  // ---- Likelihood
  double loglik = 0.0;
  for (int k = 0; k < V.K; ++k) {
    int n_k       = V.n_vk[k];
    double sum_y2 = V.sum_y2[k];
    
    if (n_k == 0) continue;
    
    double term =
      -0.5 * n_k * std::log(2.0 * M_PI * tau_candidate)
      -0.5 * (sum_y2 / tau_candidate);
      
      loglik += term;
  }
  
  // ---- Prior tau ~ Inv-Gamma(a_tau, b_tau)
  double logprior =
    a_tau * std::log(b_tau)
    - std::lgamma(a_tau)
    - (a_tau + 1.0) * std::log(tau_candidate)
    - b_tau / tau_candidate;
    
    return loglik + logprior;
}

void update_tau_v_MH() {
  for (int v = 0; v < d; ++v) {
    ViewState &V = views[v];
    
    double tau_old      = V.tau_v;
    double log_old      = log_posterior_given_tau(v, tau_old);
    
    double tau_prop     = propose_tau(tau_old);
    if (tau_prop <= 0.0) continue;
    
    double log_new      = log_posterior_given_tau(v, tau_prop);
    
    double log_acc = log_new - log_old;
    if (std::log(uniform01()) < log_acc)
      V.tau_v = tau_prop;
  }
}



// ============================================================================
// UPDATE α_v, σ_v, α_global, σ_global
// ============================================================================
void update_hyperparameters() {
  if (views.empty())
    initialize_hyperparameters();
  
  // update all tau_v
  update_tau_v_MH();
  
  // update each view’s α_v and σ_v
  for (int v = 0; v < d; ++v) {
    ViewState &V = views[v];
    
    // ---- α_v
    double alpha_old = V.alpha_v;
    double alpha_prop = propose_alpha(alpha_old);
    
    if (std::log(uniform01()) <
      log_posterior_alpha_view(v, alpha_prop)
          - log_posterior_alpha_view(v, alpha_old))
    {
      V.alpha_v = alpha_prop;
    }
    
    // ---- σ_v
    double sigma_old = V.sigma_v;
    double sigma_prop = propose_sigma(sigma_old);
    
    if (std::log(uniform01()) <
      log_posterior_sigma_view(v, sigma_prop)
          - log_posterior_sigma_view(v, sigma_old))
    {
      V.sigma_v = sigma_prop;
    }
  }
  
  // ---- α_global
  double ag_old = alpha_global;
  double ag_prop = propose_alpha(ag_old);
  
  if (std::log(uniform01()) <
    log_posterior_global_alpha(ag_prop)
        - log_posterior_global_alpha(ag_old))
  {
    alpha_global = ag_prop;
  }
  
  // ---- σ_global
  double sg_old = sigma_global;
  double sg_prop = propose_sigma(sg_old);
  
  if (std::log(uniform01()) <
    log_posterior_global_sigma(sg_prop)
        - log_posterior_global_sigma(sg_old))
  {
    sigma_global = sg_prop;
  }
}



// ============================================================================
// EPPF for a single view (Pitman–Yor partition)
// ============================================================================
double log_EPPF(int v, double alpha, double sigma) {
  if (v < 0 || v >= d) return -INFINITY;
  if (!(sigma > 0.0 && sigma < 1.0)) return -INFINITY;
  if (alpha <= -sigma) return -INFINITY;
  if (views.empty()) return -INFINITY;
  
  const ViewState &V = views[v];
  if (V.n_vk.empty()) return 0.0;
  
  // extract sizes of non-empty clusters
  std::vector<int> cluster_sizes;
  cluster_sizes.reserve(V.n_vk.size());
  
  int total_n = 0;
  for (int count : V.n_vk) {
    if (count > 0) {
      cluster_sizes.push_back(count);
      total_n += count;
    }
  }
  
  if (total_n == 0) return 0.0;
  
  // ---- Compute EPPF
  double logp = 0.0;
  
  // table creation terms
  const int K_active = static_cast<int>(cluster_sizes.size());
  for (int j = 0; j < K_active; ++j) {
    double term = alpha + j * sigma;
    if (term <= 0.0) return -INFINITY;
    logp += std::log(term);
  }
  
  // normalisation
  for (int i = 1; i < total_n; ++i) {
    double term = alpha + i;
    if (term <= 0.0) return -INFINITY;
    logp -= std::log(term);
  }
  
  // discount terms
  for (int n_k : cluster_sizes) {
    for (int m = 1; m < n_k; ++m) {
      double term = m - sigma;
      if (term <= 0.0) return -INFINITY;
      logp += std::log(term);
    }
  }
  
  return logp;
}



// ============================================================================
// PRIORS FOR α AND σ
// ============================================================================
double log_prior_alpha(double alpha) {
  if (alpha <= 0.0) return -INFINITY;
  
  const double shape = 2.0;
  const double rate  = 1.0;
  
  return (shape - 1.0) * std::log(alpha) - rate * alpha;
}

double log_prior_sigma(double sigma) {
  if (sigma <= 0.0 || sigma >= 1.0) return -INFINITY;
  
  const double a = 2.0;
  const double b = 2.0;
  
  return (a - 1.0) * std::log(sigma) + (b - 1.0) * std::log(1.0 - sigma);
}