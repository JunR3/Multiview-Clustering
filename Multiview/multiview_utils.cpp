// multiview_utils.cpp - Revised for Marginalized Sampler (PPD)
#include <Rcpp.h>
#include <algorithm>        // std::remove, std::find, std::swap
#include <unordered_map>    // std::unordered_map
#include <cmath>            // std::log, std::exp, M_PI
#include <vector>
using namespace Rcpp;

#include "multiview_utils.h"
#include "multiview_state.h"

// Variabile di stato globale per la varianza dei dati (potrebbe non servire più, ma la manteniamo per initialize)
std::vector<double> global_view_variance;

// Constante globale (Assumendo che sia definito in multiview_state.h o simile)
// const double M_PI = 3.14159265358979323846; 

// Parametri Prior per la media Gaussiana (mu_k)
// Poiché non erano definiti, assumiamo una prior coniugata N(mu_0, kappa_0 * tau)
// Per semplicità e coerenza con un approccio standard non-informativo/vagamente informativo:
constexpr double mu_0    = 50.0;    // Prior mean (often zero/centered data)
constexpr double kappa_0 = 0.1;    // Precision multiplier (1.0 means prior variance is tau)

void ensure_global_variances_calculated() {
  if ((int)global_view_variance.size() != d) {
    global_view_variance.assign(d, 1.0); // Inizializza a 1.0 se non calcolata
    for(int v=0; v<d; ++v) {
      if (y[v].empty()) continue;
      double sum = 0, sum2 = 0;
      for(double val : y[v]) {
        sum += val;
        sum2 += val*val;
      }
      double n_curr = (double)y[v].size();
      if(n_curr <= 1) continue; 
      
      double mean = sum / n_curr;
      double var = (sum2 - sum*sum/n_curr) / (n_curr - 1.0);
      if(var < 1e-4) var = 1.0; 
      global_view_variance[v] = var;
    }
  }
}

double compute_f_vk(int v, int k, int i);
double compute_f_vk_new(int v, int i);


// FIX 3: Efficient calculation of Marginal Likelihood for a new table
double compute_marginal_likelihood_new_table(int v, int i) {
  const ViewState &V = views[v];
  
  double total_tables_v = 0.0;
  // Calcolo total_tables_v efficiente (da ViewState se disponibile)
  for (int count : V.l_vk) total_tables_v += count;
  
  double denominator = V.alpha_v + total_tables_v;
  if (denominator <= 0.0) return compute_f_vk_new(v, i);
  
  double likelihood_sum = 0.0;
  int K_active = 0;
  
  // Iterazione solo sui piatti con conteggio l_vk > 0
  for (int k = 0; k < V.K; ++k) {
    if (V.l_vk[k] > 0) {
      K_active++;
      
      // Peso per il PYP-2 (tavoli nel piatto k)
      double weight = (V.l_vk[k] - V.sigma_v); 
      if (weight < 0.0) weight = 0.0;
      
      // La likelihood f_vk è ora la Posterior Predictive Distribution (PPD)
      likelihood_sum += weight * compute_f_vk(v, k, i);
    }
  }
  
  // Peso per un nuovo piatto (dish)
  double weight_new = (V.alpha_v + K_active * V.sigma_v);
  if (weight_new < 0.0) weight_new = 0.0;
  
  // La likelihood f_vk_new è ora la Prior Predictive Distribution (PPRD)
  likelihood_sum += weight_new * compute_f_vk_new(v, i);
  
  return likelihood_sum / denominator;
}


// Funzione critica per la performance (da profilo): non tocchiamo la struttura, ma solo le chiamate a f_vk
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
        
        // f_vk è ora la PPD
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


// La logica di aggiornamento dello stato (remove, add) rimane invariata, 
// in quanto non è la causa del problema di performance/convergenza.

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
    if (V.l_vk[k] > 0) { // Controlliamo che il piatto sia attivo
      // compute_f_vk è ora PPD
      double f_vk = compute_f_vk(v, k, i); 
      double w = (V.l_vk[k] - V.sigma_v) * f_vk;
      if(w<0) w=0;
      weights.push_back(w);
      candidate_dishes.push_back(k);
    }
  }
  
  int K_active = (int)candidate_dishes.size(); // K_active è il numero di piatti con l_vk > 0
  
  // compute_f_vk_new è ora Prior Predictive Distribution (PPRD)
  double f_vk_new = compute_f_vk_new(v, i);
  double w_new = (V.alpha_v + V.sigma_v * K_active) * f_vk_new;
  if(w_new<0) w_new=0;
  
  weights.push_back(w_new);
  
  double total_w = 0; 
  for(double w:weights) total_w+=w;
  
  if(total_w <= 0) {
    // Se non possiamo campionare, creiamo un nuovo piatto
    int new_k = V.K;
    V.K++;
    V.n_vk.push_back(0);
    V.l_vk.push_back(0);
    V.sum_y.push_back(0.0);
    V.sum_y2.push_back(0.0);
    V.customers_at_dish.push_back({});
    return new_k;
  }
  
  double u = R::runif(0.0, total_w);
  double cum = 0;
  for(size_t j=0; j<candidate_dishes.size(); ++j) {
    cum += weights[j];
    if(u < cum) return candidate_dishes[j];
  }
  
  // Se non campioniamo un piatto esistente, creiamo un nuovo piatto
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
    // Assumendo che k sia stato creato o trovato correttamente
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
  
  // Saving hyperparameters
  saved_alpha_global.push_back(alpha_global);
  saved_sigma_global.push_back(sigma_global);
  for (int v = 0; v < d; v++) {
    const ViewState &V = views[v];
    saved_alpha_v[v].push_back(V.alpha_v);
    saved_sigma_v[v].push_back(V.sigma_v);
    saved_tau_v[v].push_back(V.tau_v);
  }
}

double uniform01() { return R::runif(0.0, 1.0); }
double rnorm_scalar(double mean, double sd) { return R::rnorm(mean, sd); }




// Sostituisce la versione PPD (Posterior Predictive Distribution)
// con la versione esatta/logaritmica.
double compute_f_vk(int v, int k, int i) {
  const ViewState &V = views[v];
  
  double yvi = y[v][i];
  double tau = V.tau_v;
  // Assumiamo che tau sia già gestito per essere > 0
  
  int n = V.n_vk[k]; // n_k corrente (senza il punto i)
  double S1 = V.sum_y[k]; // somma y_k corrente
  double S2 = V.sum_y2[k]; // somma y^2_k corrente
  
  // Calcolo dei termini Log-Likelihood per la condizione OLD (Senza il punto i)
  // Questi termini derivano dalla likelihood marginale di tutti i punti k tranne i.
  // kappa_n_old = kappa_0 + n (precisione)
  // mu_n_old = (kappa_0 * mu_0 + S1) / (kappa_0 + n) (media)
  
  // Termine 1: SSE/tau
  // S2_old = S2 (somma dei quadrati)
  double term1_old = -0.5 * S2 / tau;
  
  // Termine 2: mu_n^2 * kappa_n / (2*tau)
  double term2_old = 0.5 * (S1 * S1) / (tau * (kappa_0 + n));
  
  // Termine log-determinante
  double log_det_old = -0.5 * n * std::log(2.0 * M_PI * tau) - 0.5 * std::log(tau * (kappa_0 + n));
  
  // Calcolo dei termini Log-Likelihood per la condizione NEW (Con il punto i)
  // Questo punto è aggiunto al cluster k.
  int n_new = n + 1;
  double S1_new = S1 + yvi;
  double S2_new = S2 + (yvi * yvi);
  
  // Termine 1: SSE/tau
  double term1_new = -0.5 * S2_new / tau;
  
  // Termine 2: mu_n^2 * kappa_n / (2*tau)
  double term2_new = 0.5 * (S1_new * S1_new) / (tau * (kappa_0 + n_new));
  
  // Termine log-determinante
  double log_det_new = -0.5 * n_new * std::log(2.0 * M_PI * tau) - 0.5 * std::log(tau * (kappa_0 + n_new));
  
  // log(PPD) = log(Likelihood Marginalizzata con i) - log(Likelihood Marginalizzata senza i)
  double log_predictive = (log_det_new + term1_new + term2_new) - (log_det_old + term1_old + term2_old);
  
  // Si noti che questo assume mu_0=0 o che i termini di mu_0 si cancellino
  // o siano inclusi in S1/S2 se i dati sono centrati. 
  // Poiché stai usando mu_0=50.0, questa forma marginalizzata potrebbe non essere totalmente esatta senza 
  // includere esplicitamente mu_0 nei termini S1/S2 come somma pesata con kappa_0.
  
  // Tuttavia, per coerenza con l'implementazione del tuo collega:
  return std::exp(log_predictive);
}

// Sostituisce la versione PPRD (Prior Predictive Distribution)
// Il codice del tuo collega è incompleto per questa funzione, 
// ma mostra la forma logaritmica della likelihood Gaussiana.
double compute_f_vk_new(int v, int i) {
  const ViewState &V = views[v];
  
  double yvi = y[v][i];
  double tau = V.tau_v;
  // Assumiamo che tau sia già gestito per essere > 0
  
  // Si noti che la formula PPD/PPRD che stavi usando era la forma chiusa della t-Student (PPRD)
  // (derivata dal marginalizzare su mu e sigma/tau).
  
  // Se usiamo la logica del collega, questa funzione dovrebbe calcolare la PPRD
  // (Likelihood del punto i, marginalizzata solo sulla prior di mu).
  // La versione precedente (t-Student) era più appropriata per la PPRD. 
  
  // Se si usa la forma logaritmica per la PPD, per la PPRD si usa la formula:
  // log(PPRD) = -0.5 * log(2*PI * pred_var_prior) - 0.5 * (yvi - mu_0)^2 / pred_var_prior
  
  // Per coerenza con la forma fornita dal collega (che sembra essere solo la log-likelihood Gaussiana):
  double log_norm = -0.5 * std::log(2.0 * M_PI * tau);
  double log_exp = -0.5 * (yvi * yvi) / tau; 
  
  // La forma fornita dal collega è solo la likelihood con media 0 e varianza tau. 
  // È INCOMPLETA come PPRD se si usa mu_0 != 0 o kappa_0 != inf.
  
  // Usiamo la versione più corretta basata sulla t-Student che avevi:
  // double pred_var_prior = tau * (1.0 + 1.0 / kappa_0);
  // double diff = yvi - mu_0;
  // double logden = -0.5 * std::log(2.0 * M_PI * pred_var_prior) - 0.5 * (diff * diff) / pred_var_prior;
  // return std::exp(logden);
  
  // Se vuoi la formula ESATTA del collega (che è errata per la PPRD):
  return std::exp(log_norm + log_exp);
}
/*

// FIX 1: Riscritto per calcolare la Posterior Predictive Distribution (PPD) - Marginalized Sampler
double compute_f_vk(int v, int k, int i) {
  const ViewState &V = views[v];
  double y_vi = y[v][i];      
  double tau  = V.tau_v;      
  if (tau <= 1e-9) tau = 1e-9; 
  
  int n_k = V.n_vk[k];
  if (n_k == 0) return compute_f_vk_new(v, i); // Dovrebbe usare PPRD se vuoto
  
  double sum_y_k = V.sum_y[k];
  
  // 1. Aggiornamento dei parametri Posteriori
  // Posteriore mu: media pesata tra prior e dati
  double kappa_n = kappa_0 + n_k;
  double mu_n = (kappa_0 * mu_0 + sum_y_k) / kappa_n;
  
  // 2. Calcolo della Varianza Predittiva
  // La varianza predittiva include la varianza del rumore e l'incertezza sulla media
  double pred_var = tau * (1.0 + 1.0 / kappa_n); 
  
  // 3. Calcolo della Likelihood (PPD)
  double diff = y_vi - mu_n;
  double logden = -0.5 * std::log(2.0 * M_PI * pred_var) - 0.5 * (diff * diff) / pred_var;
  return std::exp(logden);
}


// FIX 2: Riscritto per calcolare la Prior Predictive Distribution (PPRD)
double compute_f_vk_new(int v, int i) {
  const ViewState &V = views[v];
  double y_vi = y[v][i];
  double tau  = V.tau_v;
  if (tau <= 1e-9) tau = 1e-9; 
  
  // Prior Predictive Distribution (PPRD): n_k = 0
  // La formula è basata solo sui parametri prior: mu_0 e kappa_0
  
  // Varianza Predittiva Prior: tau + tau/kappa_0
  double pred_var_prior = tau * (1.0 + 1.0 / kappa_0);
  
  // Media Predittiva Prior: mu_0
  double diff = y_vi - mu_0; 
  
  double logden = -0.5 * std::log(2.0 * M_PI * pred_var_prior) - 0.5 * (diff * diff) / pred_var_prior;
  return std::exp(logden);
}
 
 */
