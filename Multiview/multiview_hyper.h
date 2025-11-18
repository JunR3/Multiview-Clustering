// multiview_hyper.h
#ifndef MULTIVIEW_HYPER_H
#define MULTIVIEW_HYPER_H

#include "multiview_state.h"  // per ViewState, views, d

// Aggiorna tutti i tau_v (v = 0,...,d-1) con un passo MH
void update_hyperparameters_MH();

// Proposta random–walk su log(tau)
double propose_tau(double tau_old);

// Log-posterior di tau_v dato tutto il resto (view fissata a v)
double log_posterior_given_tau(int v, double tau_candidate);

// Pitman-Yor hyperparameters
void update_alpha_sigma_MH();
double log_EPPF(int v, double alpha, double sigma);

double log_prior_alpha(double alpha);
double log_prior_sigma(double sigma);


#endif // MULTIVIEW_HYPER_H