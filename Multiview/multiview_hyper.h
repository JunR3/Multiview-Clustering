// multiview_hyper.h
#ifndef MULTIVIEW_HYPER_H
#define MULTIVIEW_HYPER_H

#include "multiview_state.h"  // per ViewState, views, d

// Da lavorare ancora:
// void update_hyperparameters_MH();

// Proposta random–walk su log(tau)
double propose_tau(double tau_old);

// Log-posterior di tau_v dato tutto il resto (view fissata a v)
double log_posterior_given_tau(int v, double tau_candidate);

// Pitman-Yor hyperparameters
// void update_alpha_sigma_MH(); understand how to separate local and global
double log_EPPF(int v, double alpha, double sigma);

double log_prior_alpha(double alpha);
double log_prior_sigma(double sigma);

#endif // MULTIVIEW_HYPER_H