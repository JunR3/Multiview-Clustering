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

#endif // MULTIVIEW_HYPER_H