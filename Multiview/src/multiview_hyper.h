// multiview_hyper.h
#ifndef MULTIVIEW_HYPER_H
#define MULTIVIEW_HYPER_H

#include "multiview_state.h"
#include <cmath>
#include <vector>

extern double alpha_global;
extern double sigma_global;
extern double a_tau;
extern double b_tau;

void initialize_hyperparameters();

double propose_tau(double tau_old);
double log_posterior_given_tau(int v, double tau_candidate);
void update_tau_v_MH();

void update_hyperparameters();

double log_EPPF(int v, double alpha, double sigma);
double log_prior_alpha(double alpha);
double log_prior_sigma(double sigma);

#endif