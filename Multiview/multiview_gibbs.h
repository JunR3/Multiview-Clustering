#ifndef MULTIVIEW_GIBBS_H
#define MULTIVIEW_GIBBS_H

#include <Rcpp.h>

Rcpp::List run_gibbs_cpp(const Rcpp::List& data_views,
                         int M, int burn_in, int thin);

void gibbs_sampler(int M, int burn_in, int thin);

void save_state(int s);
double compute_log_likelihood();

#endif