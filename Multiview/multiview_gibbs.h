//multiview_gibbs.h
#ifndef MULTIVIEW_GIBBS_H
#define MULTIVIEW_GIBBS_H

#include <Rcpp.h>
using namespace Rcpp;

Rcpp::List run_gibbs_cpp(const Rcpp::List& data_views,
                         int M, int burn_in, int thin);

void gibbs_sampler(int M, int burn_in, int thin);
static void initialize_state_from_data();
double compute_log_likelihood();

#endif