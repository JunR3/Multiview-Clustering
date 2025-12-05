#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include "multiview_functions_bayesian_project.h"
using namespace Rcpp;

/// Update

//// GAMMA ////
// [[Rcpp::export]]
double f_gamma(double gamma, double sigma, int tavoli, int n, double a, double b) {
  double vec = 0;
  for(int i = 1; i <= (tavoli - 1); ++i){
    vec += log(gamma + (i*sigma));
  }
  double l_sum_gamma = vec - ::lgamma(gamma + n) + ::lgamma(gamma + 1);
  return((l_sum_gamma + (a - 1)*log(gamma) + (-b*gamma)));
}
// [[Rcpp::export]]
double log_ratio(double gamma_star, double gamma_old, double sigma, int tavoli, int n, double a, double b ){
  return (f_gamma(gamma_star, sigma, tavoli, n, a, b) - f_gamma(gamma_old, sigma, tavoli, n, a, b) + R::dgamma(gamma_old, gamma_star, 1, 1) - R::dgamma(gamma_star, gamma_old, 1, 1));
}
// [[Rcpp::export]]
double ratio_gamma(double gamma_star, double gamma_old, double sigma, int tavoli, int n, double a, double b){
  return exp(log_ratio(gamma_star, gamma_old, sigma, tavoli, n, a, b));
}

/// SIGMA ////
// [[Rcpp::export]]
double f_sigma(DataFrame x, double gamma, double sigma, int n_tavoli, int K1, int K2, double alpha, double beta) {
  double temp = 0.0;
  List piatto = p_t_allocation(x, n_tavoli);
  NumericVector piatto1 = piatto[0];
  NumericVector piatto2 = piatto[1];
  NumericVector q_t = count_customers_table(x);
  
  for (int i = 1; i <= K1; ++i) {
    for (int j = 1; j <= K2; ++j) {
      LogicalVector index = (piatto1 == i) & (piatto2 == j);
      if (is_true(any(index))) {
        NumericVector q_t_index = q_t[index];
        for (int k = 0; k < q_t_index.size(); ++k) {
          temp += ::lgamma(q_t_index[k] - sigma);
        }
      }
    }
  }
  double l_sum_sigma = temp - n_tavoli *(::lgamma(1.0 - sigma));
  for (int i = 1; i < n_tavoli; ++i) {
    l_sum_sigma += log(gamma + i * sigma);
  }
  return (l_sum_sigma + (alpha - 1)*log(sigma) + (beta - 1)*log(1 - sigma));
}

// [[Rcpp::export]]
double log_ratio_sigma(DataFrame x, double gamma, double sigma_old, double sigma_star, int n_tavoli, int K1, int K2, double alpha, double beta){
  return(f_sigma(x, gamma, sigma_star, n_tavoli, K1, K2, alpha, beta) - f_sigma(x, gamma, sigma_old, n_tavoli, K1, K2, alpha, beta));
}

// [[Rcpp::export]]
double ratio_sigma(DataFrame x, double gamma, double sigma_old, double sigma_star, int n_tavoli, int K1, int K2, double alpha, double beta){
  return exp(log_ratio_sigma(x, gamma, sigma_old, sigma_star, n_tavoli, K1, K2, alpha, beta));
}

//// GAMMA1 ////
// [[Rcpp::export]]
double f_gamma1(DataFrame x, double gamma1, double sigma1, int n, double a, double b) {
  int n_piatti1;
  NumericVector l1 = l_1(x);
  n_piatti1 = l1.size();
  double vec = 0;
  for(int i = 1; i <= (n_piatti1 - 1); ++i){
    vec += log(gamma1 + (i*sigma1));
  }
  double l_sum_gamma1 = vec - ::lgamma(gamma1 + n) + ::lgamma(gamma1 + 1);
  return((l_sum_gamma1 + (a - 1)*::log(gamma1) - b*gamma1));
}
// [[Rcpp::export]]
double log_ratio_gamma1(DataFrame x, double gamma1_star, double gamma1_old, double sigma1, int n, double a, double b ){
  return (f_gamma1(x, gamma1_star, sigma1, n, a, b) - f_gamma1(x, gamma1_old, sigma1, n, a, b) + R::dgamma(gamma1_old, gamma1_star, 1, 1) - R::dgamma(gamma1_star, gamma1_old, 1, 1));
}
// [[Rcpp::export]]
double ratio_gamma1(DataFrame x, double gamma1_star, double gamma1_old, double sigma1, int n, double a, double b){
  return exp(log_ratio_gamma1(x, gamma1_star, gamma1_old, sigma1, n, a, b));
}


//// GAMMA 2 ////
// [[Rcpp::export]]
double f_gamma2(DataFrame x, double gamma2, double sigma2, int n, double a, double b) {
  NumericVector l2 = l_2(x);
  int n_piatti2 = l2.size();
  double vec = 0;
  for(int i = 1; i <= (n_piatti2 - 1); ++i){
    vec += log(gamma2 + (i*sigma2));
  }
  double l_sum_gamma2 = vec - ::lgamma(gamma2 + n) + ::lgamma(gamma2 + 1);
  return((l_sum_gamma2 + (a - 1)*::log(gamma2) - b*gamma2));
}
// [[Rcpp::export]]
double log_ratio2(DataFrame x, double gamma2_star, double gamma2_old, double sigma2, int n, double a, double b ){
  return f_gamma2(x, gamma2_star, sigma2, n, a, b) - f_gamma2(x, gamma2_old, sigma2, n, a, b) + R::dgamma(gamma2_old, gamma2_star, 1, 1) - R::dgamma(gamma2_star, gamma2_old, 1, 1);
}
// [[Rcpp::export]]
double ratio_gamma2(DataFrame x, double gamma2_star, double gamma2_old, double sigma2, int n, double a, double b){
  return exp(log_ratio2(x, gamma2_star, gamma2_old, sigma2, n, a, b));
}

//// SIGMA 1 ////

// [[Rcpp::export]]
double f_sigma1(DataFrame x, double gamma1, double sigma1, int K1, double alpha, double beta) {
  NumericVector n_k1 = n_people_dish1(x);
  double temp = 0.0;
  for (int i = 0; i < K1; ++i) {
    temp += ::lgamma(n_k1[i] - sigma1);
  }
  double l_sum_sigma1 = temp - K1*(::lgamma(1.0 - sigma1));
  for (int i = 1; i < K1; ++i) {
    l_sum_sigma1 += log(gamma1 + i * sigma1);
  }
  return (l_sum_sigma1 + (alpha - 1)*log(sigma1) + (beta - 1)*log(1 - sigma1));
}

// [[Rcpp::export]]
double log_ratio_sigma1(DataFrame x, double gamma1, double sigma_old, double sigma_star, int K1, double alpha, double beta){
  return(f_sigma1(x, gamma1, sigma_star, K1, alpha, beta) - f_sigma1(x, gamma1, sigma_old, K1, alpha, beta));
}

// [[Rcpp::export]]
double ratio_sigma1(DataFrame x, double gamma1, double sigma_old, double sigma_star, int K1, double alpha, double beta){
  return exp(log_ratio_sigma1(x, gamma1, sigma_old, sigma_star, K1, alpha, beta));
}

///// SIGMA 2 /////

// [[Rcpp::export]]
double f_sigma2(DataFrame x, double gamma2, double sigma2, int K2, double alpha, double beta) {
  NumericVector n_k2 = n_people_dish2(x);
  double temp = 0.0;
  for (int i = 0; i < K2; ++i) {
    temp += ::lgamma(n_k2[i] - sigma2);
  }
  double l_sum_sigma2 = temp - K2*(::lgamma(1.0 - sigma2));
  for (int i = 1; i < K2; ++i) {
    l_sum_sigma2 += log(gamma2 + i * sigma2);
  }
  return (l_sum_sigma2 + (alpha - 1)*log(sigma2) + (beta - 1)*log(1 - sigma2));
}

// [[Rcpp::export]]
double log_ratio_sigma2(DataFrame x, double gamma2, double sigma_old, double sigma_star, int K2, double alpha, double beta){
  return(f_sigma2(x, gamma2, sigma_star, K2, alpha, beta) - f_sigma2(x, gamma2, sigma_old, K2, alpha, beta));
}

// [[Rcpp::export]]
double ratio_sigma2(DataFrame x, double gamma2, double sigma_old, double sigma_star, int K2, double alpha, double beta){
  return exp(log_ratio_sigma2(x, gamma2, sigma_old, sigma_star, K2, alpha, beta));
}


////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
double update_gamma(double gamma_old, double sigma, int const tavoli, int const n, int a, int b) {
  double gamma_star;
  gamma_star = R::rgamma(gamma_old, 1);
  //gamma_provati <- 1
  while (gamma_star == 0) {
    gamma_star = R::rgamma(gamma_old, 1);
    //gamma_provati <- gamma_provati + 1
  }
  
  // Calcola il rapporto delle densità
  // a,b sono valori dei parametri per la prior per gamma: gamma ~ Gamma(a, b)
  double ratio = ratio_gamma(gamma_star, gamma_old, sigma, tavoli, n, a, b);
  NumericVector choices = NumericVector::create(1, ratio);
  double out = *std::min_element(choices.begin(), choices.end());
  double u = R::runif(0,1);
  if (u <= out) {
    return(gamma_star);
  } else {
    return(gamma_old);
  }
}

// [[Rcpp::export]]
double update_gamma1(DataFrame x, double gamma1_old, double sigma1, int const n, int a, int b) {
  double gamma1_star = R::rgamma(gamma1_old, 1);
  //gamma1_provati <- 1
  while (gamma1_star == 0) {
    gamma1_star = R::rgamma(gamma1_old, 1);
    //gamma1_provati <- gamma1_provati + 1
  }
  
  //Calcola il rapporto delle densità
  //a,b sono valori dei parametri per la prior per gamma: gamma ~ Gamma(a, b)
  double ratio = ratio_gamma1(x, gamma1_star, gamma1_old, sigma1, n, a, b);
  
  //Calcola il minimo tra 1 e il rapporto
  NumericVector choices = NumericVector::create(1, ratio);
  double out_gamma1 = *std::min_element(choices.begin(), choices.end());
  double u = R::runif(0,1);
  //Decide se accettare o rifiutare il valore proposto
  if (u <= out_gamma1) {
    return(gamma1_star);
  } else {
    return(gamma1_old);
  }
}

// [[Rcpp::export]]
double update_gamma2(DataFrame x, double gamma2_old, double sigma2, const int n, const int a, const int b) {
  double gamma2_star = R::rgamma(gamma2_old, 1);
  while (gamma2_star == 0) {
    gamma2_star = R::rgamma(gamma2_old, 1);
  }
  
  //Calcola il rapporto delle densità
  //a,b sono valori dei parametri per la prior per gamma: gamma ~ Gamma(a, b)
  double ratio = ratio_gamma2(x, gamma2_star, gamma2_old, sigma2, n, a, b);
  
  //Calcola il minimo tra 1 e il rapporto
  NumericVector choices = NumericVector::create(1, ratio);
  double out_gamma2 = *std::min_element(choices.begin(), choices.end());
  double u = R::runif(0,1);
  //Decide se accettare o rifiutare il valore proposto
  if (u <= out_gamma2) {
    return(gamma2_star);
  } else {
    return(gamma2_old);
  }
}


// [[Rcpp::export]]
double update_sigma(DataFrame x, double gamma, double sigma_old, int n_tavoli, int K1, int K2, double alpha, double beta){
  double sig_star = R::runif(0,1);
  double ratio = ratio_sigma(x, gamma, sigma_old, sig_star, n_tavoli, K1, K2, alpha, beta);
  
  NumericVector choices = NumericVector::create(1, ratio);
  double out_sig = *std::min_element(choices.begin(), choices.end());
  double u_sig = R::runif(0,1);
  if (u_sig <= out_sig) {
    return(sig_star);
  } else {
    return(sigma_old);
  }
}

// [[Rcpp::export]]
double update_sigma1(DataFrame x, double gamma1, double sigma1_old, int K1, int alpha, int beta){
  double sig1_star = R::runif(0, 1);
  double ratio = ratio_sigma1(x, gamma1, sigma1_old, sig1_star, K1, alpha, beta);
  
  NumericVector choices = NumericVector::create(1, ratio);
  double out_sig1 = *std::min_element(choices.begin(), choices.end());
  double u_sig1 = R::runif(0,1);
  
  if (u_sig1 <= out_sig1) {
    return(sig1_star);
  } else {
    return(sigma1_old);
  }
}

// [[Rcpp::export]]
double update_sigma2(DataFrame x, double gamma2, double sigma2_old, int K2, int alpha, int beta){
  double sig2_star = R::runif(0, 1);
  double ratio = ratio_sigma2(x, gamma2, sigma2_old, sig2_star, K2, alpha, beta);
  
  NumericVector choices = NumericVector::create(1, ratio);
  double out_sig2 = *std::min_element(choices.begin(), choices.end());
  double u_sig2 = R::runif(0,1);
  
  if (u_sig2 <= out_sig2) {
    return(sig2_star);
  } else {
    return(sigma2_old);
  }
}