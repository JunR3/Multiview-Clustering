// [[Rcpp::plugins(cpp11)]]

#ifndef FUNZIONI_RCPP_H
#define FUNZIONI_RCPP_H

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>

using namespace Rcpp;

// Dichiarazioni delle funzioni
// [[Rcpp::export]]
NumericVector count_customers_table(DataFrame x);

// [[Rcpp::export]]
NumericVector n_people_dish1(DataFrame x);

// [[Rcpp::export]]
NumericVector n_people_dish2(DataFrame x);

// [[Rcpp::export]]
NumericVector l_1(DataFrame x);

// [[Rcpp::export]]
NumericVector l_2(DataFrame x);

// [[Rcpp::export]]
double log_ratio1(int n_pers_same_dish, double sum_of_squares, double sum_of_observations, double prior_mean, double alpha, double beta);

// [[Rcpp::export]]
double ratio1(int n_pers_same_dish, double sum_of_squares, double sum_of_observations, double sum_of_observations_meno_i, double sum_of_squares_meno_i, double prior_mean, double alpha, double beta);

// [[Rcpp::export]]
double f(DataFrame x, int dish1, int dish2, List id_list, double prior_mean, double alpha, double beta);

// Probabilità di assegnare piatti già esistenti (Old1 - Old2)
// [[Rcpp::export]]
double prob_old1_old2(DataFrame m, double sigma1, double sigma2, int K1, int K2, List id, double prior_mean, double alpha, double beta);

// Probabilità di assegnare un nuovo primo piatto e un vecchio secondo piatto (New1 - Old2)
// [[Rcpp::export]]
double prob_new1_old2(DataFrame m, double gamma1, double sigma1, double sigma2, int K1, int K2, List id, double prior_mean, double alpha, double beta);

// Probabilità di assegnare un vecchio primo piatto e un nuovo secondo piatto (Old1 - New2)
// [[Rcpp::export]]
double prob_old1_new2(DataFrame m, double gamma2, double sigma1, double sigma2, int K1, int K2, List id, double prior_mean, double alpha, double beta);

// Probabilità di assegnare un nuovo primo e un nuovo secondo piatto (New1 - New2)
// [[Rcpp::export]]
double prob_new1_new2(DataFrame m, double gamma1, double gamma2, double sigma1, double sigma2, int K1, int K2, List id, double prior_mean, double alpha, double beta);

// [[Rcpp::export]]
List p_t_allocation(DataFrame m, int t);

// [[Rcpp::export]]
NumericVector prob_occupied_table(DataFrame data, List id, int tables, double prior_mean, double alpha, double beta);

// [[Rcpp::export]]
NumericVector prob_dish1_old(DataFrame data, int n_piatti1, List id, double prior_mean, double alpha, double beta);

// [[Rcpp::export]]
NumericVector prob_dish2_old(DataFrame data, int n_piatti2, List id, double prior_mean, double alpha, double beta);

// [[Rcpp::export]]
DataFrame relabel(int index, DataFrame data, int p1, int p2, int t);

// [[Rcpp::export]]
NumericVector prob_update_dish1(DataFrame temp, double gamma1, double sigma1, double K1_meno_i, List id, double mu0, double alpha, double beta);

// [[Rcpp::export]]
NumericVector prob_update_dish2(DataFrame temp, double gamma2, double sigma2, double K2_meno_i, List id, double mu0, double alpha, double beta);

#endif // FUNZIONI_RCPP_H