#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <string>
#include <unordered_set>

#include "multiview_functions_bayesian_project.h"

using namespace Rcpp;




/// Count how many people sit at a certain table:
// count_customers_table[0] = how many people at table labeled "1"
// count_customers_table[1] = how many people at table labeled "2" ...

// [[Rcpp::export]]
NumericVector count_customers_table(DataFrame x){
  // ARGS: x = dataframe with observations
  
  // Number of customers
  int n = x.nrow();
  
  // extract "table" variable
  NumericVector tables = x["tables"];
  int t = max(tables);
  NumericVector count(t);
  for(int i = 1; i < t + 1; i ++){
    for(int j = 1; j < n + 1; j ++){
      if(tables[j - 1] == i){
        count[i - 1] = count[i - 1] + 1;
      }
    }
  }
  return count;
}

/// How many customers eat a certain first dish
/// n_people_dish1[0] = how many are eating dish labeled "1"
/// n_people_dish1[1] = how many are eating dish lebeled "2" ...
// [[Rcpp::export]]
NumericVector n_people_dish1(DataFrame x){
  int n = x.nrow();
  
  NumericVector dish1 = x["dish1"];
  int t = max(dish1);
  NumericVector count(t);
  for(int i = 1; i < t + 1; i ++){
    for(int j = 1; j < n + 1; j ++){
      if(dish1[j - 1] == i){
        count[i - 1] = count[i - 1] + 1;
      }else{
        count[i - 1] = count[i - 1];
      }
    }
  }
  return count;
}

// How many customers eat a certain second dish
// [[Rcpp::export]]
NumericVector n_people_dish2(DataFrame x){
  int n = x.nrow();
  
  NumericVector dish2 = x["dish2"];
  int t = max(dish2);
  NumericVector count(t);
  for(int i = 1; i < t + 1; i ++){
    for(int j = 1; j < n + 1; j ++){
      if(dish2[j - 1] == i){
        count[i - 1] = count[i - 1] + 1;
      }else{
        count[i - 1] = count[i - 1];
      }
    }
  }
  return count;
}



/// For every first dish, how many tables are serving it
/// l_1[0] = how many tables are serving the first dish1
/// l_1[1] = how many tables are serving the second dish1
// [[Rcpp::export]]
NumericVector l_1(DataFrame x) {
  int n = x.nrow();
  
  const NumericVector& dish1 = x["dish1"];
  const NumericVector& tables = x["tables"];
  
  std::unordered_set<int> unique_tables;
  
  int K1 = *std::max_element(dish1.begin(), dish1.end());
  NumericVector count(K1);
  
  for (int i = 1; i <= K1; ++i) {
    unique_tables.clear(); 
    
    for (int j = 0; j < n; ++j) {
      if (dish1[j] == i) {
        unique_tables.insert(tables[j]);
      }
    }
    count[i - 1] = unique_tables.size();
  }
  return count;
}

/// For every second dish, how many tables are serving it
/// the example is the same as before
// [[Rcpp::export]]
NumericVector l_2(DataFrame x){
  int n = x.nrow();
  
  NumericVector dish2 = x["dish2"]; 
  NumericVector tables = x["tables"];
  
  int K2 = unique(dish2).length();
  
  NumericVector count(K2);
  NumericVector index(n);
  
  for(int i = 1; i < K2 + 1; i ++){
    for(int j = 1; j < n + 1; j ++){
      if(dish2[j - 1] == i){
        index[j - 1] = tables[j - 1];
      }else{
        index[j - 1] = 0;
      }
    }
    count[i - 1] = unique(index).length() - 1;
  }
  return count;
}



/// useful to compute the density of y_i fiven all the other observations in that cluster
// [[Rcpp::export]]
double log_ratio1(int n_pers_same_dish, double sum_of_squares, double sum_of_observations, double prior_mean, double alpha, double beta){
  //
  // Args:
  //   n_pers_same_dish: Number of observations.
  //   sum_of_squares: Sum of squared observations.
  //   sum_of_observations: Sum of observations.
  //   prior_mean: Prior mean.
  //   alpha, beta: inv.Gamma prior parameters
  
  double log_term1 = -(n_pers_same_dish/2.0) * log(2 * M_PI);
  double log_term2 = - 0.5 * log(n_pers_same_dish + 1.0) + log(tgamma(alpha + (n_pers_same_dish + 1.0)/2));
  double log_term3 = (alpha + (n_pers_same_dish + 1.0)/2)*log(0.5*(sum_of_squares + prior_mean*prior_mean) - 0.5*((prior_mean + sum_of_observations)*(prior_mean + sum_of_observations))/(n_pers_same_dish + 1.0) + beta);
  
  return log_term1 + log_term2 - log_term3;
}

// [[Rcpp::export]]
double ratio1(int n_pers_same_dish, double sum_of_squares, double sum_of_observations,  double sum_of_observations_meno_i, double sum_of_squares_meno_i, double prior_mean, double alpha, double beta) {
  return exp(log_ratio1(n_pers_same_dish + 1, sum_of_squares, sum_of_observations, prior_mean, alpha, beta) - log_ratio1(n_pers_same_dish, sum_of_squares_meno_i, sum_of_observations_meno_i, prior_mean, alpha, beta));
}

// computationally expensive !!!
// [[Rcpp::export]]
double f(DataFrame x, int dish1, int dish2, List id_list, double prior_mean, double alpha, double beta) {
  // Args:
  //   id_list = list containing the info about the subject we are updating
  
  int n = x.nrow();
  // extract information abut the subject
  NumericVector id(4);
  id[0] = as<double>(id_list["view1"]);
  id[1] = as<double>(id_list["view2"]);
  id[2] = as<double>(id_list["view1_sq"]);
  id[3] = as<double>(id_list["view2_sq"]);
  
  NumericVector p1 = x["dish1"];
  NumericVector p2 = x["dish2"];
  NumericVector view1 = x["view1"];
  NumericVector view2 = x["view2"];
  NumericVector view1_sq = x["view1_sq"];
  NumericVector view2_sq = x["view2_sq"];
  
  std::vector<int> dish1_indices;
  dish1_indices.reserve(n); 
  std::vector<int> dish2_indices;
  dish2_indices.reserve(n);
  
  // which customers eat the same first and second dish
  for (int i = 0; i < n; i++) {
    if (p1[i] == dish1) {
      dish1_indices.push_back(i);
    }
    if (p2[i] == dish2) {
      dish2_indices.push_back(i);
    }
  }
  // how many
  int n_dish1 = dish1_indices.size(); 
  int n_dish2 = dish2_indices.size();
  
  double sum_view1_minus_i = 0;
  double sum_view1_sq_minus_i = 0;
  
  // for dish 1 and dish 2, i need to compute these quantities to evaluate f()
  for (int index : dish1_indices) {
    sum_view1_minus_i += view1[index];
    sum_view1_sq_minus_i += view1_sq[index];
  }
  
  double sum_view1 = sum_view1_minus_i + id[0];
  double sum_view1_sq = sum_view1_sq_minus_i + id[2];
  
  double sum_view2_minus_i = 0;
  double sum_view2_sq_minus_i = 0;
  
  for (int index : dish2_indices) {
    sum_view2_minus_i += view2[index];
    sum_view2_sq_minus_i += view2_sq[index];
  }
  
  double sum_view2 = sum_view2_minus_i + id[1];
  double sum_view2_sq = sum_view2_sq_minus_i + id[3];
  

  return ratio1(n_dish1, sum_view1_sq, sum_view1, sum_view1_minus_i, sum_view1_sq_minus_i, prior_mean, alpha,  beta)*ratio1(n_dish2, sum_view2_sq, sum_view2, sum_view2_minus_i, sum_view2_sq_minus_i, prior_mean, alpha,  beta);
}


// Probability of assigning a first and a second dish that are already served at the restaurant
// [[Rcpp::export]]
double prob_old1_old2(DataFrame m, double sigma1, double sigma2, int K1, int K2, List id, double prior_mean, double alpha, double beta) {
  NumericVector l1 = l_1(m);
  NumericVector l2 = l_2(m);
  double vec = 0.0;
  
  for (int i = 0; i < K1; ++i) {
    for (int j = 0; j < K2; ++j) {
      vec += (l1[i] - sigma1) * (l2[j] - sigma2) * f(m, i + 1, j + 1, id, prior_mean, alpha, beta);
    }
  }
  
  return vec;
}

// Probability of assigning a new first dish and a second dish that is already served at the restaurant
// [[Rcpp::export]]
double prob_new1_old2(DataFrame m, double gamma1, double sigma1, double sigma2, int K1, int K2, List id, double prior_mean, double alpha, double beta) {
  NumericVector l2 = l_2(m);
  double vec = 0.0;
  
  for (int j = 0; j < K2; ++j) {
    vec += (l2[j] - sigma2) * f(m, K1 + 1, j + 1, id, prior_mean, alpha, beta);
  }
  return (gamma1 + sigma1*K1)*vec;
}

// Probability of assigning an old first dish and a new second dish 
// [[Rcpp::export]]
double prob_old1_new2(DataFrame m, double gamma2, double sigma1, double sigma2, int K1, int K2, List id, double prior_mean, double alpha, double beta) {
  NumericVector l1 = l_1(m);
  double vec = 0.0;
  
  for (int i = 0; i < K1; ++i) {
    vec += (l1[i] - sigma1) * f(m, i + 1, K2 + 1, id, prior_mean, alpha, beta); // +1 per l'indicizzazione R
  }
  
  return (gamma2 + sigma2*K2)*vec;
}

// Probability of assigning a new first and second dish 
// [[Rcpp::export]]
double prob_new1_new2(DataFrame m, double gamma1, double gamma2, double sigma1, double sigma2, int K1, int K2, List id, double prior_mean, double alpha, double beta) {
  return (gamma1 + sigma1*K1)*(gamma2 + sigma2*K2)*f(m,  (K1 + 1), (K2 + 1), id, prior_mean, alpha, beta);
}


/// List containing the pair of dishes associated with each table
// [[Rcpp::export]]
List p_t_allocation(DataFrame m, int t) {
  
  NumericVector dish1(t);
  NumericVector dish2(t);
  
  NumericVector tables = m["tables"];
  NumericVector dish1_m = m["dish1"];
  NumericVector dish2_m = m["dish2"];
  
  for (int i = 0; i < t; ++i) {
    for (int j = 0; j < m.nrows(); ++j) {
      if (tables[j] == (i + 1)) {
        dish1[i] = dish1_m[j];
        dish2[i] = dish2_m[j];
        break; 
      }
    }
  }
  
  return List::create(Named("dish1") = dish1, Named("dish2") = dish2);
}

// probability of sitting at an old table
// [[Rcpp::export]]
NumericVector prob_occupied_table(DataFrame data, List id, int tables, double prior_mean, double alpha, double beta){
  NumericVector p_vec(tables);
  List piatti = p_t_allocation(data, tables);
  NumericVector dish1 = piatti[0];
  NumericVector dish2 = piatti[1];
  for(int j = 0; j < tables; ++j){
    p_vec[j] = f(data, dish1[j], dish2[j], id, prior_mean, alpha, beta);
  }
  
  return p_vec;
}

// [[Rcpp::export]]
NumericVector prob_dish1_old(DataFrame data, int n_piatti1, List id, double prior_mean, double alpha, double beta){
  NumericVector p_vec(n_piatti1);
  for(int j = 0; j < n_piatti1; ++j){
    p_vec[j] = f(data, j + 1, 0, id, prior_mean, alpha, beta);
  }
  
  return p_vec;
}

// [[Rcpp::export]]
NumericVector prob_dish2_old(DataFrame data, int n_piatti2, List id, double prior_mean, double alpha, double beta){
  NumericVector p_vec(n_piatti2);
  for(int j = 0; j < n_piatti2; ++j){
    p_vec[j] = f(data, 0, j + 1, id, prior_mean, alpha, beta);
  }
  
  return p_vec;
}

// relabel step
// [[Rcpp::export]]
DataFrame relabel(int index, DataFrame data, int p1, int p2, int t){
  NumericVector n1 = n_people_dish1(data);
  NumericVector n2 = n_people_dish2(data);
  
  NumericVector q_t = count_customers_table(data);
  
  int a = q_t[t - 1];
  int b = n1[p1 - 1];
  int c = n2[p2 - 1];
  NumericVector tables = data["tables"];
  NumericVector dish1 = data["dish1"];
  NumericVector dish2 = data["dish2"];
  
  if (a == 1) { // i.e. when I remove observation "index", no more people are sitting at that table
    for (int i = 0; i < tables.size(); ++i) {
      // decrease by one all tables' labels numerated after that one
      if (tables[i] > tables[index - 1]) {
        tables[i] = tables[i] - 1;
      }
    }
  }
  
  if (b == 1) { // i.e. no cumstomer is eating that first dish once I remove observation "index"
    for (int i = 0; i < dish1.size(); ++i) {
      // decrease by one all dishes' labels numerated after that one
      if (dish1[i] > dish1[index - 1]) {
        dish1[i] = dish1[i] - 1;
      }
    }
  }
  
  if (c == 1) { // i.e. no cumstomer is eating that second dish once I remove observation "index"
    for (int i = 0; i < dish2.size(); ++i) {
      // decrease by one all dishes' labels numerated after that one
      if (dish2[i] > dish2[index - 1]) {
        dish2[i] = dish2[i] - 1;
      }
    }
  }
  
  data["tables"] = tables;
  data["dish1"] = dish1;
  data["dish2"] = dish2;
  
  return data;
}

// vector of length (number of dishes + 1), which in position 0 computes teh probability of a new first dish
// and in position 1, 2... computes the probability of assigning an old dish, corresponding to the index position 
// i.e. prob_update_dish1[0] = probability of "new dish" as first dish
//      prob_update_dish1[1] = probability of assigning "1" as first dish
//      prob_update_dish1[2] = probability of assigning "2" as first dish ...
// [[Rcpp::export]]
NumericVector prob_update_dish1(DataFrame temp, double gamma1, double sigma1, double K1_meno_i, List id, double mu0, double alpha, double beta){
  NumericVector prob_k1(K1_meno_i + 1);
  NumericVector l1 = l_1(temp);
  prob_k1[0] = (gamma1 + sigma1*K1_meno_i)*f(temp, K1_meno_i + 1, 0, id, mu0, alpha, beta);
  for(int i = 1; i <= K1_meno_i; ++i){
    prob_k1[i] = (l1[i - 1] - sigma1)*prob_dish1_old(temp, K1_meno_i, id, mu0, alpha, beta)[i - 1];
  }
  return prob_k1;
}

// this function is the same as the previous one
// [[Rcpp::export]]
NumericVector prob_update_dish2(DataFrame temp, double gamma2, double sigma2, double K2_meno_i, List id, double mu0, double alpha, double beta){
  NumericVector prob_k2(K2_meno_i + 1);
  NumericVector l2 = l_2(temp);
  prob_k2[0] = (gamma2 + sigma2*K2_meno_i)*f(temp, 0, K2_meno_i + 1, id, mu0, alpha, beta);
  for(int i = 1; i <= K2_meno_i; ++i){
    prob_k2[i] = (l2[i - 1] - sigma2)*prob_dish2_old(temp, K2_meno_i, id, mu0, alpha, beta)[i - 1];
  }
  return prob_k2;
}