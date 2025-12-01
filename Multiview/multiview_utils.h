// multiview_utils.h
#ifndef MULTIVIEW_UTILS_H
#define MULTIVIEW_UTILS_H

#include <vector>
#include <unordered_map> // Importante: aggiungi questo include!
#include "multiview_state.h"

// --- Global Statistics Helpers ---
void ensure_global_variances_calculated();

// --- Likelihood helpers ---
double compute_f_vk(int v, int k, int i);
double compute_f_vk_new(int v, int i);

// --- Table assignment helpers ---

// MODIFICA QUI: Aggiunto l'ultimo argomento (cache_fvk) per il workspace
void compute_table_probs_with_cache(
    int i,
    std::vector<double> &prob_existing,
    double &prob_new,
    std::vector<std::unordered_map<int, double>> &cache_fvk 
);

// --- Customer movement ---
void remove_customer(int i);
void add_customer_to_existing_table(int i, int t);

int create_empty_table();
void add_customer_to_new_table(int i, int t_new);

// --- Dish assignment for new tables ---
int sample_dish_for_new_table(int v, int i);
void assign_dishes_new_table(int i, int t_new);

// --- Utils & Random ---
void save_state();
double uniform01();
double rnorm_scalar(double mean, double sd);

#endif // MULTIVIEW_UTILS_H