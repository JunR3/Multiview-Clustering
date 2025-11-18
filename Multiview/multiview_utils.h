//
// Created by rafam on 18/11/2025.
//

#ifndef MULTIVIEW_UTILS_H
#define MULTIVIEW_UTILS_H

#include <vector>
#include "multiview_state.h"

// --- Likelihood helpers ---
double compute_f_vk(int v, int k, int i);
double compute_f_vk_new(int v, int i);

// --- Table assignment helpers ---
void compute_table_probs_with_cache(
    int i,
    std::vector<double> &prob_existing,
    double &prob_new
);

int sample_table(int i,
                 const std::vector<double> &prob_existing,
                 double prob_new);

// --- Customer movement ---
void remove_customer(int i);
void add_customer_to_existing_table(int i, int t);

void create_empty_table(int t_new);
void add_customer_to_new_table(int i, int t_new);

// --- Dish assignment for new tables ---
int sample_dish_for_new_table(int v, int i);
void assign_dishes_new_table(int i, int t_new);

void save_state();

#endif // MULTIVIEW_UTILS_H
