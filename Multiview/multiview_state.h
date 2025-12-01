//multiview_state.h
#ifndef MULTIVIEW_STATE_H
#define MULTIVIEW_STATE_H

#include <vector>

struct ViewState {
  int K = 0;                                
  std::vector<int> n_vk;                    
  std::vector<int> l_vk;                    
  std::vector<double> sum_y;                
  std::vector<double> sum_y2;               
  std::vector<std::vector<int>> customers_at_dish; 
  
  double alpha_v = 1.0;                     // α_v (local concentration)
  double sigma_v = 0.5;                     // σ_v (local discount)
  double tau_v   = 1.0;                     // τ_v (kernel variance)
};


extern int n, d;                            // n = # clienti, d = # views
extern std::vector<std::vector<double>> y;  // y[v][i] = dato cliente i nella view v

extern double alpha_global;                 // α globale (franchise)
extern double sigma_global;                 // σ globale (franchise)

extern int T;                              
extern std::vector<int> table_of;          // table_of[i] = tavolo del cliente i
extern std::vector<int> n_t;               // n_t[t] = # clienti al tavolo t
extern std::vector<std::vector<int>> customers_at_table; // clienti per tavolo
extern std::vector<std::vector<int>> dish_of;            // dish_of[v][t] = dish alla view v, tavolo t

extern std::vector<ViewState> views;       

// Saved states (per output R)
extern std::vector<std::vector<int>> saved_table_of;
extern std::vector<std::vector<std::vector<int>>> saved_dish_of;
extern std::vector<double> saved_loglik;

#endif // MULTIVIEW_STATE_H