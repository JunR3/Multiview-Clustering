# Multiview Clustering Simulation

This repository contains a C++ implementation of a Multiview Gibbs Sampler for clustering, along with R scripts to interface with and run the simulations. 

## Project Structure

The project has recently been refactored into the following structure:

- `src/`: Contains all the raw C++ `.cpp` and `.h` source files.
- `scripts/`: Contains the R scripts (e.g., `CPP_Simulation.R`, `real.R`) used to prepare data, run the C++ sampler, and analyze the results.
- `dataset/`: Contains the datasets used by the simulations.

## Prerequisites

To run the simulations, you will need:
- R installed on your system.
- The `Rcpp` package installed in R: `install.packages("Rcpp")`
- Various other R libraries for clustering and visualization, such as `dplyr`, `mcclust`, `ggplot2`, etc. (Check individual simulation scripts for full library requirements).

## How to Run Your First Simulation

1. **Set your Working Directory:**
   Everything runs from the `Multiview/` folder. Ensure your R session is working out of this directory:
   ```r
   setwd("path/to/Multiview-Clustering/Multiview")
   ```

2. **Run an R Script:**
   You can run a simulation script directly using `Rscript` from the command line, or by executing the file in RStudio. 

   From the command line (while in the `Multiview/` folder):
   ```bash
   Rscript scripts/CPP_Simulation.R
   ```

   The R scripts are programmed to compile the C++ source files dynamically on the fly via `Rcpp::sourceCpp("src/multiview_gibbs.cpp")`.

## Adjusting Hyperparameters

The `run_gibbs_cpp` function inside the simulation scripts accepts initial hyperparameters that configure the behavior of the sampler from the R interface:

```r
res_gibbs <- run_gibbs_cpp(
    data_views = data_views,
    M          = 10000,
    burn_in    = 1000,
    thin       = 5,
    alpha_global_init = 1.0,
    sigma_global_init = 0.6,
    alpha_v_init = c(1.0, 1.5), # Vector for View configuration
    sigma_v_init = c(0.5, 0.4), # Vector for View configuration
    a_tau_prior = 2.0,          # Tau prior control
    b_tau_prior = 1.0
)
```

**Available Hyperparameters:**
- `alpha_global_init`: Initial value for the global concentration parameter. Defaults to `1.0`.
- `sigma_global_init`: Initial value for the global discount parameter. Defaults to `0.6`.
- `alpha_v_init`: Vector of initial values for view-specific concentration. Defaults to `1.0`.
- `sigma_v_init`: Vector of initial values for view-specific discount. Defaults to `0.5`.
- `tau_v_init`: Vector of initial values for view-specific precision. Defaults to an empirically driven estimate.
- `a_tau_prior`: Shape parameter for the Inverse-Gamma prior on tau. Defaults to `2.0`.
- `b_tau_prior`: Scale parameter for the Inverse-Gamma prior on tau. Defaults to `1.0`.
- `K_init_tables`: Initial number of global tables. Defaults to `4`.
- `K_init_dishes`: Initial number of dishes per view. Defaults to `2`.