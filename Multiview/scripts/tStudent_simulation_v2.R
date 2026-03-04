# ==============================================================================
# 0. REQUIRED LIBRARIES
# ==============================================================================
# Note: Ensure mcclust.ext is installed.
# If not: install.packages("remotes"); remotes::install_github("sarawade/mcclust.ext")
library(dplyr)
library(ggplot2)
library(Rcpp)
library(mcclust) # For comp.psm and arandi
library(mcclust.ext) # For minVI
library(gridExtra)
Rcpp::sourceCpp("src/multiview_gibbs.cpp")
# ==============================================================================
# 1. HELPER FUNCTIONS: ROBUST minVI EXTRACTION
# ==============================================================================

# Internal helper to ensure labels are consecutive integers starting from 1
# This is a strict requirement for the mcclust/mcclust.ext functions
relabel_consecutive <- function(vec) {
  return(as.integer(as.factor(vec)))
}

# Extracts all posterior samples from the Gibbs output for a specific view
extract_all_samples <- function(res_gibbs, view_idx) {
  n_iters <- length(res_gibbs$table_of)
  n_obs <- length(res_gibbs$table_of[[1]])

  # Matrix structure: Rows = Observations, Columns = MCMC Iterations
  samples_mat <- matrix(NA, nrow = n_obs, ncol = n_iters)

  for (i in 1:n_iters) {
    # Convert C++ 0-based indexing to R 1-based indexing
    tables_r <- res_gibbs$table_of[[i]] + 1
    dishes_view <- res_gibbs$dish_of[[i]][[view_idx]]

    # Extract raw labels and normalize them
    raw_labels <- dishes_view[tables_r]
    samples_mat[, i] <- relabel_consecutive(raw_labels)
  }
  return(samples_mat)
}

# Finds the optimal partition by minimizing the Variation of Information (VI) loss
get_minVI_clustering <- function(res_gibbs, view_idx) {
  cat(sprintf("--- Calculating minVI for View %d ---\n", view_idx))

  # 1. Retrieve all posterior samples
  post_samples <- extract_all_samples(res_gibbs, view_idx)

  # 2. Compute the Posterior Similarity Matrix (PSM)
  # Rows must be iterations for comp.psm, so we transpose: t()
  psm <- mcclust::comp.psm(t(post_samples))

  # 3. Find the partition that minimizes VI using mcclust.ext
  # This function summarizes the posterior uncertainty into a single "best" clustering
  vi_result <- mcclust.ext::minVI(psm, post_samples)

  return(vi_result$cl)
}

# ==============================================================================
# 2. DATA GENERATION
# ==============================================================================

set.seed(2024)
n_per_cluster <- 100
df_student <- 8
dist_val <- 8

# --- Define Dependency Structure ---
rho <- 0.5

# Construct a 3x3 Correlation Matrix for the 3 Views
# 1s on the diagonal, rho on the off-diagonals
sigma_dep <- matrix(
  c(
    1, rho, rho,
    rho, 1, rho,
    rho, rho, 1
  ),
  nrow = 3, ncol = 3
)

# --- Generate Clusters with Shared Covariance Structure ---
# Now, View 1, 2, and 3 are generated simultaneously for each point
# effectively creating a linear dependency between them.

# Cluster 1: Shifted in View 3
c1 <- mvtnorm::rmvt(
  n = n_per_cluster, sigma = sigma_dep, df = df_student,
  delta = c(0, 0, dist_val)
)

# Cluster 2: Shifted in View 2
c2 <- mvtnorm::rmvt(
  n = n_per_cluster, sigma = sigma_dep, df = df_student,
  delta = c(0, dist_val, 0)
)

# Cluster 3: Shifted in View 1
c3 <- mvtnorm::rmvt(
  n = n_per_cluster, sigma = sigma_dep, df = df_student,
  delta = c(dist_val, 0, 0)
)

# Combine data
data_matrix <- rbind(c1, c2, c3)
x_data <- data.frame(
  view1 = data_matrix[, 1],
  view2 = data_matrix[, 2],
  view3 = data_matrix[, 3]
)

true_labels <- c(rep(1, n_per_cluster), rep(2, n_per_cluster), rep(3, n_per_cluster))

# ==============================================================================
# VISUALIZING THE DEPENDENCY
# ==============================================================================
# If views are dependent, a scatter plot of View 1 vs View 2 should
# show an elongated/diagonal shape instead of a circular cloud.
plot(x_data$view1, x_data$view2,
  col = true_labels, pch = 19,
  main = paste("Cross-View Dependency (rho =", rho, ")"),
  xlab = "View 1", ylab = "View 2"
)
# ==============================================================================
# 3. SCENARIO EXECUTION
# ==============================================================================
nsim <- 2000
burn_in <- 500
thin <- 1

# --- SCENARIO 1: SINGLE VIEW ---
cat("\n>>> RUNNING SCENARIO 1: SINGLE VIEW...\n")
res_single <- run_gibbs_cpp(list(as.vector(x_data$view1)), M = nsim, burn_in = burn_in, thin = thin)
clusters_single_vi <- get_minVI_clustering(res_single, view_idx = 1)

# --- SCENARIO 2: MULTI VIEW ---
cat("\n>>> RUNNING SCENARIO 2: MULTI VIEW...\n")
res_multi <- run_gibbs_cpp(list(as.vector(x_data$view1), as.vector(x_data$view2), as.vector(x_data$view3)),
  M = nsim, burn_in = burn_in, thin = thin
)
clusters_multi_vi <- get_minVI_clustering(res_multi, view_idx = 1)

# ==============================================================================
# 4. RESULTS & VISUALIZATION
# ==============================================================================

# Calculate Adjusted Rand Index (ARI) to measure performance against ground truth
ari_single <- mcclust::arandi(clusters_single_vi, true_labels)
ari_multi <- mcclust::arandi(clusters_multi_vi, true_labels)

cat("\n==============================================\n")
cat(sprintf("ARI (Single View): %.4f\n", ari_single))
cat(sprintf("ARI (Multi View):  %.4f\n", ari_multi))
cat("==============================================\n")

# Prepare data for plotting
plot_df <- data.frame(
  x = x_data$view1,
  y = x_data$view2,
  Truth = as.factor(true_labels),
  Single_VI = as.factor(clusters_single_vi),
  Multi_VI = as.factor(clusters_multi_vi)
)

p1 <- ggplot(plot_df, aes(x, y, color = Truth)) +
  geom_point(alpha = 0.5) +
  labs(title = "Ground Truth") +
  theme_minimal()

p2 <- ggplot(plot_df, aes(x, y, color = Single_VI)) +
  geom_point(alpha = 0.5) +
  labs(title = "Single View Prediction (minVI)") +
  theme_minimal()

p3 <- ggplot(plot_df, aes(x, y, color = Multi_VI)) +
  geom_point(alpha = 0.5) +
  labs(title = "Multi View Prediction (minVI)") +
  theme_minimal()

grid.arrange(p1, p2, p3, nrow = 1)

# ==============================================================================
# 5.TRACEPLOTS
# ==============================================================================
res_gibbs <- res_multi

df_global <- data.frame(
  iter = seq_along(res_gibbs$alpha_global),
  alpha_global = res_gibbs$alpha_global,
  sigma_global = res_gibbs$sigma_global
)

p_alpha_g <- ggplot(df_global, aes(iter, alpha_global)) +
  geom_line() +
  theme_minimal() +
  labs(title = "alpha_global trace")

p_sigma_g <- ggplot(df_global, aes(iter, sigma_global)) +
  geom_line() +
  theme_minimal() +
  labs(title = "sigma_global trace")

if (is.list(res_gibbs$alpha_v)) {
  n_views <- length(res_gibbs$alpha_v)
  n_iter <- length(res_gibbs$alpha_v[[1]])
  alpha_df <- bind_rows(lapply(seq_len(n_views), function(v) {
    data.frame(
      iter  = seq_len(n_iter),
      view  = v,
      alpha = res_gibbs$alpha_v[[v]],
      sigma = res_gibbs$sigma_v[[v]],
      tau   = res_gibbs$tau_v[[v]]
    )
  }))
} else if (is.matrix(res_gibbs$alpha_v)) {
  n_iter <- nrow(res_gibbs$alpha_v)
  n_views <- ncol(res_gibbs$alpha_v)
  alpha_df <- bind_rows(lapply(seq_len(n_views), function(v) {
    data.frame(
      iter  = seq_len(n_iter),
      view  = v,
      alpha = res_gibbs$alpha_v[, v],
      sigma = res_gibbs$sigma_v[, v],
      tau   = res_gibbs$tau_v[, v]
    )
  }))
} else {
  stop("alpha_v has an unsupported type.")
}

p_alpha_v <- ggplot(alpha_df, aes(iter, alpha, colour = factor(view))) +
  geom_line() +
  theme_minimal() +
  labs(title = "alpha_v by view", colour = "view")

p_sigma_v <- ggplot(alpha_df, aes(iter, sigma, colour = factor(view))) +
  geom_line() +
  theme_minimal() +
  labs(title = "sigma_v by view", colour = "view")

p_tau_v <- ggplot(alpha_df, aes(iter, tau, colour = factor(view))) +
  geom_line() +
  theme_minimal() +
  labs(title = "tau_v by view", colour = "view")

param_df <- tidyr::pivot_longer(
  alpha_df,
  cols      = c(alpha, sigma, tau),
  names_to  = "param",
  values_to = "value"
)

p_by_view <- ggplot(param_df, aes(iter, value, colour = param)) +
  geom_line() +
  theme_minimal() +
  facet_wrap(~view, scales = "free_y") +
  labs(title = "View-specific hyperparameters", colour = "param")

print(p_alpha_g)
print(p_sigma_g)
print(p_alpha_v)
print(p_sigma_v)
print(p_tau_v)
print(p_by_view)


window <- 500

final_alpha_global <- tail(res_gibbs$alpha_global, 1)
final_sigma_global <- tail(res_gibbs$sigma_global, 1)

mean_alpha_global <- mean(tail(res_gibbs$alpha_global, window))
sd_alpha_global <- sd(tail(res_gibbs$alpha_global, window))

mean_sigma_global <- mean(tail(res_gibbs$sigma_global, window))
sd_sigma_global <- sd(tail(res_gibbs$sigma_global, window))

cat("\n================ GLOBAL HYPERPARAMETERS ================\n")
cat("Final alpha_global:", final_alpha_global, "\n")
cat("Mean last", window, "alpha_global:", mean_alpha_global, " (sd =", sd_alpha_global, ")\n")
cat("Final sigma_global:", final_sigma_global, "\n")
cat("Mean last", window, "sigma_global:", mean_sigma_global, " (sd =", sd_sigma_global, ")\n")


# View-specific parameters
if (is.list(res_gibbs$alpha_v)) {
  n_views <- length(res_gibbs$alpha_v)
} else {
  n_views <- ncol(res_gibbs$alpha_v)
}


for (v in 1:n_views) {
  if (is.list(res_gibbs$alpha_v)) {
    a_vec <- res_gibbs$alpha_v[[v]]
    s_vec <- res_gibbs$sigma_v[[v]]
    t_vec <- res_gibbs$tau_v[[v]]
  } else {
    a_vec <- res_gibbs$alpha_v[, v]
    s_vec <- res_gibbs$sigma_v[, v]
    t_vec <- res_gibbs$tau_v[, v]
  }

  cat("\n--- View", v, "---\n")
  cat(
    "Final alpha_v:", tail(a_vec, 1),
    " | Mean last", window, "=", mean(tail(a_vec, window)), "\n"
  )
  cat(
    "Final sigma_v:", tail(s_vec, 1),
    " | Mean last", window, "=", mean(tail(s_vec, window)), "\n"
  )
  cat(
    "Final tau_v:", tail(t_vec, 1),
    " | Mean last", window, "=", mean(tail(t_vec, window)), "\n"
  )
}
