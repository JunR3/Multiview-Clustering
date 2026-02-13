library(dplyr)
library(ggplot2)
library(cowplot)
library(Rcpp)
library(mcclust)
library(mclust)
library(mcclust.ext)
library(gridExtra)
library(clue)
library(tidyr)

Rcpp::sourceCpp("MultiView/lib/multiview_gibbs.cpp")

set.seed(2024)

df_raw <- read.csv("dataset/CPP_dataset.csv", check.names = FALSE)

if (colnames(df_raw)[1] == "") {
    colnames(df_raw)[1] <- "Index"
}

df_filtered <- df_raw %>% filter(gest < 42)

n_sample <- 1000
if (nrow(df_filtered) > n_sample) {
    df_sample <- df_filtered %>% sample_n(n_sample)
} else {
    df_sample <- df_filtered
}

print(paste("Data loaded. Filtered n =", nrow(df_filtered), ". Sampled n =", nrow(df_sample)))

y_true <- as.factor(df_sample$smoke)

v1_raw <- df_sample$weight
v2_raw <- df_sample$gest
v3_raw <- df_sample$dde

v1_scaled <- as.vector(scale(v1_raw))
v2_scaled <- as.vector(scale(v2_raw))

data_views <- list(
    v1_scaled,
    v2_scaled
)

view_names <- c("Weight", "Gestation")

true_clust_list <- list(
    as.integer(y_true),
    as.integer(y_true)
)

nsim <- 100000
burn_in <- 10000
thin <- 5

print("Running CPP Simulation (Dombowsky Replication)...")
print(paste("M =", nsim, ", Burn-in =", burn_in, ", Thin =", thin))

res_gibbs <- run_gibbs_cpp(
    data_views = data_views,
    M          = nsim,
    burn_in    = burn_in,
    thin       = thin
)

get_final_clusters <- function(res_gibbs, n_views) {
    last_iter_idx <- length(res_gibbs$table_of)
    raw_tables <- res_gibbs$table_of[[last_iter_idx]]
    raw_dishes <- res_gibbs$dish_of[[last_iter_idx]]
    tables_r_index <- raw_tables + 1
    n_customers <- length(tables_r_index)

    cluster_matrix <- matrix(NA, nrow = n_customers, ncol = n_views)
    colnames(cluster_matrix) <- paste0("View_", 1:n_views)

    for (v in 1:n_views) {
        dishes_for_view <- raw_dishes[[v]]
        cluster_matrix[, v] <- dishes_for_view[tables_r_index]
    }
    return(cluster_matrix)
}

my_clusters <- get_final_clusters(res_gibbs, length(data_views))

print("\n--- Computing Consensus Partition (minVI) ---")

compute_psm <- function(saved_dish_of, saved_table_of, view_idx) {
    n_saved <- length(saved_table_of)
    n_cust <- length(saved_table_of[[1]])

    psm <- matrix(0, nrow = n_cust, ncol = n_cust)

    for (iter in 1:n_saved) {
        tables <- saved_table_of[[iter]] + 1
        dishes <- saved_dish_of[[iter]][[view_idx]]

        clust <- dishes[tables]

        for (i in 1:(n_cust - 1)) {
            for (j in (i + 1):n_cust) {
                if (clust[i] == clust[j]) {
                    psm[i, j] <- psm[i, j] + 1
                    psm[j, i] <- psm[j, i] + 1
                }
            }
        }
    }

    psm <- psm / n_saved
    diag(psm) <- 1
    return(psm)
}

consensus_clusters <- matrix(NA, nrow = nrow(df_sample), ncol = length(data_views))
colnames(consensus_clusters) <- paste0("View_", 1:length(data_views))

for (v in seq_along(data_views)) {
    cat(paste("Computing PSM for", view_names[v], "...\n"))
    psm_v <- compute_psm(res_gibbs$dish_of, res_gibbs$table_of, v)

    minVI_result <- mcclust.ext::minVI(psm_v, method = "avg", max.k = 10)
    consensus_clusters[, v] <- minVI_result$cl

    cat(paste(view_names[v], ": minVI found", length(unique(minVI_result$cl)), "clusters\n"))
}

print("\n--- Final Iteration Cluster Counts per View ---")
for (v in seq_along(data_views)) {
    n_clust <- length(unique(my_clusters[, v]))
    print(paste(view_names[v], ":", n_clust, "clusters found."))
}

print("\n--- Consensus (minVI) Cluster Counts per View ---")
for (v in seq_along(data_views)) {
    n_clust <- length(unique(consensus_clusters[, v]))
    print(paste(view_names[v], ":", n_clust, "clusters found."))
}

print("\n--- ARI Scores vs Smoking Status (Consensus) ---")
ari_scores_consensus <- sapply(seq_along(data_views), function(v) mcclust::arandi(consensus_clusters[, v], true_clust_list[[v]]))
names(ari_scores_consensus) <- view_names
print(ari_scores_consensus)

print("\n--- ARI Scores vs Smoking Status (Final Iteration) ---")
ari_scores <- sapply(seq_along(data_views), function(v) mcclust::arandi(my_clusters[, v], true_clust_list[[v]]))
names(ari_scores) <- view_names
print(ari_scores)

print("\n--- Confusion Matrices (Consensus vs Smoking) ---")
for (v in seq_along(data_views)) {
    cat(paste("\nView:", view_names[v], "\n"))
    print(table(Predicted = consensus_clusters[, v], Smoke = true_clust_list[[v]]))
}

df_global <- data.frame(
    iter = seq_along(res_gibbs$alpha_global),
    alpha_global = res_gibbs$alpha_global,
    sigma_global = res_gibbs$sigma_global
)

p_alpha_g <- ggplot(df_global, aes(iter, alpha_global)) +
    geom_line() +
    theme_minimal() +
    labs(title = "Alpha Global Trace")
p_sigma_g <- ggplot(df_global, aes(iter, sigma_global)) +
    geom_line() +
    theme_minimal() +
    labs(title = "Sigma Global Trace")

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
} else {
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
}

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


df_plot <- data.frame(
    Weight = v1_scaled,
    Gestation = v2_scaled,
    Pred_Weight = as.factor(my_clusters[, 1]),
    Pred_Gest = as.factor(my_clusters[, 2]),
    Smoke = as.factor(df_sample$smoke)
)

p_scatter_weight <- ggplot(df_plot, aes(x = Weight, y = Gestation, color = Pred_Weight)) +
    geom_point(alpha = 0.6) +
    theme_bw() +
    labs(title = "Weight vs Gestation", subtitle = "Colored by Weight Clusters")

p_scatter_gest <- ggplot(df_plot, aes(x = Weight, y = Gestation, color = Pred_Gest)) +
    geom_point(alpha = 0.6) +
    theme_bw() +
    labs(title = "Weight vs Gestation", subtitle = "Colored by Gestation Clusters")

p_scatter_smoke <- ggplot(df_plot, aes(x = Weight, y = Gestation, color = Smoke)) +
    geom_point(alpha = 0.6) +
    theme_bw() +
    labs(title = "Weight vs Gestation", subtitle = "Colored by Smoking (Ground Truth)")


p_interaction <- ggplot(df_plot, aes(
    x = Weight, y = Gestation,
    shape = Pred_Weight,
    color = Pred_Gest
)) +
    geom_point(size = 3.5, alpha = 0.8) +
    labs(
        title = "Predicted Interaction: Weight vs Gestation",
        subtitle = "Shape: Weight Clusters | Color: Gestation Clusters",
        x = "Weight (Scaled)",
        y = "Gestation (Scaled)",
        shape = "Weight Clust",
        color = "Gest Clust"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

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

if (is.list(res_gibbs$alpha_v)) {
    n_views <- length(res_gibbs$alpha_v)
} else {
    n_views <- ncol(res_gibbs$alpha_v)
}

for (v in seq_len(n_views)) {
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

pdf("CPP_Simulation_Plots.pdf", width = 10, height = 7)
print(p_alpha_g)
print(p_sigma_g)
print(p_by_view)
print(p_scatter_weight)
print(p_scatter_gest)
print(p_interaction)
print(p_scatter_smoke)
dev.off()
print("Plots saved to CPP_Simulation_Plots.pdf")
