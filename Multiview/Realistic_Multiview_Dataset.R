# Realistic Multiview Dataset: Hospital Patient Risk Stratification
#
# USAGE: cd Multiview && Rscript Realistic_Multiview_Dataset.R
#
# Generates a multiview clustering dataset with 5 clinical views:
#   View 1: Vital Signs (3 clusters - Low/Moderate/High risk)
#   View 2: Metabolic Panel (2 clusters - Normal/Elevated)
#   View 3: Lifestyle (2 clusters - Active/Sedentary)
#   View 4: Mental Health (3 clusters - Good/Moderate/Poor)
#   View 5: Medical History (2 clusters - Simple/Complex)

library(dplyr)
library(ggplot2)
library(cowplot)
library(Rcpp)
library(mcclust)
library(mclust)
library(mcclust.ext)
library(gridExtra)
library(tidyr)
library(MASS)

Rcpp::sourceCpp("multiview_gibbs.cpp")

set.seed(42)

patient_counts <- c(
    healthy_active = 40,
    healthy_sedentary = 35,
    prediabetic = 30,
    cardiovascular = 35,
    chronic_complex = 30,
    mental_health = 30
)

n_total <- sum(patient_counts)
cat(sprintf("Total patients: %d\n", n_total))

generate_vital_signs_view <- function(patient_counts) {
    n_total <- sum(patient_counts)

    low_risk <- c(
        rnorm(patient_counts["healthy_active"], mean = 0, sd = 0.8),
        rnorm(patient_counts["mental_health"], mean = 0.5, sd = 0.9)
    )

    moderate_risk <- c(
        rnorm(patient_counts["healthy_sedentary"], mean = 3, sd = 0.8),
        rnorm(patient_counts["prediabetic"], mean = 3.5, sd = 0.9)
    )

    high_risk <- c(
        rnorm(patient_counts["cardiovascular"], mean = 6.5, sd = 0.8),
        rnorm(patient_counts["chronic_complex"], mean = 7, sd = 0.9)
    )

    data <- c(low_risk, moderate_risk, high_risk)

    labels <- c(
        rep(1, patient_counts["healthy_active"] + patient_counts["mental_health"]),
        rep(2, patient_counts["healthy_sedentary"] + patient_counts["prediabetic"]),
        rep(3, patient_counts["cardiovascular"] + patient_counts["chronic_complex"])
    )

    list(data = data, labels = labels, name = "Vital Signs", n_clusters = 3)
}

generate_metabolic_view <- function(patient_counts) {
    normal <- c(
        rnorm(patient_counts["healthy_active"], mean = -2, sd = 1.0),
        rnorm(patient_counts["healthy_sedentary"], mean = -1, sd = 1.1),
        rnorm(patient_counts["mental_health"], mean = -1.5, sd = 1.0),
        rnorm(patient_counts["cardiovascular"], mean = 0, sd = 1.2)
    )

    elevated <- c(
        rnorm(patient_counts["prediabetic"], mean = 4, sd = 1.0),
        rnorm(patient_counts["chronic_complex"], mean = 5, sd = 1.2)
    )

    data <- c(normal, elevated)

    labels <- c(
        rep(1, patient_counts["healthy_active"] + patient_counts["healthy_sedentary"] +
            patient_counts["mental_health"] + patient_counts["cardiovascular"]),
        rep(2, patient_counts["prediabetic"] + patient_counts["chronic_complex"])
    )

    list(data = data, labels = labels, name = "Metabolic Panel", n_clusters = 2)
}

generate_lifestyle_view <- function(patient_counts) {
    active <- c(
        rnorm(patient_counts["healthy_active"], mean = -3, sd = 0.9),
        rnorm(patient_counts["cardiovascular"], mean = -2, sd = 1.0)
    )

    sedentary <- c(
        rnorm(patient_counts["healthy_sedentary"], mean = 2.5, sd = 1.0),
        rnorm(patient_counts["prediabetic"], mean = 3, sd = 1.1),
        rnorm(patient_counts["chronic_complex"], mean = 3.5, sd = 1.0),
        rnorm(patient_counts["mental_health"], mean = 2, sd = 1.2)
    )

    data <- c(active, sedentary)

    labels <- c(
        rep(1, patient_counts["healthy_active"] + patient_counts["cardiovascular"]),
        rep(2, patient_counts["healthy_sedentary"] + patient_counts["prediabetic"] +
            patient_counts["chronic_complex"] + patient_counts["mental_health"])
    )

    list(data = data, labels = labels, name = "Lifestyle", n_clusters = 2)
}

generate_mental_health_view <- function(patient_counts) {
    good <- rnorm(patient_counts["healthy_active"], mean = -4, sd = 0.8)

    moderate <- c(
        rnorm(patient_counts["healthy_sedentary"], mean = 0, sd = 0.9),
        rnorm(patient_counts["cardiovascular"], mean = 0.5, sd = 1.0),
        rnorm(patient_counts["prediabetic"], mean = 0, sd = 0.9)
    )

    poor <- c(
        rnorm(patient_counts["chronic_complex"], mean = 4, sd = 0.9),
        rnorm(patient_counts["mental_health"], mean = 5, sd = 1.0)
    )

    data <- c(good, moderate, poor)

    labels <- c(
        rep(1, patient_counts["healthy_active"]),
        rep(2, patient_counts["healthy_sedentary"] + patient_counts["cardiovascular"] +
            patient_counts["prediabetic"]),
        rep(3, patient_counts["chronic_complex"] + patient_counts["mental_health"])
    )

    list(data = data, labels = labels, name = "Mental Health", n_clusters = 3)
}

generate_medical_history_view <- function(patient_counts) {
    simple <- c(
        rnorm(patient_counts["healthy_active"], mean = -3, sd = 0.8),
        rnorm(patient_counts["healthy_sedentary"], mean = -2.5, sd = 0.9),
        rnorm(patient_counts["mental_health"], mean = -2, sd = 1.0)
    )

    complex <- c(
        rnorm(patient_counts["prediabetic"], mean = 2, sd = 1.0),
        rnorm(patient_counts["cardiovascular"], mean = 3, sd = 0.9),
        rnorm(patient_counts["chronic_complex"], mean = 4, sd = 0.8)
    )

    data <- c(simple, complex)

    labels <- c(
        rep(1, patient_counts["healthy_active"] + patient_counts["healthy_sedentary"] +
            patient_counts["mental_health"]),
        rep(2, patient_counts["prediabetic"] + patient_counts["cardiovascular"] +
            patient_counts["chronic_complex"])
    )

    list(data = data, labels = labels, name = "Medical History", n_clusters = 2)
}

cat("\n=======================================================\n")
cat("GENERATING REALISTIC MULTIVIEW PATIENT DATASET\n")
cat("=======================================================\n\n")

view1 <- generate_vital_signs_view(patient_counts)
view2 <- generate_metabolic_view(patient_counts)
view3 <- generate_lifestyle_view(patient_counts)
view4 <- generate_mental_health_view(patient_counts)
view5 <- generate_medical_history_view(patient_counts)

data_views <- list(
    view1$data,
    view2$data,
    view3$data,
    view4$data,
    view5$data
)

true_labels <- list(
    view1$labels,
    view2$labels,
    view3$labels,
    view4$labels,
    view5$labels
)

view_info <- data.frame(
    view = 1:5,
    name = c(view1$name, view2$name, view3$name, view4$name, view5$name),
    n_clusters = c(
        view1$n_clusters, view2$n_clusters, view3$n_clusters,
        view4$n_clusters, view5$n_clusters
    )
)

cat("Dataset Summary:\n")
print(view_info)
cat(sprintf("\nTotal patients: %d\n", n_total))

patient_df <- data.frame(
    patient_id = 1:n_total,
    vital_signs = view1$data,
    metabolic = view2$data,
    lifestyle = view3$data,
    mental_health = view4$data,
    medical_history = view5$data,
    true_vital = view1$labels,
    true_metabolic = view2$labels,
    true_lifestyle = view3$labels,
    true_mental = view4$labels,
    true_history = view5$labels
)

p1 <- ggplot(patient_df, aes(
    x = vital_signs, y = metabolic,
    color = factor(true_vital), shape = factor(true_metabolic)
)) +
    geom_point(size = 2.5, alpha = 0.7) +
    labs(title = "Vitals vs Metabolic", color = "Vital Risk", shape = "Metabolic") +
    theme_minimal() +
    theme(legend.position = "bottom")

p2 <- ggplot(patient_df, aes(
    x = lifestyle, y = mental_health,
    color = factor(true_lifestyle), shape = factor(true_mental)
)) +
    geom_point(size = 2.5, alpha = 0.7) +
    labs(title = "Lifestyle vs Mental Health", color = "Lifestyle", shape = "Mental") +
    theme_minimal() +
    theme(legend.position = "bottom")

p3 <- ggplot(patient_df, aes(
    x = vital_signs, y = medical_history,
    color = factor(true_vital), shape = factor(true_history)
)) +
    geom_point(size = 2.5, alpha = 0.7) +
    labs(title = "Vitals vs Medical History", color = "Vital Risk", shape = "Complexity") +
    theme_minimal() +
    theme(legend.position = "bottom")

p_data_overview <- plot_grid(p1, p2, p3, ncol = 2, labels = "AUTO")
print(p_data_overview)

density_plots <- lapply(1:5, function(v) {
    df <- data.frame(value = data_views[[v]], cluster = factor(true_labels[[v]]))
    ggplot(df, aes(x = value, fill = cluster)) +
        geom_density(alpha = 0.5) +
        labs(title = view_info$name[v], x = "Value", fill = "Cluster") +
        theme_minimal()
})

p_densities <- plot_grid(plotlist = density_plots, ncol = 3)
print(p_densities)

cat("\n=======================================================\n")
cat("RUNNING MULTIVIEW GIBBS SAMPLER\n")
cat("=======================================================\n\n")

nsim <- 10000
burn_in <- 8000
thin <- 1

cat(sprintf(
    "MCMC settings: %d iterations, %d burn-in, thin = %d\n\n",
    nsim, burn_in, thin
))

res_gibbs <- run_gibbs_cpp(
    data_views = data_views,
    M = nsim,
    burn_in = burn_in,
    thin = thin
)

get_final_clusters <- function(res_gibbs) {
    last_iter_idx <- length(res_gibbs$table_of)
    raw_tables <- res_gibbs$table_of[[last_iter_idx]]
    raw_dishes <- res_gibbs$dish_of[[last_iter_idx]]

    tables_r_index <- raw_tables + 1
    n_customers <- length(tables_r_index)
    n_views <- length(raw_dishes)

    cluster_matrix <- matrix(NA, nrow = n_customers, ncol = n_views)
    for (v in 1:n_views) {
        cluster_matrix[, v] <- raw_dishes[[v]][tables_r_index]
    }

    cluster_matrix
}

predicted_clusters <- get_final_clusters(res_gibbs)

ari_scores <- sapply(1:5, function(v) {
    mcclust::arandi(predicted_clusters[, v], true_labels[[v]])
})

cat("\n=======================================================\n")
cat("EVALUATION RESULTS\n")
cat("=======================================================\n\n")

results_df <- data.frame(
    View = 1:5,
    Name = view_info$name,
    True_K = view_info$n_clusters,
    Pred_K = sapply(1:5, function(v) length(unique(predicted_clusters[, v]))),
    ARI = round(ari_scores, 4)
)

print(results_df)

cat(sprintf("\nMean ARI across all views: %.4f\n", mean(ari_scores)))

cat("\n--- Confusion Matrices ---\n")
for (v in 1:5) {
    cat(sprintf("\n%s (View %d):\n", view_info$name[v], v))
    print(table(Predicted = predicted_clusters[, v], Truth = true_labels[[v]]))
}

cat("\n=======================================================\n")
cat("CONVERGENCE DIAGNOSTICS\n")
cat("=======================================================\n\n")

df_global <- data.frame(
    iter = seq_along(res_gibbs$alpha_global),
    alpha_global = res_gibbs$alpha_global,
    sigma_global = res_gibbs$sigma_global
)

p_alpha <- ggplot(df_global, aes(iter, alpha_global)) +
    geom_line(alpha = 0.7) +
    labs(title = "α_global trace", x = "Iteration", y = "α") +
    theme_minimal()

p_sigma <- ggplot(df_global, aes(iter, sigma_global)) +
    geom_line(alpha = 0.7) +
    labs(title = "σ_global trace", x = "Iteration", y = "σ") +
    theme_minimal()

p_traces <- plot_grid(p_alpha, p_sigma, ncol = 2)
print(p_traces)

if (is.list(res_gibbs$alpha_v)) {
    n_views <- length(res_gibbs$alpha_v)
    n_iter <- length(res_gibbs$alpha_v[[1]])
    alpha_df <- bind_rows(lapply(1:n_views, function(v) {
        data.frame(
            iter = 1:n_iter,
            view = v,
            alpha = res_gibbs$alpha_v[[v]],
            sigma = res_gibbs$sigma_v[[v]],
            tau = res_gibbs$tau_v[[v]]
        )
    }))
} else {
    n_iter <- nrow(res_gibbs$alpha_v)
    n_views <- ncol(res_gibbs$alpha_v)
    alpha_df <- bind_rows(lapply(1:n_views, function(v) {
        data.frame(
            iter = 1:n_iter,
            view = v,
            alpha = res_gibbs$alpha_v[, v],
            sigma = res_gibbs$sigma_v[, v],
            tau = res_gibbs$tau_v[, v]
        )
    }))
}

p_tau <- ggplot(alpha_df, aes(iter, tau, color = factor(view))) +
    geom_line(alpha = 0.6) +
    labs(title = "τ_v (variance) by view", color = "View") +
    theme_minimal()

print(p_tau)

window <- 500
cat("\n--- Final Hyperparameter Estimates ---\n")
cat(sprintf(
    "α_global: %.3f (mean last %d: %.3f)\n",
    tail(res_gibbs$alpha_global, 1),
    window, mean(tail(res_gibbs$alpha_global, window))
))
cat(sprintf(
    "σ_global: %.3f (mean last %d: %.3f)\n",
    tail(res_gibbs$sigma_global, 1),
    window, mean(tail(res_gibbs$sigma_global, window))
))

ggsave("patient_data_overview.png", p_data_overview, width = 12, height = 8, dpi = 150)
ggsave("patient_view_densities.png", p_densities, width = 12, height = 6, dpi = 150)
ggsave("patient_hyperparameter_traces.png", p_traces, width = 10, height = 4, dpi = 150)

sink("patient_clustering_results.txt")
cat("=======================================================\n")
cat("REALISTIC MULTIVIEW CLUSTERING RESULTS\n")
cat("Hospital Patient Risk Stratification Scenario\n")
cat("=======================================================\n\n")
cat("View Definitions:\n")
print(view_info)
cat("\n\nPatient Archetype Distribution:\n")
print(patient_counts)
cat("\n\nClustering Performance:\n")
print(results_df)
cat(sprintf("\nOverall Mean ARI: %.4f\n", mean(ari_scores)))
sink()

cat("\n=======================================================\n")
cat("FILES SAVED:\n")
cat("=======================================================\n")
cat("  - patient_data_overview.png\n")
cat("  - patient_view_densities.png\n")
cat("  - patient_hyperparameter_traces.png\n")
cat("  - patient_clustering_results.txt\n")

cat("\n*** ANALYSIS COMPLETE ***\n")
