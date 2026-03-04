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

Rcpp::sourceCpp("../src/multiview_gibbs.cpp")

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

# =============================================================================
# TRANSITION PROBABILITY MATRICES (View 1 -> Views 2-5)
# =============================================================================
# Each matrix: rows = View 1 clusters (3), cols = target view clusters
# Row i, Col j = P(Target View = j | View 1 = i)
# Rows must sum to 1.0

# View 1 (3 clusters: Low/Moderate/High risk) -> View 2 (2 clusters: Normal/Elevated)
transition_v1_to_v2 <- matrix(c(
    0.85, 0.15, # V1-Low risk     -> 85% Normal, 15% Elevated
    0.50, 0.50, # V1-Moderate     -> 50% Normal, 50% Elevated
    0.20, 0.80 # V1-High risk    -> 20% Normal, 80% Elevated
), nrow = 3, byrow = TRUE)

# View 1 -> View 3 (2 clusters: Active/Sedentary)
transition_v1_to_v3 <- matrix(c(
    0.75, 0.25, # V1-Low risk     -> 75% Active, 25% Sedentary
    0.40, 0.60, # V1-Moderate     -> 40% Active, 60% Sedentary
    0.30, 0.70 # V1-High risk    -> 30% Active, 70% Sedentary
), nrow = 3, byrow = TRUE)

# View 1 -> View 4 (3 clusters: Good/Moderate/Poor mental health)
transition_v1_to_v4 <- matrix(c(
    0.70, 0.25, 0.05, # V1-Low risk     -> 70% Good, 25% Moderate, 5% Poor
    0.20, 0.60, 0.20, # V1-Moderate     -> 20% Good, 60% Moderate, 20% Poor
    0.05, 0.35, 0.60 # V1-High risk    -> 5% Good, 35% Moderate, 60% Poor
), nrow = 3, byrow = TRUE)

# View 1 -> View 5 (2 clusters: Simple/Complex history)
transition_v1_to_v5 <- matrix(c(
    0.80, 0.20, # V1-Low risk     -> 80% Simple, 20% Complex
    0.45, 0.55, # V1-Moderate     -> 45% Simple, 55% Complex
    0.15, 0.85 # V1-High risk    -> 15% Simple, 85% Complex
), nrow = 3, byrow = TRUE)

# Helper function: assign cluster based on View 1 cluster and transition matrix
assign_dependent_cluster <- function(v1_labels, transition_matrix) {
    n <- length(v1_labels)
    n_target_clusters <- ncol(transition_matrix)
    target_labels <- integer(n)

    for (i in 1:n) {
        v1_cluster <- v1_labels[i]
        probs <- transition_matrix[v1_cluster, ]
        target_labels[i] <- sample(1:n_target_clusters, size = 1, prob = probs)
    }

    return(target_labels)
}

generate_vital_signs_view <- function(patient_counts) {
    # CLOSER CLUSTERS: means at 0, 1.5, 3 (was 0, 3, 6.5)
    low_risk <- c(
        rnorm(patient_counts["healthy_active"], mean = 0, sd = 0.8),
        rnorm(patient_counts["mental_health"], mean = 0.25, sd = 0.9)
    )

    moderate_risk <- c(
        rnorm(patient_counts["healthy_sedentary"], mean = 1.5, sd = 0.8),
        rnorm(patient_counts["prediabetic"], mean = 1.75, sd = 0.9)
    )

    high_risk <- c(
        rnorm(patient_counts["cardiovascular"], mean = 3.0, sd = 0.8),
        rnorm(patient_counts["chronic_complex"], mean = 3.25, sd = 0.9)
    )

    data <- c(low_risk, moderate_risk, high_risk)

    labels <- c(
        rep(1, patient_counts["healthy_active"] + patient_counts["mental_health"]),
        rep(2, patient_counts["healthy_sedentary"] + patient_counts["prediabetic"]),
        rep(3, patient_counts["cardiovascular"] + patient_counts["chronic_complex"])
    )

    list(data = data, labels = labels, name = "Vital Signs", n_clusters = 3)
}

# Cluster parameters for each view (used with dependent generation)
# CLOSER CLUSTERS: separation reduced
metabolic_params <- list(
    cluster_1 = list(mean = -0.75, sd = 1.1), # Normal (was -1.5)
    cluster_2 = list(mean = 2.0, sd = 1.1) # Elevated (was 4.5)
)

generate_metabolic_view <- function(v1_labels, transition_matrix) {
    # Assign clusters based on View 1 dependency
    labels <- assign_dependent_cluster(v1_labels, transition_matrix)
    n <- length(labels)

    # Generate data based on assigned clusters
    data <- numeric(n)
    for (i in 1:n) {
        cluster <- labels[i]
        params <- metabolic_params[[paste0("cluster_", cluster)]]
        data[i] <- rnorm(1, mean = params$mean, sd = params$sd)
    }

    list(data = data, labels = labels, name = "Metabolic Panel", n_clusters = 2)
}

# CLOSER CLUSTERS
lifestyle_params <- list(
    cluster_1 = list(mean = -1.0, sd = 1.0), # Active (was -2.5)
    cluster_2 = list(mean = 1.25, sd = 1.0) # Sedentary (was 2.75)
)

generate_lifestyle_view <- function(v1_labels, transition_matrix) {
    labels <- assign_dependent_cluster(v1_labels, transition_matrix)
    n <- length(labels)

    data <- numeric(n)
    for (i in 1:n) {
        cluster <- labels[i]
        params <- lifestyle_params[[paste0("cluster_", cluster)]]
        data[i] <- rnorm(1, mean = params$mean, sd = params$sd)
    }

    list(data = data, labels = labels, name = "Lifestyle", n_clusters = 2)
}

# CLOSER CLUSTERS: means at -1.5, 0, 1.5 (was -4, 0, 4.5)
mental_health_params <- list(
    cluster_1 = list(mean = -1.5, sd = 0.8), # Good
    cluster_2 = list(mean = 0, sd = 0.9), # Moderate
    cluster_3 = list(mean = 1.5, sd = 0.9) # Poor
)

generate_mental_health_view <- function(v1_labels, transition_matrix) {
    labels <- assign_dependent_cluster(v1_labels, transition_matrix)
    n <- length(labels)

    data <- numeric(n)
    for (i in 1:n) {
        cluster <- labels[i]
        params <- mental_health_params[[paste0("cluster_", cluster)]]
        data[i] <- rnorm(1, mean = params$mean, sd = params$sd)
    }

    list(data = data, labels = labels, name = "Mental Health", n_clusters = 3)
}

# CLOSER CLUSTERS
medical_history_params <- list(
    cluster_1 = list(mean = -1.0, sd = 0.9), # Simple (was -2.5)
    cluster_2 = list(mean = 1.25, sd = 0.9) # Complex (was 3)
)

generate_medical_history_view <- function(v1_labels, transition_matrix) {
    labels <- assign_dependent_cluster(v1_labels, transition_matrix)
    n <- length(labels)

    data <- numeric(n)
    for (i in 1:n) {
        cluster <- labels[i]
        params <- medical_history_params[[paste0("cluster_", cluster)]]
        data[i] <- rnorm(1, mean = params$mean, sd = params$sd)
    }

    list(data = data, labels = labels, name = "Medical History", n_clusters = 2)
}

cat("\n=======================================================\n")
cat("GENERATING REALISTIC MULTIVIEW PATIENT DATASET\n")
cat("WITH VIEW DEPENDENCIES (View 1 -> Views 2-5)\n")
cat("=======================================================\n\n")

# Generate View 1 first (independent base view)
view1 <- generate_vital_signs_view(patient_counts)
cat(sprintf("View 1 generated: %d patients in %d clusters\n", length(view1$data), view1$n_clusters))

# Generate Views 2-5 with dependencies on View 1
view2 <- generate_metabolic_view(view1$labels, transition_v1_to_v2)
view3 <- generate_lifestyle_view(view1$labels, transition_v1_to_v3)
view4 <- generate_mental_health_view(view1$labels, transition_v1_to_v4)
view5 <- generate_medical_history_view(view1$labels, transition_v1_to_v5)

cat("\nDependency-based cluster distributions:\n")
cat("\nView 2 (Metabolic) cluster counts:\n")
print(table(view2$labels))
cat("View 3 (Lifestyle) cluster counts:\n")
print(table(view3$labels))
cat("View 4 (Mental Health) cluster counts:\n")
print(table(view4$labels))
cat("View 5 (Medical History) cluster counts:\n")
print(table(view5$labels))

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

# Save dataset to CSV
csv_path <- "../dataset/dependent_multiview_dataset.csv"
write.csv(patient_df, csv_path, row.names = FALSE)
cat(sprintf("\nDataset saved to: %s\n", csv_path))

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

# =============================================================================
# CONSENSUS CLUSTERING WITH minVI
# =============================================================================
# Build posterior similarity matrices and find optimal partition using minVI

cat("\n=======================================================\n")
cat("COMPUTING CONSENSUS CLUSTERING (minVI)\n")
cat("=======================================================\n\n")

n_samples <- length(res_gibbs$table_of)
n_patients <- length(res_gibbs$table_of[[1]])
n_views <- 5

cat(sprintf("Using %d MCMC samples for consensus\n", n_samples))

# Function to get cluster assignments for all samples
get_cluster_matrix_per_view <- function(res_gibbs, view_idx) {
    n_samples <- length(res_gibbs$table_of)
    n_patients <- length(res_gibbs$table_of[[1]])

    cluster_mat <- matrix(NA, nrow = n_samples, ncol = n_patients)

    for (s in 1:n_samples) {
        tables <- res_gibbs$table_of[[s]] + 1 # Convert to 1-indexed
        dishes <- res_gibbs$dish_of[[s]][[view_idx]]
        raw_labels <- dishes[tables]

        # Remap to consecutive integers 1:K (required by comp.psm)
        unique_labels <- unique(raw_labels)
        label_map <- setNames(seq_along(unique_labels), unique_labels)
        cluster_mat[s, ] <- label_map[as.character(raw_labels)]
    }

    cluster_mat
}

# Compute minVI consensus for each view
predicted_clusters <- matrix(NA, nrow = n_patients, ncol = n_views)

for (v in 1:n_views) {
    cat(sprintf("Processing View %d (%s)...\n", v, view_info$name[v]))

    # Get cluster assignments across all samples
    cluster_samples <- get_cluster_matrix_per_view(res_gibbs, v)

    # Compute posterior similarity matrix (PSM)
    psm <- mcclust::comp.psm(cluster_samples)

    # Find optimal partition using minVI (avg method)
    minvi_result <- mcclust.ext::minVI(psm, method = "avg")

    predicted_clusters[, v] <- minvi_result$cl

    cat(sprintf("  Found %d clusters (minVI)\n", length(unique(minvi_result$cl))))
}

cat("\nConsensus clustering complete.\n")

# Compute ARI scores
ari_scores <- sapply(1:5, function(v) {
    mcclust::arandi(predicted_clusters[, v], true_labels[[v]])
})

cat("\n=======================================================\n")
cat("EVALUATION RESULTS (minVI Consensus)\n")
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
