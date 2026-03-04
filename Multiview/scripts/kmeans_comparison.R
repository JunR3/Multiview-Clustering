library(dplyr)
library(mclust)

data <- read.csv("../dataset/dependent_multiview_dataset.csv")

cat("=======================================================\n")
cat("K-MEANS VS MULTIVIEW CLUSTERING COMPARISON\n")
cat("=======================================================\n\n")

cat(sprintf("Dataset: %d patients, %d views\n", nrow(data), 5))

views <- list(
    vital_signs = data$vital_signs,
    metabolic = data$metabolic,
    lifestyle = data$lifestyle,
    mental_health = data$mental_health,
    medical_history = data$medical_history
)

true_labels <- list(
    vital_signs = data$true_vital,
    metabolic = data$true_metabolic,
    lifestyle = data$true_lifestyle,
    mental_health = data$true_mental,
    medical_history = data$true_history
)

view_names <- c("Vital Signs", "Metabolic Panel", "Lifestyle", "Mental Health", "Medical History")
true_k <- c(3, 2, 2, 3, 2)

cat("\n--- APPROACH 1: K-means + Silhouette K Selection (1D) ---\n\n")

library(cluster)

find_optimal_k_silhouette <- function(data, k_range = 2:6) {
    data_mat <- matrix(data, ncol = 1)

    sil_scores <- sapply(k_range, function(k) {
        km <- kmeans(data_mat, centers = k, nstart = 25)
        sil <- silhouette(km$cluster, dist(data_mat))
        mean(sil[, 3])
    })

    best_k <- k_range[which.max(sil_scores)]
    list(best_k = best_k, scores = sil_scores)
}

kmeans_1d_results <- data.frame(
    View = view_names,
    True_K = true_k,
    Selected_K = integer(5),
    ARI = numeric(5)
)

for (i in 1:5) {
    view_data <- views[[i]]

    opt <- find_optimal_k_silhouette(view_data, k_range = 2:6)
    selected_k <- opt$best_k

    km_result <- kmeans(matrix(view_data, ncol = 1), centers = selected_k, nstart = 25)
    ari <- adjustedRandIndex(km_result$cluster, true_labels[[i]])

    kmeans_1d_results$Selected_K[i] <- selected_k
    kmeans_1d_results$ARI[i] <- round(ari, 4)

    cat(sprintf(
        "%s: Silhouette selected K=%d (true K=%d)\n",
        view_names[i], selected_k, true_k[i]
    ))
}

cat("\n")
print(kmeans_1d_results)
cat(sprintf("\nMean ARI (K-means + Silhouette): %.4f\n", mean(kmeans_1d_results$ARI)))

cat("\n--- APPROACH 2: K-means on concatenated views (5D data) ---\n\n")

combined_data <- cbind(
    data$vital_signs,
    data$metabolic,
    data$lifestyle,
    data$mental_health,
    data$medical_history
)

combined_scaled <- scale(combined_data)

cat("Testing different K values on combined 5D data:\n\n")

kmeans_combined_results <- data.frame(
    View = view_names,
    True_K = true_k,
    ARI_K2 = numeric(5),
    ARI_K3 = numeric(5),
    ARI_K4 = numeric(5),
    ARI_K5 = numeric(5)
)

for (k in 2:5) {
    km_combined <- kmeans(combined_scaled, centers = k, nstart = 25)

    for (i in 1:5) {
        ari <- adjustedRandIndex(km_combined$cluster, true_labels[[i]])
        col_name <- paste0("ARI_K", k)
        kmeans_combined_results[i, col_name] <- round(ari, 4)
    }
}

print(kmeans_combined_results)

cat("\nBest combined K-means (mean ARI per K):\n")
for (k in 2:5) {
    col_name <- paste0("ARI_K", k)
    mean_ari <- mean(kmeans_combined_results[[col_name]])
    cat(sprintf("  K=%d: Mean ARI = %.4f\n", k, mean_ari))
}

cat("\n--- APPROACH 3: Multiview Clustering (from previous run) ---\n\n")

multiview_results <- data.frame(
    View = view_names,
    True_K = true_k,
    Pred_K = c(2, 2, 2, 2, 2),
    ARI = c(0.4086, 0.6707, 0.2668, 0.1964, 0.3215)
)

print(multiview_results)
cat(sprintf("\nMean ARI (Multiview Gibbs): %.4f\n", mean(multiview_results$ARI)))

cat("\n=======================================================\n")
cat("SUMMARY COMPARISON\n")
cat("=======================================================\n\n")

best_k_combined <- which.max(sapply(2:5, function(k) {
    mean(kmeans_combined_results[[paste0("ARI_K", k)]])
})) + 1

comparison <- data.frame(
    Method = c(
        "K-means + Silhouette (1D)",
        sprintf("K-means (combined K=%d)", best_k_combined),
        "Multiview Gibbs + minVI"
    ),
    Mean_ARI = c(
        mean(kmeans_1d_results$ARI),
        mean(kmeans_combined_results[[paste0("ARI_K", best_k_combined)]]),
        mean(multiview_results$ARI)
    )
)

comparison$Mean_ARI <- round(comparison$Mean_ARI, 4)

print(comparison)

cat("\n--- Per-View ARI Comparison ---\n\n")
per_view_comparison <- data.frame(
    View = view_names,
    Kmeans_1D = kmeans_1d_results$ARI,
    Kmeans_5D = kmeans_combined_results[[paste0("ARI_K", best_k_combined)]],
    Multiview = multiview_results$ARI
)
print(per_view_comparison)

cat("\n--- Analysis ---\n")
cat(sprintf(
    "Best overall method: %s (Mean ARI = %.4f)\n",
    comparison$Method[which.max(comparison$Mean_ARI)],
    max(comparison$Mean_ARI)
))

cat("\nPer-view winners:\n")
for (i in 1:5) {
    methods <- c("K-means 1D", "K-means 5D", "Multiview")
    aris <- c(
        per_view_comparison$Kmeans_1D[i],
        per_view_comparison$Kmeans_5D[i],
        per_view_comparison$Multiview[i]
    )
    winner <- methods[which.max(aris)]
    cat(sprintf("  %s: %s (ARI = %.4f)\n", view_names[i], winner, max(aris)))
}

cat("\n*** COMPARISON COMPLETE ***\n")
