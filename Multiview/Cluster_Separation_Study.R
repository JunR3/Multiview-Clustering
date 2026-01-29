library(dplyr)
library(ggplot2)
library(cowplot)
library(Rcpp)
library(mcclust)
library(mclust)
library(mcclust.ext)
library(gridExtra)
library(tidyr)

Rcpp::sourceCpp("multiview_gibbs.cpp")

set.seed(2024)

generate_two_cluster_view <- function(n_per_cluster = 100, separation = 3.0, sd = 1.0) {
    mean1 <- -separation * sd / 2
    mean2 <- separation * sd / 2

    data <- c(
        rnorm(n_per_cluster, mean = mean1, sd = sd),
        rnorm(n_per_cluster, mean = mean2, sd = sd)
    )

    labels <- c(rep(1, n_per_cluster), rep(2, n_per_cluster))

    list(data = data, labels = labels)
}

generate_three_cluster_view <- function(n_cluster1 = 50, n_cluster2 = 100,
                                        n_cluster3 = 50, separation = 3.0, sd = 1.0) {
    mean1 <- -separation * sd
    mean2 <- 0
    mean3 <- separation * sd

    data <- c(
        rnorm(n_cluster1, mean = mean1, sd = sd),
        rnorm(n_cluster2, mean = mean2, sd = sd),
        rnorm(n_cluster3, mean = mean3, sd = sd)
    )

    labels <- c(rep(1, n_cluster1), rep(2, n_cluster2), rep(3, n_cluster3))

    list(data = data, labels = labels)
}

generate_multiview_dataset <- function(n_obs = 200, separation = 3.0,
                                       n_views = 5, view2_three_cluster = TRUE) {
    n_per_cluster <- n_obs / 2

    data_views <- list()
    true_labels <- list()

    for (v in 1:n_views) {
        if (v == 2 && view2_three_cluster) {
            result <- generate_three_cluster_view(
                n_cluster1 = n_obs / 4,
                n_cluster2 = n_obs / 2,
                n_cluster3 = n_obs / 4,
                separation = separation,
                sd = 1.0
            )
        } else {
            result <- generate_two_cluster_view(
                n_per_cluster = n_per_cluster,
                separation = separation,
                sd = 1.0
            )
        }

        data_views[[v]] <- result$data
        true_labels[[v]] <- result$labels
    }

    list(
        data_views = data_views,
        true_labels = true_labels,
        n_obs = n_obs,
        separation = separation,
        n_views = n_views
    )
}

get_final_clusters <- function(res_gibbs, n_samples = 100) {
    cat("Getting final clusters\n")
    # Use last n_samples iterations to build posterior similarity matrix
    n_iters <- length(res_gibbs$table_of)
    start_idx <- max(1, n_iters - n_samples + 1)
    sample_indices <- start_idx:n_iters

    # Get dimensions from first sample
    first_tables <- res_gibbs$table_of[[sample_indices[1]]]
    n_customers <- length(first_tables)
    n_views <- length(res_gibbs$dish_of[[sample_indices[1]]])

    cluster_matrix <- matrix(NA, nrow = n_customers, ncol = n_views)
    colnames(cluster_matrix) <- paste0("View_", 1:n_views)

    # For each view, build PSM and find optimal partition via minVI
    cat("Building PSM and finding optimal partition via minVI\n")
    for (v in 1:n_views) {
        # Collect cluster assignments across samples for this view
        cluster_samples <- matrix(NA, nrow = length(sample_indices), ncol = n_customers)

        for (i in seq_along(sample_indices)) {
            cat("Processing sample", i, "\n")
            idx <- sample_indices[i]
            tables_r_index <- res_gibbs$table_of[[idx]] + 1
            dishes_for_view <- res_gibbs$dish_of[[idx]][[v]]
            cluster_samples[i, ] <- dishes_for_view[tables_r_index]
        }

        # Normalize labels to be consecutive integers starting from 1
        # comp.psm requires labels in 1:nobs format
        for (i in 1:nrow(cluster_samples)) {
            cluster_samples[i, ] <- as.integer(as.factor(cluster_samples[i, ]))
        }

        # Build posterior similarity matrix
        cat("Building PSM\n")
        psm <- mcclust::comp.psm(cluster_samples)

        # Find partition minimizing Variation of Information
        cat("Finding partition minimizing Variation of Information\n")
        vi_result <- mcclust.ext::minVI(psm, method = "avg")
        cat("Found partition minimizing Variation of Information\n")
        cluster_matrix[, v] <- vi_result$cl
        cat("Done finding partition minimizing Variation of Information\n")
    }

    cat("Done getting final clusters\n")
    cluster_matrix
}

compute_ari_scores <- function(predicted, true_labels) {
    n_views <- ncol(predicted)
    ari_scores <- sapply(1:n_views, function(v) {
        mcclust::arandi(predicted[, v], true_labels[[v]])
    })
    names(ari_scores) <- paste0("View_", 1:n_views)
    ari_scores
}

run_separation_sweep <- function(separations = seq(6, 0.5, by = -0.5),
                                 n_obs = 200,
                                 nsim = 5000,
                                 burn_in = 4000,
                                 thin = 1) {
    results <- data.frame()

    cat("\n=======================================================\n")
    cat("CLUSTER SEPARATION STUDY\n")
    cat("=======================================================\n")
    cat(sprintf(
        "Testing %d separation values from %.1f to %.1f sigma\n",
        length(separations), max(separations), min(separations)
    ))
    cat(sprintf("Sample size: %d, MCMC iterations: %d\n\n", n_obs, nsim))

    for (sep in separations) {
        cat(sprintf("\n--- Testing separation = %.1f sigma ---\n", sep))

        dataset <- generate_multiview_dataset(
            n_obs = n_obs,
            separation = sep,
            n_views = 5,
            view2_three_cluster = TRUE
        )

        res_gibbs <- run_gibbs_cpp(
            data_views = dataset$data_views,
            M = nsim,
            burn_in = burn_in,
            thin = thin
        )

        predicted <- get_final_clusters(res_gibbs)
        cat("Got final clusters\n")
        ari_scores <- compute_ari_scores(predicted, dataset$true_labels)

        n_clusters_found <- sapply(1:5, function(v) length(unique(predicted[, v])))
        n_clusters_true <- sapply(dataset$true_labels, function(x) length(unique(x)))

        cat(sprintf("  ARI scores: %s\n", paste(round(ari_scores, 3), collapse = ", ")))
        cat(sprintf(
            "  Clusters found: %s (true: %s)\n",
            paste(n_clusters_found, collapse = ", "),
            paste(n_clusters_true, collapse = ", ")
        ))

        for (v in 1:5) {
            results <- rbind(results, data.frame(
                separation = sep,
                view = v,
                ari = ari_scores[v],
                n_clusters_found = n_clusters_found[v],
                n_clusters_true = n_clusters_true[v],
                cluster_match = n_clusters_found[v] == n_clusters_true[v]
            ))
        }
    }

    results
}

plot_ari_curve <- function(results) {
    summary_df <- results %>%
        group_by(separation) %>%
        summarize(
            mean_ari = mean(ari),
            sd_ari = sd(ari),
            min_ari = min(ari),
            max_ari = max(ari),
            .groups = "drop"
        )

    threshold_row <- summary_df %>%
        filter(mean_ari < 0.9) %>%
        slice_max(separation, n = 1)

    threshold <- if (nrow(threshold_row) > 0) threshold_row$separation else NA

    p <- ggplot(summary_df, aes(x = separation, y = mean_ari)) +
        geom_ribbon(aes(ymin = min_ari, ymax = max_ari),
            fill = "steelblue", alpha = 0.2
        ) +
        geom_line(linewidth = 1.2, color = "steelblue") +
        geom_point(size = 3, color = "steelblue") +
        geom_hline(yintercept = 0.9, linetype = "dashed", color = "red", linewidth = 0.8) +
        geom_hline(yintercept = 0.8, linetype = "dotted", color = "orange", linewidth = 0.8) +
        {
            if (!is.na(threshold)) {
                geom_vline(
                    xintercept = threshold,
                    linetype = "dashed", color = "darkred"
                )
            }
        } +
        annotate("text",
            x = max(results$separation) - 0.5, y = 0.92,
            label = "ARI = 0.9 threshold", color = "red", hjust = 1
        ) +
        labs(
            title = "Cluster Separation Threshold Analysis",
            subtitle = "Mean ARI across views (shaded = min/max range)",
            x = "Cluster Separation (sigma units)",
            y = "Adjusted Rand Index (ARI)"
        ) +
        scale_x_reverse() +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
        theme_minimal(base_size = 14) +
        theme(
            plot.title = element_text(face = "bold"),
            panel.grid.minor = element_blank()
        )

    if (!is.na(threshold)) {
        cat(sprintf(
            "\n*** CRITICAL THRESHOLD: ARI drops below 0.9 at separation = %.1f sigma ***\n",
            threshold
        ))
    }

    p
}

plot_ari_by_view <- function(results) {
    results$view_label <- paste("View", results$view)
    results$view_label[results$view == 2] <- "View 2 (3 clusters)"

    ggplot(results, aes(x = separation, y = ari, color = view_label)) +
        geom_line(linewidth = 1) +
        geom_point(size = 2) +
        geom_hline(yintercept = 0.9, linetype = "dashed", color = "gray50") +
        labs(
            title = "ARI by View Across Separations",
            x = "Cluster Separation (sigma units)",
            y = "Adjusted Rand Index (ARI)",
            color = "View"
        ) +
        scale_x_reverse() +
        scale_y_continuous(limits = c(0, 1)) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "bottom")
}

plot_scatter_panels <- function(separations = c(6, 3, 2, 1), n_obs = 200) {
    plots <- list()

    for (i in seq_along(separations)) {
        sep <- separations[i]
        dataset <- generate_multiview_dataset(n_obs = n_obs, separation = sep)

        df <- data.frame(
            x = dataset$data_views[[1]],
            y = dataset$data_views[[2]],
            cluster_v1 = as.factor(dataset$true_labels[[1]]),
            cluster_v2 = as.factor(dataset$true_labels[[2]])
        )

        plots[[i]] <- ggplot(df, aes(x = x, y = y, color = cluster_v1, shape = cluster_v2)) +
            geom_point(size = 2.5, alpha = 0.7) +
            labs(
                title = sprintf("Sep = %.1f sigma", sep),
                x = "View 1", y = "View 2",
                color = "V1 Cluster", shape = "V2 Cluster"
            ) +
            theme_minimal(base_size = 10) +
            theme(legend.position = "none")
    }

    plot_grid(plotlist = plots, ncol = 2, labels = "AUTO")
}

cat("\n\n")
cat("============================================================\n")
cat("   MULTIVIEW GIBBS SAMPLER - CLUSTER SEPARATION STUDY\n")
cat("============================================================\n\n")

separations <- seq(6, 0.5, by = -0.5)

results <- run_separation_sweep(
    separations = separations,
    n_obs = 200,
    nsim = 5000,
    burn_in = 4000,
    thin = 1
)

cat("\nGenerating visualizations...\n")

p_ari_curve <- plot_ari_curve(results)
p_ari_by_view <- plot_ari_by_view(results)
p_scatter <- plot_scatter_panels(separations = c(6, 3, 2, 1))

print(p_ari_curve)
print(p_ari_by_view)
print(p_scatter)

ggsave("separation_ari_curve.png", p_ari_curve, width = 10, height = 6, dpi = 150)
ggsave("separation_ari_by_view.png", p_ari_by_view, width = 10, height = 6, dpi = 150)
ggsave("cluster_scatter_panels.png", p_scatter, width = 10, height = 8, dpi = 150)

cat("\n\n=======================================================\n")
cat("SUMMARY STATISTICS\n")
cat("=======================================================\n")

summary_stats <- results %>%
    group_by(separation) %>%
    summarize(
        mean_ari = round(mean(ari), 3),
        min_ari = round(min(ari), 3),
        correct_k = sum(cluster_match),
        .groups = "drop"
    )

print(summary_stats)

critical_sep <- summary_stats %>%
    filter(mean_ari < 0.9) %>%
    slice_max(separation, n = 1)

if (nrow(critical_sep) > 0) {
    cat(sprintf(
        "\n*** CONCLUSION: Model starts struggling (ARI < 0.9) at separation = %.1f sigma ***\n",
        critical_sep$separation
    ))
} else {
    cat("\n*** CONCLUSION: Model maintains ARI >= 0.9 across all tested separations ***\n")
}

cat("\nPlots saved to:\n")
cat("  - separation_ari_curve.png\n")
cat("  - separation_ari_by_view.png\n")
cat("  - cluster_scatter_panels.png\n")
