library(ggplot2)
library(cowplot)
library(Rcpp)
library(dplyr)
library(tidyr)
library(mcclust)

Rcpp::sourceCpp("lib/multiview_gibbs.cpp")
set.seed(1999)
run_one_hospital <- function(csv_path = NULL,
                             df_all = NULL,
                             hospital_id = 1,
                             view_cols = c("gest", "dde", "weight"),
                             nsim = 1000,
                             burn_in = 900,
                             thin = 1,
                             bins = 30) {
  if (is.null(df_all)) {
    df_all <- read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)
    bad_names <- which(is.na(names(df_all)) | names(df_all) == "")
    if (length(bad_names) > 0) names(df_all)[bad_names] <- paste0("X", bad_names)
  }

  stopifnot("hosp" %in% names(df_all))
  stopifnot(all(view_cols %in% names(df_all)))

  df_hosp <- df_all[df_all$hosp == hospital_id, , drop = FALSE]
  n <- nrow(df_hosp)
  stopifnot(n > 1)

  views_df <- df_hosp[, view_cols, drop = FALSE]
  views_df <- as.data.frame(lapply(views_df, function(x) as.numeric(as.character(x))))

  for (j in seq_len(ncol(views_df))) {
    if (anyNA(views_df[[j]])) {
      m <- mean(views_df[[j]], na.rm = TRUE)
      views_df[[j]][is.na(views_df[[j]])] <- m
    }
  }

  stats <- data.frame(
    view = view_cols,
    mean = sapply(views_df, mean, na.rm = TRUE),
    sd   = sapply(views_df, sd, na.rm = TRUE),
    min  = sapply(views_df, min, na.rm = TRUE),
    max  = sapply(views_df, max, na.rm = TRUE)
  )

  df_long <- tidyr::pivot_longer(
    dplyr::mutate(views_df, .id = seq_len(n)),
    cols = all_of(view_cols),
    names_to = "view",
    values_to = "value"
  )

  p_dist <- ggplot(df_long, aes(x = value)) +
    geom_histogram(aes(y = after_stat(density)), bins = bins, alpha = 0.5) +
    geom_density(linewidth = 1.2) +
    facet_wrap(~view, scales = "free_x", nrow = 1) +
    theme_bw() +
    labs(
      title = paste0("Hospital ", hospital_id, " — raw distributions"),
      x = "value",
      y = "density"
    )

  sds <- stats$sd
  sds[is.na(sds) | sds == 0] <- 1
  views_df_scaled <- views_df
  for (j in seq_along(view_cols)) views_df_scaled[[j]] <- views_df_scaled[[j]] / sds[j]

  data_views <- lapply(view_cols, function(v) as.vector(views_df_scaled[[v]]))

  res_gibbs <- run_gibbs_cpp(
    data_views = data_views,
    M = nsim,
    burn_in = burn_in,
    thin = thin
  )

  get_final_clusters <- function(res) {
    last_iter_idx <- length(res$table_of)
    raw_tables <- res$table_of[[last_iter_idx]]
    raw_dishes <- res$dish_of[[last_iter_idx]]
    tables_r_index <- raw_tables + 1
    n_customers <- length(tables_r_index)
    n_views <- length(raw_dishes)
    cluster_matrix <- matrix(NA_integer_, nrow = n_customers, ncol = n_views)
    colnames(cluster_matrix) <- paste0("View_", seq_len(n_views))
    for (v in seq_len(n_views)) {
      dishes_for_view <- raw_dishes[[v]]
      cluster_matrix[, v] <- dishes_for_view[tables_r_index]
    }
    cluster_matrix
  }

  clusters <- get_final_clusters(res_gibbs)
  k_per_view <- sapply(seq_len(ncol(clusters)), function(j) length(unique(clusters[, j])))
  names(k_per_view) <- paste0("View_", seq_len(ncol(clusters)))

  df_global <- data.frame(
    iter = seq_along(res_gibbs$alpha_global),
    alpha_global = res_gibbs$alpha_global,
    sigma_global = res_gibbs$sigma_global
  )

  p_alpha_g <- ggplot(df_global, aes(iter, alpha_global)) +
    geom_line(alpha = 0.8) +
    theme_minimal() +
    labs(title = "alpha_global trace", x = "iteration", y = "alpha_global")

  p_sigma_g <- ggplot(df_global, aes(iter, sigma_global)) +
    geom_line(alpha = 0.8) +
    theme_minimal() +
    labs(title = "sigma_global trace", x = "iteration", y = "sigma_global")

  if (is.list(res_gibbs$alpha_v)) {
    n_views <- length(res_gibbs$alpha_v)
    n_iter <- length(res_gibbs$alpha_v[[1]])
    alpha_df <- dplyr::bind_rows(lapply(seq_len(n_views), function(v) {
      data.frame(
        iter = seq_len(n_iter), view = v,
        alpha = res_gibbs$alpha_v[[v]],
        sigma = res_gibbs$sigma_v[[v]],
        tau = res_gibbs$tau_v[[v]]
      )
    }))
  } else if (is.matrix(res_gibbs$alpha_v)) {
    n_iter <- nrow(res_gibbs$alpha_v)
    n_views <- ncol(res_gibbs$alpha_v)
    alpha_df <- dplyr::bind_rows(lapply(seq_len(n_views), function(v) {
      data.frame(
        iter = seq_len(n_iter), view = v,
        alpha = res_gibbs$alpha_v[, v],
        sigma = res_gibbs$sigma_v[, v],
        tau = res_gibbs$tau_v[, v]
      )
    }))
  } else {
    stop("alpha_v has an unsupported type.")
  }

  p_alpha_v <- ggplot(alpha_df, aes(iter, alpha, colour = factor(view))) +
    geom_line(alpha = 0.7) +
    theme_minimal() +
    labs(title = "alpha_v by view", x = "iteration", y = "alpha_v", colour = "view")

  p_sigma_v <- ggplot(alpha_df, aes(iter, sigma, colour = factor(view))) +
    geom_line(alpha = 0.7) +
    theme_minimal() +
    labs(title = "sigma_v by view", x = "iteration", y = "sigma_v", colour = "view")

  p_tau_v <- ggplot(alpha_df, aes(iter, tau, colour = factor(view))) +
    geom_line(alpha = 0.7) +
    theme_minimal() +
    labs(title = "tau_v by view", x = "iteration", y = "tau_v", colour = "view")

  list(
    hospital_id = hospital_id,
    n = n,
    view_cols = view_cols,
    data = views_df,
    stats = stats,
    res_gibbs = res_gibbs,
    clusters = clusters,
    k_per_view = k_per_view,
    plots = list(
      dist = p_dist,
      alpha_global = p_alpha_g,
      sigma_global = p_sigma_g,
      alpha_v = p_alpha_v,
      sigma_v = p_sigma_v,
      tau_v = p_tau_v
    )
  )
}

run_many_hospitals <- function(csv_path,
                               hospital_ids,
                               view_cols = c("gest", "dde", "weight"),
                               nsim = 1000,
                               burn_in = 900,
                               thin = 1,
                               bins = 30) {
  df_all <- read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)
  bad_names <- which(is.na(names(df_all)) | names(df_all) == "")
  if (length(bad_names) > 0) names(df_all)[bad_names] <- paste0("X", bad_names)

  outs <- lapply(hospital_ids, function(hid) {
    run_one_hospital(
      df_all = df_all,
      hospital_id = hid,
      view_cols = view_cols,
      nsim = nsim,
      burn_in = burn_in,
      thin = thin,
      bins = bins
    )
  })
  names(outs) <- paste0("hosp_", hospital_ids)
  outs
}


outs <- run_many_hospitals(
  csv_path = "CPP_dataset.csv",
  hospital_ids = c(1, 3, 6, 10),
  view_cols = c("gest", "dde", "weight"),
  nsim = 5000,
  burn_in = 4000
)

sapply(outs, function(o) o$k_per_view)
