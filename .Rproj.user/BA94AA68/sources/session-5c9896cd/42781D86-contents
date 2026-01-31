library(ggplot2)
library(cowplot)
library(Rcpp)
library(dplyr)
library(tidyr)

set.seed(1999)

load_cpp_data <- function(csv_path = NULL, df_all = NULL) {
  if (is.null(df_all)) {
    stopifnot(!is.null(csv_path))
    df_all <- read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)
  }
  bad_names <- which(is.na(names(df_all)) | names(df_all) == "")
  if (length(bad_names) > 0) names(df_all)[bad_names] <- paste0("X", bad_names)
  df_all
}

subset_cpp_data <- function(df_all,
                            hospital_ids = 1,
                            smoke_level = NULL) {
  stopifnot("hosp" %in% names(df_all))
  if (!is.null(smoke_level)) stopifnot("smoke" %in% names(df_all))
  
  df_sub <- df_all[df_all$hosp %in% hospital_ids, , drop = FALSE]
  if (!is.null(smoke_level)) df_sub <- df_sub[df_sub$smoke %in% smoke_level, , drop = FALSE]
  
  stopifnot(nrow(df_sub) > 1)
  df_sub
}

make_views <- function(df_sub,
                       view_cols = c("gest", "dde", "weight"),
                       log_dde = TRUE,
                       divide_by_sd = TRUE) {
  stopifnot(all(view_cols %in% names(df_sub)))
  
  views_raw <- df_sub[, view_cols, drop = FALSE]
  views_raw <- as.data.frame(lapply(views_raw, function(x) as.numeric(as.character(x))))
  
  for (j in seq_len(ncol(views_raw))) {
    if (anyNA(views_raw[[j]])) {
      m <- mean(views_raw[[j]], na.rm = TRUE)
      views_raw[[j]][is.na(views_raw[[j]])] <- m
    }
  }
  
  views <- views_raw
  view_names <- view_cols
  
  if (log_dde && "dde" %in% view_cols) {
    dde_idx <- match("dde", view_cols)
    x <- views[[dde_idx]]
    if (any(x <= 0, na.rm = TRUE)) {
      views[[dde_idx]] <- log1p(x)
      view_names[dde_idx] <- "log1p(dde)"
    } else {
      views[[dde_idx]] <- log(x)
      view_names[dde_idx] <- "log(dde)"
    }
  }
  
  stats <- data.frame(
    view = view_names,
    mean = sapply(views, mean, na.rm = TRUE),
    sd   = sapply(views, sd,   na.rm = TRUE),
    min  = sapply(views, min,  na.rm = TRUE),
    max  = sapply(views, max,  na.rm = TRUE)
  )
  
  views_scaled <- views
  if (divide_by_sd) {
    sds <- stats$sd
    sds[is.na(sds) | sds == 0] <- 1
    for (j in seq_len(ncol(views_scaled))) views_scaled[[j]] <- views_scaled[[j]] / sds[j]
  }
  
  list(
    views_raw = views,
    views_scaled = views_scaled,
    view_names = view_names,
    stats = stats
  )
}

plot_distributions <- function(views_df, view_names, bins = 30, title = "raw distributions") {
  df_long <- tidyr::pivot_longer(
    dplyr::mutate(views_df, .id = seq_len(nrow(views_df))),
    cols = everything(),
    names_to = "view",
    values_to = "value"
  )
  
  map_names <- setNames(view_names, names(views_df))
  df_long$view <- factor(df_long$view, levels = names(map_names), labels = unname(map_names))
  
  ggplot(df_long, aes(x = value)) +
    geom_histogram(aes(y = after_stat(density)), bins = bins, alpha = 0.5) +
    geom_density(linewidth = 1.2) +
    facet_wrap(~ view, scales = "free_x", nrow = 1) +
    theme_bw() +
    labs(title = title, x = "value", y = "density")
}

plot_scatter <- function(views_df, view_names, i = 1, j = 2, title = NULL) {
  stopifnot(i != j, i >= 1, j >= 1, i <= ncol(views_df), j <= ncol(views_df))
  
  dfp <- data.frame(
    x = views_df[[i]],
    y = views_df[[j]]
  )
  
  if (is.null(title)) title <- paste0(view_names[i], " vs ", view_names[j])
  
  ggplot(dfp, aes(x = x, y = y)) +
    geom_point(alpha = 0.6) +
    theme_bw() +
    labs(title = title, x = view_names[i], y = view_names[j])
}

compute_correlations <- function(views_df, view_names, method = "pearson") {
  M <- cor(views_df, use = "pairwise.complete.obs", method = method)
  colnames(M) <- view_names
  rownames(M) <- view_names
  M
}

step1_prepare <- function(csv_path = NULL,
                          df_all = NULL,
                          hospital_ids = 1,
                          smoke_level = NULL,
                          view_cols = c("gest", "dde", "weight"),
                          log_dde = TRUE,
                          divide_by_sd = TRUE,
                          bins = 30) {
  df_all <- load_cpp_data(csv_path = csv_path, df_all = df_all)
  df_sub <- subset_cpp_data(df_all, hospital_ids = hospital_ids, smoke_level = smoke_level)
  
  v <- make_views(
    df_sub = df_sub,
    view_cols = view_cols,
    log_dde = log_dde,
    divide_by_sd = divide_by_sd
  )
  
  title_raw <- paste0(
    "Hospitals: ", paste(hospital_ids, collapse = ", "),
    if (!is.null(smoke_level)) paste0(" | smoke: ", paste(smoke_level, collapse = ",")) else "",
    " — raw distributions"
  )
  title_scaled <- paste0(
    "Hospitals: ", paste(hospital_ids, collapse = ", "),
    if (!is.null(smoke_level)) paste0(" | smoke: ", paste(smoke_level, collapse = ",")) else "",
    " — scaled by sd distributions"
  )
  
  p_raw <- plot_distributions(v$views_raw, v$view_names, bins = bins, title = title_raw)
  p_scaled <- plot_distributions(v$views_scaled, v$view_names, bins = bins, title = title_scaled)
  
  cor_raw <- compute_correlations(v$views_raw, v$view_names, method = "pearson")
  cor_scaled <- compute_correlations(v$views_scaled, v$view_names, method = "pearson")
  
  list(
    df_all = df_all,
    df_sub = df_sub,
    hospital_ids = hospital_ids,
    smoke_level = smoke_level,
    view_cols = view_cols,
    view_names = v$view_names,
    views_raw = v$views_raw,
    views_scaled = v$views_scaled,
    stats = v$stats,
    cor_raw = cor_raw,
    cor_scaled = cor_scaled,
    plots = list(
      dist_raw = p_raw,
      dist_scaled = p_scaled
    )
  )
}

s1 <- step1_prepare(
  csv_path = "CPP_dataset.csv",
  hospital_ids = c(1, 3, 6, 10),
  smoke_level = NULL,
  view_cols = c("gest", "dde", "weight"),
  log_dde = TRUE,
  divide_by_sd = TRUE,
  bins = 30
)

print(s1$stats)
print(s1$cor_raw)
print(s1$cor_scaled)

print(s1$plots$dist_raw)
print(s1$plots$dist_scaled)

print(plot_scatter(s1$views_raw, s1$view_names, 1, 3))
print(plot_scatter(s1$views_raw, s1$view_names, 1, 2))
print(plot_scatter(s1$views_raw, s1$view_names, 2, 3))