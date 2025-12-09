library(dplyr)
library(ggplot2)
library(cowplot)
library(Rcpp)
library(mcclust)
library(mclust)
library(mcclust.ext)
library(gridExtra)
library(clue)   
Rcpp::sourceCpp("multiview_gibbs.cpp")

set.seed(1999)


x <- data.frame(
  view1 = c(rnorm(100, mean = 10, sd = 1), 
            rnorm(100, mean = -10, sd = 1)), 
  view2 = c(rnorm(50, mean = 0, sd = 1), 
            rnorm(100, mean = -10, sd = 1), 
            rnorm(50, mean = 10, sd = 1)), 
  view3 = c(rnorm(100, mean = 10, sd = 1), 
            rnorm(100, mean = -10, sd = 1)),
  view4 = c(rnorm(100, mean = 10, sd = 1), 
            rnorm(100, mean = -10, sd = 1)),
  view5 = c(rnorm(100, mean = -5, sd = 1), 
            rnorm(100, mean = 5, sd = 1)),
  tables = sample(c(1:19), size = 200, replace = TRUE)
)


x <- data.frame(
  view1 = c(rnorm(100, mean = 3, sd = 1), 
            rnorm(100, mean = -3, sd = 1)), 
  view2 = c(rnorm(50, mean = 0, sd = 1), 
            rnorm(100, mean = -5, sd = 1), 
            rnorm(50, mean = 5, sd = 1)), 
  view3 = c(rnorm(100, mean = 3, sd = 1), 
            rnorm(100, mean = -3, sd = 1)),
  view4 = c(rnorm(100, mean = 3, sd = 1), 
            rnorm(100, mean = -3, sd = 1)),
  view5 = c(rnorm(100, mean = -3, sd = 1), 
            rnorm(100, mean = 3, sd = 1)),
  tables = sample(c(1:19), size = 200, replace = TRUE)
)


x <- data.frame(
  view1 = c(rnorm(100, mean = 3, sd = 1.3), 
            rnorm(100, mean = -3, sd = 1.3)), 
  view2 = c(rnorm(50, mean = 0, sd = 1.3), 
            rnorm(100, mean = -5, sd = 1.3), 
            rnorm(50, mean = 5, sd = 1.3)), 
  view3 = c(rnorm(100, mean = 3, sd = 1.3), 
            rnorm(100, mean = -3, sd = 1.3)),
  view4 = c(rnorm(100, mean = 3, sd = 1.3), 
            rnorm(100, mean = -3, sd = 1.3)),
  view5 = c(rnorm(100, mean = -3, sd = 1.3), 
            rnorm(100, mean = 3, sd = 1.3)),
  tables = sample(c(1:19), size = 200, replace = TRUE)
)


true_labels_2_clusters <- c(rep(1, 100), rep(2, 100))
true_labels_3_clusters <- c(rep(1, 50), rep(2, 100), rep(3, 50))

print(paste("Total rows:", nrow(x)))
print(head(x))

get_labels_for_view <- function(view_name) {
  if (view_name == "view2") {
    return(as.factor(true_labels_3_clusters))
  } else {
    return(as.factor(true_labels_2_clusters))
  }
}

plot_shape_color <- function(data, view_x_name, view_y_name) {
  x_data <- data[[view_x_name]]
  y_data <- data[[view_y_name]]
  labels_for_shape <- get_labels_for_view(view_x_name)
  labels_for_color <- get_labels_for_view(view_y_name)
  df_plot <- data.frame(
    x_val = x_data,
    y_val = y_data,
    Shape_Cluster = labels_for_shape, 
    Color_Cluster = labels_for_color
  )
  p <- ggplot(df_plot, aes(x = x_val, y = y_val, 
                           shape = Shape_Cluster, 
                           color = Color_Cluster)) +
    geom_point(size = 3.5, alpha = 0.8) +
    labs(
      title = paste("Interaction:", view_x_name, "vs", view_y_name),
      subtitle = paste("Shape by", view_x_name, "| Color by", view_y_name),
      x = view_x_name,
      y = view_y_name,
      shape = paste("Clusters", view_x_name),
      color = paste("Clusters", view_y_name)
    ) +
    theme_bw()
  return(p)
}

p1 <- plot_shape_color(x, "view2", "view3")
print(p1)

data_views <- list(
  as.vector(x$view1),
  as.vector(x$view2),
  as.vector(x$view3),
  as.vector(x$view4),
  as.vector(x$view5)
)

true_clust_list <- list(
  true_labels_2_clusters,
  true_labels_3_clusters,
  true_labels_2_clusters,
  true_labels_2_clusters,
  true_labels_2_clusters
)

V <- 5
nsim    <- 10000
burn_in <- 9000
thin    <- 1

res_gibbs <- run_gibbs_cpp(
  data_views = data_views,
  M          = nsim,
  burn_in    = burn_in,
  thin       = thin
)

get_final_clusters <- function(res_gibbs) {
  last_iter_idx <- length(res_gibbs$table_of)
  raw_tables <- res_gibbs$table_of[[last_iter_idx]]
  raw_dishes <- res_gibbs$dish_of[[last_iter_idx]]
  tables_r_index <- raw_tables + 1
  n_customers <- length(tables_r_index)
  n_views     <- length(raw_dishes)
  cluster_matrix <- matrix(NA, nrow = n_customers, ncol = n_views)
  colnames(cluster_matrix) <- paste0("View_", 1:n_views)
  for (v in 1:n_views) {
    dishes_for_view <- raw_dishes[[v]]
    cluster_matrix[, v] <- dishes_for_view[tables_r_index]
  }
  return(cluster_matrix)
}

my_clusters <- get_final_clusters(res_gibbs)
print("First 6 customers and their assigned clusters per view:")
print(head(my_clusters))
print(paste("Dimensions:", paste(dim(my_clusters), collapse = " x ")))

plot_predicted_interaction <- function(data, cluster_matrix, view_x_name, view_y_name) {
  idx_x <- as.numeric(gsub("view", "", view_x_name))
  idx_y <- as.numeric(gsub("view", "", view_y_name))
  x_data <- data[[view_x_name]]
  y_data <- data[[view_y_name]]
  pred_labels_x <- as.factor(cluster_matrix[, idx_x])
  pred_labels_y <- as.factor(cluster_matrix[, idx_y])
  df_plot <- data.frame(
    x_val = x_data,
    y_val = y_data,
    Shape_Cluster = pred_labels_x,
    Color_Cluster = pred_labels_y
  )
  p <- ggplot(df_plot, aes(x = x_val, y = y_val, 
                           shape = Shape_Cluster, 
                           color = Color_Cluster)) +
    geom_point(size = 3.5, alpha = 0.8) +
    labs(
      title = paste("Predicted: ", view_x_name, "vs", view_y_name),
      subtitle = paste("Shape:", view_x_name, "| Color:", view_y_name),
      x = view_x_name,
      y = view_y_name,
      shape = "Pred. Cluster (X)",
      color = "Pred. Cluster (Y)"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
  return(p)
}

p1 <- plot_predicted_interaction(x, my_clusters, "view2", "view3")
print(p1)

ari_scores <- sapply(1:5, function(v) mcclust::arandi(my_clusters[,v], true_clust_list[[v]]))
names(ari_scores) <- paste0("View_", 1:5)
print(ari_scores)

print(table(Predicted = my_clusters[,1], Truth = true_clust_list[[1]]))
print(table(Predicted = my_clusters[,2], Truth = true_clust_list[[2]]))
print(table(Predicted = my_clusters[,3], Truth = true_clust_list[[3]]))
print(table(Predicted = my_clusters[,4], Truth = true_clust_list[[4]]))
print(table(Predicted = my_clusters[,5], Truth = true_clust_list[[5]]))

df_global <- data.frame(
  iter = seq_along(res_gibbs$alpha_global),
  alpha_global = res_gibbs$alpha_global,
  sigma_global = res_gibbs$sigma_global
)

p_alpha_g <- ggplot(df_global, aes(iter, alpha_global)) +
  geom_line() + theme_minimal() +
  labs(title = "alpha_global trace")

p_sigma_g <- ggplot(df_global, aes(iter, sigma_global)) +
  geom_line() + theme_minimal() +
  labs(title = "sigma_global trace")

if (is.list(res_gibbs$alpha_v)) {
  n_views <- length(res_gibbs$alpha_v)
  n_iter  <- length(res_gibbs$alpha_v[[1]])
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
  n_iter  <- nrow(res_gibbs$alpha_v)
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
  geom_line() + theme_minimal() +
  labs(title = "alpha_v by view", colour = "view")

p_sigma_v <- ggplot(alpha_df, aes(iter, sigma, colour = factor(view))) +
  geom_line() + theme_minimal() +
  labs(title = "sigma_v by view", colour = "view")

p_tau_v <- ggplot(alpha_df, aes(iter, tau, colour = factor(view))) +
  geom_line() + theme_minimal() +
  labs(title = "tau_v by view", colour = "view")

param_df <- tidyr::pivot_longer(
  alpha_df,
  cols      = c(alpha, sigma, tau),
  names_to  = "param",
  values_to = "value"
)

p_by_view <- ggplot(param_df, aes(iter, value, colour = param)) +
  geom_line() + theme_minimal() +
  facet_wrap(~ view, scales = "free_y") +
  labs(title = "View-specific hyperparameters", colour = "param")

print(p_alpha_g); print(p_sigma_g)
print(p_alpha_v); print(p_sigma_v); print(p_tau_v)
print(p_by_view)


window <- 500  

final_alpha_global <- tail(res_gibbs$alpha_global, 1)
final_sigma_global <- tail(res_gibbs$sigma_global, 1)

mean_alpha_global <- mean(tail(res_gibbs$alpha_global, window))
sd_alpha_global   <- sd(tail(res_gibbs$alpha_global, window))

mean_sigma_global <- mean(tail(res_gibbs$sigma_global, window))
sd_sigma_global   <- sd(tail(res_gibbs$sigma_global, window))

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
  cat("Final alpha_v:", tail(a_vec, 1), 
      " | Mean last", window, "=", mean(tail(a_vec, window)), "\n")
  cat("Final sigma_v:", tail(s_vec, 1), 
      " | Mean last", window, "=", mean(tail(s_vec, window)), "\n")
  cat("Final tau_v:", tail(t_vec, 1), 
      " | Mean last", window, "=", mean(tail(t_vec, window)), "\n")
}

