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

# --- Data Generation Step ---

x <- data.frame(
  # View 1: 2 Clusters (100 data points each). Means at 10 and -10.
  view1 = c(rnorm(100, mean = 10, sd = 1), 
            rnorm(100, mean = -10, sd = 1)), 
  
  # View 2: 3 Clusters. Split is 50/100/50 to sum 200. Means at 0, -12, 12.
  view2 = c(rnorm(50, mean = 0, sd = 1), 
            rnorm(100, mean = -10, sd = 1), 
            rnorm(50, mean = 10, sd = 1)), 
  
  # View 3: 2 Clusters (100 data points each). Means at 10 and -10.
  view3 = c(rnorm(100, mean = 10, sd = 1), 
            rnorm(100, mean = -10, sd = 1)),
  
  # View 4: 2 Clusters. Means at 9 and -9.
  view4 = c(rnorm(100, mean = 10, sd = 1), 
   rnorm(100, mean = -10, sd = 1)),
  # 
  # # View 5: 2 Clusters. Means at 1 and 8.
  view5 = c(rnorm(100, mean = -5, sd = 1), 
   rnorm(100, mean = 5, sd = 1)),
  
  # Categorical variable for context (Tables)
  tables = sample(c(1:19), size = 200, replace = TRUE)
)

# --- Ground Truth Labels ---
# We store the "real" cluster assignments to check accuracy later.

# Truth for Views 1, 3, 4, 5 (2 clusters: first 100 are group 1, next 100 are group 2)
true_labels_2_clusters <- c(rep(1, 100), rep(2, 100))

# Truth for View 2 (3 clusters: 50 group 1, 100 group 2, 50 group 3)
true_labels_3_clusters <- c(rep(1, 50), rep(2, 100), rep(3, 50))
#true_labels_3_clusters <- c(rep(1, 100), rep(2, 100))

# --- Verification ---
print(paste("Total rows:", nrow(x)))
print(head(x))

# --- Helper Function to select correct Ground Truth ---
# Since we generated data manually:
# - "view2" uses 'true_clust_3groups'
# - All other views use 'true_clust_2groups'
get_labels_for_view <- function(view_name) {
  if (view_name == "view2") {
    return(as.factor(true_labels_3_clusters))
    #return(as.factor(true_labels_3_clusters))
  } else {
    return(as.factor(true_labels_2_clusters))
  }
}

# --- Main Plotting Function ---
plot_shape_color <- function(data, view_x_name, view_y_name) {
  
  # 1. Get the data for X and Y axes
  x_data <- data[[view_x_name]]
  y_data <- data[[view_y_name]]
  
  # 2. Get the corresponding True Labels for each view
  # Labels for X-axis view -> determine SHAPE
  labels_for_shape <- get_labels_for_view(view_x_name)
  
  # Labels for Y-axis view -> determine COLOR
  labels_for_color <- get_labels_for_view(view_y_name)
  
  # 3. Create a temporary dataframe for plotting
  df_plot <- data.frame(
    x_val = x_data,
    y_val = y_data,
    Shape_Cluster = labels_for_shape, 
    Color_Cluster = labels_for_color
  )
  
  # 4. Generate Plot
  p <- ggplot(df_plot, aes(x = x_val, y = y_val, 
                           shape = Shape_Cluster, 
                           color = Color_Cluster)) +
    geom_point(size = 3.5, alpha = 0.8) + # Large points to see shapes clearly
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
  as.vector(x$view1), # View 1
  as.vector(x$view2), # View 2
  as.vector(x$view3), # View 3
  as.vector(x$view4), # View 4
  as.vector(x$view5)  # View 5
)

# We map the ground truth vectors we created earlier to each view.
# Remember: View 2 has 3 clusters, the rest have 2.

true_clust_list <- list(
  true_labels_2_clusters, # Truth for View 1
  true_labels_3_clusters, # Truth for View 2 (3 clusters)
  true_labels_2_clusters, # Truth for View 3
  true_labels_2_clusters, # Truth for View 4
  true_labels_2_clusters  # Truth for View 5
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
  
  # 1. Identify the index of the last iteration
  # We assume the last iteration represents the converged state
  last_iter_idx <- length(res_gibbs$table_of)
  
  # 2. Extract raw assignments from the Gibbs output
  # 'table_of': Vector of length N (customers). Maps Customer -> Table ID.
  raw_tables <- res_gibbs$table_of[[last_iter_idx]]
  
  # 'dish_of': List of length V (views). Each element is a vector of length T (tables).
  raw_dishes <- res_gibbs$dish_of[[last_iter_idx]]
  
  # 3. Adjust indices from 0-based (C++) to 1-based (R)
  # We need 1-based indices to use them as lookup positions in R vectors.
  # We assume raw_tables contains values like 0, 1, 2...
  tables_r_index <- raw_tables + 1
  
  # 4. Initialize the result matrix
  n_customers <- length(tables_r_index)
  n_views     <- length(raw_dishes)
  
  # Matrix: Rows = Customers, Columns = Views
  cluster_matrix <- matrix(NA, nrow = n_customers, ncol = n_views)
  colnames(cluster_matrix) <- paste0("View_", 1:n_views)
  
  # 5. Perform the mapping: Customer -> Table -> Dish
  for (v in 1:n_views) {
    # Get the specific dish configuration for view 'v'
    # logical structure: dishes_for_view[table_id] = dish_id
    dishes_for_view <- raw_dishes[[v]]
    
    # Map each customer to their cluster in this view using their table assignment
    cluster_matrix[, v] <- dishes_for_view[tables_r_index]
  }
  
  return(cluster_matrix)
}

# --- USAGE EXAMPLE ---

# Extract clusters
my_clusters <- get_final_clusters(res_gibbs)

# Inspect the first few rows
print("First 6 customers and their assigned clusters per view:")
print(head(my_clusters))

# Check dimensions 
print(paste("Dimensions:", paste(dim(my_clusters), collapse = " x ")))

# --- Function to Plot Predicted Clusters ---
plot_predicted_interaction <- function(data, cluster_matrix, view_x_name, view_y_name) {
  
  # 1. Parse view indices from names (e.g., "view1" -> 1, "view2" -> 2)
  idx_x <- as.numeric(gsub("view", "", view_x_name))
  idx_y <- as.numeric(gsub("view", "", view_y_name))
  
  # 2. Extract data coordinates
  x_data <- data[[view_x_name]]
  y_data <- data[[view_y_name]]
  
  # 3. Extract PREDICTED labels from your results matrix
  # We convert to factor to ensure discrete colors/shapes
  pred_labels_x <- as.factor(cluster_matrix[, idx_x])
  pred_labels_y <- as.factor(cluster_matrix[, idx_y])
  
  # 4. Create temporary dataframe for ggplot
  df_plot <- data.frame(
    x_val = x_data,
    y_val = y_data,
    Shape_Cluster = pred_labels_x,
    Color_Cluster = pred_labels_y
  )
  
  # 5. Generate Plot
  p <- ggplot(df_plot, aes(x = x_val, y = y_val, 
                           shape = Shape_Cluster, 
                           color = Color_Cluster)) +
    geom_point(size = 3.5, alpha = 0.8) + 
    labs(
      title = paste("Predicted: ", view_x_name, "vs", view_y_name),
      subtitle = paste("Shape: Clusters found in", view_x_name, "| Color: Clusters found in", view_y_name),
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

print(table(Predicted = my_clusters[,2], Truth = true_clust_list[[2]]))

################################################################################
#
#                      TRACE PLOTS OF HYPERPARAMETERS
#
################################################################################



# Global traces
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

# Per-view: alpha_v/sigma_v/tau_v are list(view) -> numeric(iter)
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

# One plot per parameter, lines colored by view
p_alpha_v <- ggplot(alpha_df, aes(iter, alpha, colour = factor(view))) +
  geom_line() + theme_minimal() +
  labs(title = "alpha_v by view", colour = "view")
p_sigma_v <- ggplot(alpha_df, aes(iter, sigma, colour = factor(view))) +
  geom_line() + theme_minimal() +
  labs(title = "sigma_v by view", colour = "view")
p_tau_v <- ggplot(alpha_df, aes(iter, tau, colour = factor(view))) +
  geom_line() + theme_minimal() +
  labs(title = "tau_v by view", colour = "view")

# All params per view in one panel
param_df <- pivot_longer(alpha_df, cols = c(alpha, sigma, tau),
                         names_to = "param", values_to = "value")
p_by_view <- ggplot(param_df, aes(iter, value, colour = param)) +
  geom_line() + theme_minimal() +
  facet_wrap(~ view, scales = "free_y") +
  labs(title = "View-specific hyperparameters", colour = "param")

# Print our graphs
print(p_alpha_g); print(p_sigma_g)
print(p_alpha_v); print(p_sigma_v); print(p_tau_v)
# or
print(p_by_view)
