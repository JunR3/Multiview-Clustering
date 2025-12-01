library(dplyr)
library(ggplot2)
library(cowplot)
library(Rcpp)
library(mcclust)
library(mcclust.ext)
library(gridExtra)
library(clue)   
Rcpp::sourceCpp("multiview_gibbs.cpp")


simulate_multiview_data <- function(n = 200, V) {
  
  set.seed(1)
  
  # sizes of the three latent true clusters
  n1 <- round(0.40*n)
  n2 <- round(0.30*n)
  n3 <- n - n1 - n2
  
  # global latent cluster assignment (3 clusters)
  z  <- sample(c(rep(1,n1), rep(2,n2), rep(3,n3)))
  
  X_mat <- matrix(NA, n, V)
  true_clust_list <- vector("list", V)
  
  # ---- View 1: 2 clusters very well separated (+12 vs -12)
  true_clust_list[[1]] <- ifelse(z==1,1L,2L)
  X_mat[true_clust_list[[1]]==1,1] <- rnorm(sum(true_clust_list[[1]]==1), +12, 1)
  X_mat[true_clust_list[[1]]==2,1] <- rnorm(sum(true_clust_list[[1]]==2), -12, 1)
  
  # ---- View 2: 3 clusters very well separated (+12, 0, -12)
  true_clust_list[[2]] <- z
  X_mat[z==1,2] <- rnorm(sum(z==1), +12,1)
  X_mat[z==2,2] <- rnorm(sum(z==2),   0,1)
  X_mat[z==3,2] <- rnorm(sum(z==3), -12,1)
  
  # ---- Views >= 3: increasing separation (+15 vs -15), (+20 vs -20), ...
  if (V>=3) {
    for (v in 3:V) {
      mu <- 10 + 5*(v - 2)   # separations: 15, 20, 25, ...
      
      # Clustering for views >=3: {1,2} vs {3}
      true_clust_list[[v]] <- ifelse(z %in% c(1,2), 1L, 2L)
      
      X_mat[true_clust_list[[v]]==1, v] <- rnorm(sum(true_clust_list[[v]]==1), +mu,1)
      X_mat[true_clust_list[[v]]==2, v] <- rnorm(sum(true_clust_list[[v]]==2), -mu,1)
    }
  }
  
  x <- as.data.frame(X_mat)
  colnames(x) <- paste0("view",1:V)
  
  # list of vectors, one per view
  data_views <- lapply(seq_len(V), function(j) x[[j]])
  
  list(
    x = x,
    true_clust_list = true_clust_list,
    data_views = data_views
  )
}


# ---- Simulate 5 views
V <- 5
sim <- simulate_multiview_data(n = 200, V = V)

x               <- sim$x
true_clust_list <- sim$true_clust_list
data_views      <- sim$data_views

# ---- Run Gibbs sampler
nsim    <- 10000
burn_in <- 5000
thin    <- 1

res_gibbs <- run_gibbs_cpp(
  data_views = data_views,
  M          = nsim,
  burn_in    = burn_in,
  thin       = thin
)

n <- nrow(x)
S <- length(res_gibbs$table_of)

# allocate list of cluster matrices (one per view)
cluster_list <- lapply(seq_len(V),
                       function(v) matrix(NA_integer_, n, S))

# relabel clusters so labels start at 1,2,3,...
relabel_to_1K <- function(z) {
  z <- as.integer(z)
  if (min(z)<1L) z <- z - min(z) + 1L
  as.integer(as.numeric(factor(z)))
}

# ---- Extract cluster labels from Gibbs samples
for (s in seq_len(S)) {
  tab_s  <- res_gibbs$table_of[[s]]   # table index for each observation
  dish_s <- res_gibbs$dish_of[[s]]    # list of dish labels per view
  
  lab_mat <- matrix(NA_integer_, n, V)
  
  for (i in seq_len(n)) {
    t_i <- tab_s[i] + 1L
    for (v in seq_len(V))
      lab_mat[i,v] <- dish_s[[v]][t_i] + 1L
  }
  
  for (v in seq_len(V))
    cluster_list[[v]][,s] <- relabel_to_1K(lab_mat[,v])
}

# ---- Compare any pair of views
view_idx1 <- 2
view_idx2 <- 3  # could choose 4 or 5 as well

cl1_t <- t(cluster_list[[view_idx1]])
cl2_t <- t(cluster_list[[view_idx2]])

# pairwise similarity matrices
psm1 <- comp.psm(cl1_t)
psm2 <- comp.psm(cl2_t)

# VI-optimal partition
cluster1_def_raw <- minVI(psm1)$cl
cluster2_def_raw <- minVI(psm2)$cl


# Align estimated labels with true labels using Hungarian algorithm
align_to_true <- function(z_hat, z_true) {
  z_hat  <- as.integer(as.numeric(factor(z_hat)))
  z_true <- as.integer(as.numeric(factor(z_true)))
  
  K_hat  <- max(z_hat)
  K_true <- max(z_true)
  
  K <- max(K_hat, K_true)
  
  M <- matrix(0L, nrow = K, ncol = K)
  
  # fill contingency matrix
  for (i in seq_along(z_hat)) {
    M[z_hat[i], z_true[i]] <- M[z_hat[i], z_true[i]] + 1L
  }
  
  # convert to cost matrix
  cost <- max(M) - M
  
  # Hungarian algorithm
  perm <- solve_LSAP(cost)
  
  # apply permutation
  z_new <- z_hat
  for (k in seq_len(K_hat)) {
    z_new[z_hat == k] <- perm[k]
  }
  z_new
}

true1 <- true_clust_list[[view_idx1]]
true2 <- true_clust_list[[view_idx2]]

cluster1_def <- align_to_true(cluster1_def_raw, true1)
cluster2_def <- align_to_true(cluster2_def_raw, true2)

# add results to dataframe
x$cluster_v1_def <- cluster1_def
x$cluster_v2_def <- cluster2_def
x$true_v1        <- true1
x$true_v2        <- true2

# misclassified points
mis1 <- which(x$true_v1 != x$cluster_v1_def)
mis2 <- which(x$true_v2 != x$cluster_v2_def)

cat("Table view", view_idx1,":\n")
print(table(cluster1_def, true1))

cat("\nTable view", view_idx2,":\n")
print(table(cluster2_def, true2))

# ---- Scatterplot comparing two views
g_scatter <- ggplot(x, aes_string(x=paste0("view",view_idx1),
                                  y=paste0("view",view_idx2))) +
  geom_point(aes(colour=factor(cluster_v1_def),
                 shape=factor(cluster_v2_def)),
             size=3) +
  geom_point(colour="black",size=1.2) +
  geom_point(data=x[mis1,], colour="red",shape=1,size=4,stroke=1) +
  geom_point(data=x[mis2,], colour="red",shape=1,size=4,stroke=1) +
  theme_minimal()

print(g_scatter)

# ---- Histograms for all views
plots <- lapply(seq_len(V), function(v){
  ggplot(data.frame(val=x[[paste0("view",v)]]), aes(val)) +
    geom_histogram(bins=40, fill="lightblue", colour="darkblue") +
    ggtitle(paste("View",v)) +
    theme_minimal()
})

grid.arrange(grobs=plots, ncol=2)

