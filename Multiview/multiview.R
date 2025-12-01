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
  
  n1 <- round(0.40*n)
  n2 <- round(0.30*n)
  n3 <- n - n1 - n2
  z  <- sample(c(rep(1,n1), rep(2,n2), rep(3,n3)))
  
  X_mat <- matrix(NA, n, V)
  true_clust_list <- vector("list", V)
  
  true_clust_list[[1]] <- ifelse(z==1,1L,2L)
  X_mat[true_clust_list[[1]]==1,1] <- rnorm(sum(true_clust_list[[1]]==1), +5, 1)
  X_mat[true_clust_list[[1]]==2,1] <- rnorm(sum(true_clust_list[[1]]==2), -5, 1)
  
  true_clust_list[[2]] <- z
  X_mat[z==1,2] <- rnorm(sum(z==1), +6,1)
  X_mat[z==2,2] <- rnorm(sum(z==2),  0,1)
  X_mat[z==3,2] <- rnorm(sum(z==3), -6,1)
  
  if (V>=3) {
    for (v in 3:V) {
      mu <- 4 + 2*(v-1)
      true_clust_list[[v]] <- ifelse(z==1,1L,2L)
      X_mat[true_clust_list[[v]]==1, v] <- rnorm(sum(true_clust_list[[v]]==1), +mu,1)
      X_mat[true_clust_list[[v]]==2, v] <- rnorm(sum(true_clust_list[[v]]==2), -mu,1)
    }
  }
  
  x <- as.data.frame(X_mat)
  colnames(x) <- paste0("view",1:V)
  
  data_views <- lapply(seq_len(V), function(j) x[[j]])
  
  list(
    x = x,
    true_clust_list = true_clust_list,
    data_views = data_views
  )
}



V <- 5
sim <- simulate_multiview_data(n=200, V=V)

x               <- sim$x
true_clust_list <- sim$true_clust_list
data_views      <- sim$data_views

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

cluster_list <- lapply(seq_len(V),
                       function(v) matrix(NA_integer_, n, S))

relabel_to_1K <- function(z) {
  z <- as.integer(z)
  if (min(z)<1L) z <- z - min(z) + 1L
  as.integer(as.numeric(factor(z)))
}

for (s in seq_len(S)) {
  tab_s  <- res_gibbs$table_of[[s]]
  dish_s <- res_gibbs$dish_of[[s]]
  
  lab_mat <- matrix(NA_integer_, n, V)
  
  for (i in seq_len(n)) {
    t_i <- tab_s[i] + 1L
    for (v in seq_len(V))
      lab_mat[i,v] <- dish_s[[v]][t_i] + 1L
  }
  
  for (v in seq_len(V))
    cluster_list[[v]][,s] <- relabel_to_1K(lab_mat[,v])
}



view_idx1 <- 1
view_idx2 <- 3

cl1_t <- t(cluster_list[[view_idx1]])
cl2_t <- t(cluster_list[[view_idx2]])

psm1 <- comp.psm(cl1_t)
psm2 <- comp.psm(cl2_t)

cluster1_def_raw <- minVI(psm1)$cl
cluster2_def_raw <- minVI(psm2)$cl



align_to_true <- function(z_hat, z_true) {
  z_hat  <- as.integer(z_hat)
  z_true <- as.integer(z_true)
  
  z_hat  <- as.integer(as.numeric(factor(z_hat)))
  z_true <- as.integer(as.numeric(factor(z_true)))
  
  K_hat  <- max(z_hat)
  K_true <- max(z_true)
  

  M <- matrix(0L, nrow = K_hat, ncol = K_true)
  for (i in seq_along(z_hat)) {
    M[z_hat[i], z_true[i]] <- M[z_hat[i], z_true[i]] + 1L
  }
  
  cost <- max(M) - M
  perm <- solve_LSAP(cost)   
  
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


x$cluster_v1_def <- cluster1_def
x$cluster_v2_def <- cluster2_def
x$true_v1        <- true1
x$true_v2        <- true2

mis1 <- which(x$true_v1 != x$cluster_v1_def)
mis2 <- which(x$true_v2 != x$cluster_v2_def)

cat("Tabella view", view_idx1,":\n")
print(table(cluster1_def, true1))

cat("\nTabella view", view_idx2,":\n")
print(table(cluster2_def, true2))

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


plots <- lapply(seq_len(V), function(v){
  ggplot(data.frame(val=x[[paste0("view",v)]]), aes(val)) +
    geom_histogram(bins=40, fill="lightblue", colour="darkblue") +
    ggtitle(paste("View",v)) +
    theme_minimal()
})

grid.arrange(grobs=plots, ncol=2)

