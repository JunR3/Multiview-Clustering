############################################################
## 0. Packages and C++ code
############################################################

library(Rcpp)
library(ggplot2)
library(mcclust)
library(mcclust.ext)

## Set working directory where the C++ files live
setwd("/Users/jr/Desktop/Code_BP/Multiview")   # <-- change this path

## Compile and load the C++ code
## multiview_state.cpp includes multiview_state.h
## multiview_gibbs.cpp includes multiview_gibbs.h e usa le funzioni di state
sourceCpp("multiview_state.cpp")
sourceCpp("multiview_gibbs.cpp")

## C++ side:
## run_gibbs_cpp(data_views, M, burn_in, thin) must return a list, e.g.:
## list(
##   table_of_samples = saved_table_of,          # list length S, each length n
##   dish_of_samples  = saved_dish_of,           # list length S, each [[v]][t]
##   loglik_trace     = saved_loglik             # numeric length S
## )


############################################################
## 1. Simulate multiview data (generic d views)
############################################################

set.seed(1999)

n <- 200      # number of customers
d <- 5        # number of views  <-- change only this

## Simple example: each view is a Gaussian with different mean
means <- seq(-2, 2, length.out = d)
X <- sapply(
  means,
  function(m) c(rnorm(n/2, m - 2, 1),
                rnorm(n/2, m + 2, 1))
)
X <- as.matrix(X)   # n x d
colnames(X) <- paste0("view", 1:d)
View(X)
## Build list of views: data_views[[v]] numeric vector length n
data_views <- lapply(1:d, function(j) as.numeric(X[, j]))
names(data_views) <- paste0("view", 1:d)


############################################################
## 2. MCMC settings
############################################################

M       <- 2000   # total Gibbs iterations
burn_in <- 1000   # burn-in
thin    <- 1      # thinning factor


############################################################
## 3. Run multiview Gibbs sampler in C++
############################################################

res_gibbs <- run_gibbs_cpp(
  data_views = data_views,
  M          = M,
  burn_in    = burn_in,
  thin       = thin
)

str(res_gibbs)  # per controllare i nomi effettivi


############################################################
## 4. Build per-view cluster label matrix (for mcclust)
############################################################

## Adjust these names if your C++ uses slightly different ones
table_samples <- res_gibbs$table_of_samples   # list length S, each length n
dish_samples  <- res_gibbs$dish_of_samples    # list length S, each is list per view

S <- length(table_samples)          # number of saved iterations
n <- length(table_samples[[1]])     # number of customers
d <- length(dish_samples[[1]])      # number of views (from C++ perspective)

## Helper: build S x n matrix of cluster labels for a given view v
build_clusters_for_view <- function(v, table_samples, dish_samples) {
  S <- length(table_samples)
  n <- length(table_samples[[1]])
  
  cl_mat <- matrix(NA_integer_, nrow = S, ncol = n)
  
  for (s in 1:S) {
    table_s <- table_samples[[s]]   # length n, table index for each i
    dish_s  <- dish_samples[[s]]    # list length d, each integer vector of length #tables
    ## We assume:
    ## - table_s[i] is 0-based or 1-based according to C++; adjust shift
    ## - dish_s[[v]][t_index] is cluster label for table t_index in view v
    
    ## Convert to 1-based indices for R if C++ stored 0-based:
    ## If you KNOW C++ used 0..T-1, uncomment the +1 line:
    # table_s_R <- table_s + 1L
    # cl_v <- integer(n)
    # for (i in 1:n) {
    #   t  <- table_s_R[i]
    #   cl_v[i] <- dish_s[[v]][t]
    # }
    
    ## If C++ already stores 1..T in table_of and dish_of, use directly:
    cl_v <- integer(n)
    for (i in 1:n) {
      t  <- table_s[i]               # table label in 1..T
      cl_v[i] <- dish_s[[v]][t]      # dish index for that table in view v
    }
    
    ## mcclust prefers labels >= 1, so if there are 0s add +1
    if (min(cl_v) == 0L) {
      cl_v <- cl_v + 1L
    }
    
    cl_mat[s, ] <- cl_v
  }
  
  cl_mat
}

## Example: choose a view to analyse
v_analyze <- 1
clusters_v <- build_clusters_for_view(v_analyze, table_samples, dish_samples)
dim(clusters_v)  # S x n


############################################################
## 5. Posterior Similarity Matrix (PSM) and representative partition
############################################################

psm_v    <- comp.psm(clusters_v)
cl_v_def <- minVI(psm_v)$cl   # length n, final clusters for view v_analyze


############################################################
## 6. Visualisation options
############################################################

## 6.1 Scatter plot for a chosen pair of views (vx, vy)
vx <- 1
vy <- 2

df_pair <- data.frame(
  x  = data_views[[vx]],
  y  = data_views[[vy]],
  cl = factor(cl_v_def)
)

p_pair <- ggplot(df_pair, aes(x = x, y = y, colour = cl)) +
  geom_point(size = 2) +
  labs(
    x = paste0("View ", vx),
    y = paste0("View ", vy),
    colour = paste0("Clusters (view ", v_analyze, ")")
  ) +
  theme_minimal()

print(p_pair)


## 6.2 Generic helper: plot clusters for any pair of views (d views)
plot_cluster_pair <- function(vx, vy, cl_view, data_views) {
  df <- data.frame(
    x  = data_views[[vx]],
    y  = data_views[[vy]],
    cl = factor(cl_view)
  )
  ggplot(df, aes(x = x, y = y, colour = cl)) +
    geom_point(size = 2) +
    labs(
      x = paste0("View ", vx),
      y = paste0("View ", vy),
      colour = "Cluster"
    ) +
    theme_minimal()
}

## Example: view 2 vs view 5, colored by clusters from view 1
if (d >= 5) {
  p_2_5 <- plot_cluster_pair(2, 5, cl_v_def, data_views)
  print(p_2_5)
}


## 6.3 (Optional) Star / radar plot for a single observation across all d views

Ymat <- do.call(cbind, data_views)         # n x d
colnames(Ymat) <- paste0("view", 1:d)

## Normalise each view to [0,1]
Y_scaled <- apply(Ymat, 2, function(z) (z - min(z)) / (max(z) - min(z)))

## Radar plot for observation i0
i0 <- 1
radar_df <- rbind(
  rep(1, d),
  rep(0, d),
  Y_scaled[i0, ]
)
rownames(radar_df) <- c("max", "min", paste0("obs", i0))

## If you want the radar:
## library(fmsb)
## radarchart(radar_df)