library(Rcpp)
library(ggplot2)
library(mcclust)
library(mcclust.ext)

setwd("/Users/jr/Desktop/Code_BP/Multiview") 


## multiview_state.cpp includes multiview_state.h
## multiview_gibbs.cpp includes multiview_gibbs.h 

sourceCpp("multiview_state.cpp")
sourceCpp("multiview_gibbs.cpp")

############################################################
## 1. Simulate multiview data for generic d views
############################################################

set.seed(1999)

n <- 200      # number of customers
d <- 5        # number of views 

## each view is a Gaussian with different mean
means <- seq(-2, 2, length.out = d)
X <- sapply(
  means,
  function(m) c(rnorm(n/2, m - 2, 1),
                rnorm(n/2, m + 2, 1))
)
X <- as.matrix(X)   # n x d
colnames(X) <- paste0("view", 1:d)
View(X)

## build list of views
data_views <- lapply(1:d, function(j) as.numeric(X[, j]))
names(data_views) <- paste0("view", 1:d)


############################################################
## 2. MCMC settings
############################################################

M       <- 2000   # Gibbs iterations
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
## 4. Build per-view cluster label matrix 
############################################################

table_samples <- res_gibbs$table_of_samples   # list length S, each length n
dish_samples  <- res_gibbs$dish_of_samples    # list length S, each is list per view

S <- length(table_samples)          # number of saved iterations
n <- length(table_samples[[1]])     # number of customers
d <- length(dish_samples[[1]])      # number of views 

## build S x n matrix of cluster labels for a given view v
build_clusters_for_view <- function(v, table_samples, dish_samples) {
  S <- length(table_samples)
  n <- length(table_samples[[1]])
  
  cl_mat <- matrix(NA_integer_, nrow = S, ncol = n)
  
  for (s in 1:S) {
    table_s <- table_samples[[s]]   # length n, table index for each i
    dish_s  <- dish_samples[[s]]    # list length d, each integer vector of length #tables
    
    ## - dish_s[[v]][t_index] is cluster label for table t_index in view v
    
    cl_v <- integer(n)
    for (i in 1:n) {
      t  <- table_s[i]               # table label in 1..T
      cl_v[i] <- dish_s[[v]][t]      # dish index for that table in view v
    }
    
    ## mcclust prefers labels >= 1
    if (min(cl_v) == 0L) {
      cl_v <- cl_v + 1L
    }
    
    cl_mat[s, ] <- cl_v
  }
  
  cl_mat
}


## an example: choose a view to analyse
v_analyze <- 1
clusters_v <- build_clusters_for_view(v_analyze, table_samples, dish_samples)
dim(clusters_v)  # S x n


############################################################
## 5. Posterior Similarity Matrix and representative partition
############################################################

psm_v    <- comp.psm(clusters_v)
cl_v_def <- minVI(psm_v)$cl   # length n, final clusters for view v_analyze

