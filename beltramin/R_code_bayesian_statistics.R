
library(dplyr)

library(ggplot2)
library(cowplot)
library(Rcpp)
library(tictoc)
#setwd('/Users/jr/Desktop/Bayesian Statistics/Project/Code')
Rcpp::sourceCpp('src/multiview_functions_bayesian_project.cpp')
Rcpp::sourceCpp('src/MH_update_bayesian_statistics.cpp')


#############################
### Simulate data ###########
#############################
set.seed(1999)

x <- data.frame(view1 = c(rnorm(100,2,1), rnorm(100,-2,1)), 
                view2 = c(rnorm(50,0,1), rnorm(100,-4,1), rnorm(50,4,1)), 
                tables = sample(c(1:19), size = 200, replace = T))

x <- x %>%
  group_by(tables)%>%
  mutate(dish1 = sample(c(1:4), size = 1, replace = T))%>%
  mutate(dish2 = sample(c(1:4), size = 1, replace = T))%>%
  mutate(view1_sq = view1^2) %>%
  mutate(view2_sq = view2^2)
# Beware not to leave gaps in the enumeration of dishes and tables (i.e, if I have table "3" I must also have
# tables "1" and "2" in my restaurant)
x <- as.data.frame(x)
n <- nrow(x)
true_clust1 <- c(rep(1,100), rep(2,100))
true_clust2 <- c(rep(1,75), rep(2,50), rep(3,75))

##############################
##############################
##############################

## MCMC parameters
nsim <- 200
cluster1 <- matrix(NA, ncol = nsim, nrow = n)
cluster2 <- matrix(NA, ncol = nsim, nrow = n)

## Hyperparameters
gamma_vec <- gamma1_vec <- gamma2_vec <- numeric()
gamma1_vec[1] <- gamma2_vec[1] <- .4 
gamma_vec[1] <- 2   # initial values

sigma_vec <- sigma1_vec <- sigma2_vec <- numeric()
sigma_vec[1] <- sigma1_vec[1] <- sigma2_vec[1] <- .1

# prior per theta|sigma ~ N(mu0, sigma)
mu0 <- 0
##

# sigma ~ Inv.Gamma(alpha0, beta0)
alpha0 <- 2
beta0 <- 1

for (j in 1:nsim) {
  if(j %in% seq(100,nsim, by = 100)){
    cat("iterations = ",j, "\n")
  }
  # tic()
  for (i in 1:n) {
    # select latest hyperparametrs
    gamma <- gamma_vec[j]
    gamma1 <- gamma1_vec[j]
    gamma2 <- gamma2_vec[j]
    sigma <- sigma_vec[j]
    sigma1 <- sigma1_vec[j]
    sigma2 <- sigma2_vec[j]
    
    
    id <- as.vector(x[i,])
    # 
    # relabeling step
    x <- relabel(i, data = x, p1 = id$dish1, p2 = id$dish2, t = id$tables)
    # remove i-th observation
    temp <- x[-i,]

    # how many people sitting at each table
    q_t_meno_i <- count_customers_table(temp)
    # 
    # number of tables
    tavoli <- length(q_t_meno_i)

    # number of unique first and second dishes
    K1_meno_i <- max(temp$dish1)
    K2_meno_i <- max(temp$dish2)
    # 
    # ####
    # # update table i # 
    # # probability of new table
    prob_vec <- numeric()
    f_tav_nuovo <- prob_old1_old2(temp, sigma1, sigma2, K1_meno_i, K2_meno_i, id, mu0, alpha0, beta0) + prob_old1_new2(temp,gamma2,sigma1,sigma2,K1_meno_i,K2_meno_i,id, mu0, alpha0,  beta0) 
    +  prob_new1_old2(temp,gamma1,sigma1,sigma2,K1_meno_i,K2_meno_i,id, mu0, alpha0,  beta0) + prob_new1_new2(temp,gamma1,gamma2,sigma1,sigma2,K1_meno_i,K2_meno_i,id, mu0, alpha0, beta0)
    # 
    prob_vec[1] <- (gamma + sigma*tavoli)/((tavoli + gamma1)*(tavoli + gamma2))*(f_tav_nuovo)
    # 
    # # probability old table
    prob_vec <- c(prob_vec, (q_t_meno_i - sigma)*prob_occupied_table(temp, id, tavoli, mu0, alpha0, beta0))
    # 
    prob_vec <- prob_vec/sum(prob_vec)
    # 
    
    index <- sample.int(tavoli + 1, size = 1, prob = prob_vec)
    
    if(index == 1){ 
      # if I select a new table, increase the number of tables in the restaurant
      x$tables[i] <- tavoli + 1
      
      
      # choose first dish for new table
      prob_k1 <- prob_update_dish1(temp, gamma1, sigma1, K1_meno_i, id, mu0, alpha0, beta0)
      prob_k1 <- prob_k1/sum(prob_k1)
      index_1 <- sample.int(K1_meno_i + 1, size = 1, prob = prob_k1)
      if(index_1 == 1){
        x$dish1[i] <- K1_meno_i + 1 # new
      }else{
        x$dish1[i] <- index_1 - 1 # old
      }
      
      
      
      # choose second dish for new table
      prob_k2 <- prob_update_dish2(temp, gamma2, sigma2, K2_meno_i, id, mu0, alpha0, beta0)
      prob_k2 <- prob_k2/sum(prob_k2)
      index_2 <- sample.int(K2_meno_i + 1, size = 1, prob = prob_k2)
      if(index_2 == 1){
        x$dish2[i] <- K2_meno_i + 1 # new
      }else{
        x$dish2[i] <- index_2 - 1   # old
      }
    }else{
      # If I choose an old table, I inherit its dishes
      piatto <- p_t_allocation(temp, tavoli)
      piatto1 <- piatto$dish1
      piatto2 <- piatto$dish2
      x$tables[i] <- index - 1 # remember, I if sample index = 2 it means I am selecting table labeled "1"
      x$dish1[i] <- piatto1[index - 1]
      x$dish2[i] <- piatto2[index - 1]
    }
    cluster1[,j] <- x$dish1
    cluster2[,j] <- x$dish2
  }
  ###
  # update hyperparameters
 
  ###
  # update sigma
  t <- max(x$tables)
  K1 <- max(x$dish1)
  K2 <- max(x$dish2)
  
  # prior for: sigma ~ Beta(2, 5), sigma1 ~ Beta(2, 5), sigma2 ~ Beta(2, 5)
  sigma_vec[j+1] <- update_sigma(x, gamma, sigma, t, K1, K2, 2, 5)
  sigma1_vec[j+1] <- update_sigma1(x, gamma1, sigma1, K1, 2, 5)
  sigma2_vec[j+1] <- update_sigma2(x, gamma2, sigma2, K2, 2, 5)
  
  ##
  
  # update gamma
  # prior for gamma: gamma ~ Gamma(1,1), gamma1 ~ Gamma(1,1), gamma2 ~ Gamma(1,1)
  gamma_vec[j+1] <- update_gamma(gamma_vec[j], sigma_vec[j+1], t, n, 1, 1)
  
  gamma1_vec[j+1] <- update_gamma1(x, gamma1_vec[j], sigma1_vec[j+1], n, 1, 1)
  
  gamma2_vec[j+1] <- update_gamma2(x, gamma2_vec[j], sigma2_vec[j+1], n, 1, 1)
  
}


n_burnin <- 15000
library(coda)
mcmc_object_sig <- mcmc(ts(sigma_vec[(n_burnin+1):nsim]))
mcmc_object_gamma <- mcmc(ts(gamma_vec[(n_burnin+1):nsim]))
mcmc_object_gamma1 <- mcmc(ts(gamma1_vec[(n_burnin+1):nsim]))
mcmc_object_gamma2 <- mcmc(ts(gamma2_vec[(n_burnin+1):nsim]))
mcmc_object_sig1 <- mcmc(ts(sigma1_vec[(n_burnin+1):nsim]))
mcmc_object_sig2 <- mcmc(ts(sigma2_vec[(n_burnin+1):nsim]))

mcmc_obj <- mcmc.list(list(mcmc_object_gamma, mcmc_object_sig))
gelman.diag(mcmc_obj)
effectiveSize(mcmc_object_gamma)
effectiveSize(mcmc_object_sig)
effectiveSize(mcmc_object_gamma1)
effectiveSize(mcmc_object_gamma2)
autocorr.plot(mcmc_object_gamma2)
autocorr.plot(mcmc_object_sig)

plot(ts(gamma_vec[(n_burnin+1):nsim]))
acf(ts(gamma_vec[(n_burnin+1):nsim]))
plot(ts(sigma1_vec[(n_burnin+1):nsim]))
acf(ts(sigma_vec[n_burnin+1:nsim]))
cluster1_1 <- cluster1[,seq((n_burnin + 1), nsim, by=2)]
cluster2_2 <- cluster2[,seq((n_burnin + 1), nsim, by=2)]


library(mcclust)
nsim <- 30000
n_burnin <- 20000
cluster1 <- t(cluster1[,(n_burnin:nsim)])
psm1 <- comp.psm(cluster1)
plotpsm(psm1)
cluster1_def <- minVI(psm1)$cl
help("minVI")
# help("comp.psm")
cluster2 <- t(cluster2[,(n_burnin:nsim)])
psm2 <- comp.psm(cluster2)
cluster2_def <- minVI(psm2)$cl
plotpsm(psm2)
help("minVI")
x$cluster1_def <- cluster1_def
table(cluster1_def, true_clust1)
x$cluster2_def <- cluster2_def
table(cluster2_def, true_clust2)
x$cluster1_def <- ifelse(x$cluster1_def == 1, 2, 1)
x$cluster2_def <- ifelse(x$cluster2_def == 1, 2, 1)

x$true_clust1 <- true_clust1
x$true_clust2 <- true_clust2
g1 <- ggplot(data = x, mapping = aes(x = view1, y = view2)) +
  geom_point(aes(shape = as.factor(cluster2_def), colour = as.factor(cluster1_def)), size = 4) +
  geom_point(colour = "black", size = 2) +
  geom_point(data = x[which(x$true_clust1 != x$cluster1_def),], color = "red", size = 5, shape = 1, stroke = 2) +
  geom_point(data = x[which(x$true_clust2 != x$cluster2_def),], color = "red", size = 5, shape = 1, stroke = 2) +
  guides(
    colour = guide_legend(
      title = "Estimated cluster view 1",
      order = 1,
      position = "top"
    ),
    shape = guide_legend(
      title = "Estimated cluster view 2",
      position = "top"
    )
  ) +
  theme_minimal()

cowplot::plot_grid(g1, g2)
