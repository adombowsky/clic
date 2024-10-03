## varying sample size and dimension ##
source("rfuncts/samp_size_and_dim.R")
require(mclust)
require(mcclust)
require(mcclust.ext)
require(Rcpp)
require(RcppArmadillo)
library(foreach)
library(doParallel)

# securing cores
n.cores <- 9
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

# running sims
p1 <- 2 # dimension of view 1
P2 <- c(2,10,25) # dimension of view 2
N <- c(100,200,500) # sample size
sims <- foreach (i = 1:length(N)) %:%
  foreach (j = 1:length(P2),
           .packages = c("Rcpp", "RcppArmadillo", "mclust", "mcclust", "mcclust.ext", "progress", "mvtnorm", "LaplacesDemon"),
           .noexport = c("gibbs_markov_multivariate", "gibbs_vanilla_multivariate"))  %dopar% {
  samp_size_and_dim(n=N[i], # sample size
                    p1 = p1, # dimension (view 1)
                    p2 = P2[j], # dimension (view 2)
                    omega0 = 1/3, # proportion in cluster 1
                    R=30000, # iterations
                    B=5000, # burn-in
                    Th=2, # thin
                    eta=sqrt(0.4), # overlap
                    seed=1) # random number seed
           }
# computing results
# ARI matrix
ari_mat <- function(cmat) {
  c0 = cmat[,1]
  cmat = cmat[,-1]
  ari <- c()
  for (i in 1:ncol(cmat)) {
    ari[i] <- adjustedRandIndex(c0,cmat[,i])
  }
  return(round(ari,3))
}

# ARIs
## n = 100
ari_mat(sims[[1]][[1]]$c1) # p = 2
ari_mat(sims[[1]][[1]]$c2)

ari_mat(sims[[1]][[2]]$c1) # p = 10
ari_mat(sims[[1]][[2]]$c2)

ari_mat(sims[[1]][[3]]$c1) # p = 25
ari_mat(sims[[1]][[3]]$c2)
## n = 200
ari_mat(sims[[2]][[1]]$c1) # p = 2
ari_mat(sims[[2]][[1]]$c2)

ari_mat(sims[[2]][[2]]$c1) # p = 10
ari_mat(sims[[2]][[2]]$c2)

ari_mat(sims[[2]][[3]]$c1) # p = 25
ari_mat(sims[[2]][[3]]$c2)
## n = 500
ari_mat(sims[[3]][[1]]$c1) # p = 2
ari_mat(sims[[3]][[1]]$c2)

ari_mat(sims[[3]][[2]]$c1) # p = 10
ari_mat(sims[[3]][[2]]$c2)

ari_mat(sims[[3]][[3]]$c1) # p = 25
ari_mat(sims[[3]][[3]]$c2)

# Cluster sizes
##n=100 
apply(sims[[1]][[1]]$c1[,-1], 2, function(x) length(table(x))) # p = 2
apply(sims[[1]][[1]]$c2[,-1], 2, function(x) length(table(x)))

apply(sims[[1]][[2]]$c1[,-1], 2, function(x) length(table(x))) # p = 10
apply(sims[[1]][[2]]$c2[,-1], 2, function(x) length(table(x)))

apply(sims[[1]][[3]]$c1[,-1], 2, function(x) length(table(x))) # p = 25
apply(sims[[1]][[3]]$c2[,-1], 2, function(x) length(table(x)))
##n=200 
apply(sims[[2]][[1]]$c1[,-1], 2, function(x) length(table(x))) # p = 2
apply(sims[[2]][[1]]$c2[,-1], 2, function(x) length(table(x)))

apply(sims[[2]][[2]]$c1[,-1], 2, function(x) length(table(x))) # p = 10
apply(sims[[2]][[2]]$c2[,-1], 2, function(x) length(table(x)))

apply(sims[[2]][[3]]$c1[,-1], 2, function(x) length(table(x))) # p = 25
apply(sims[[2]][[3]]$c2[,-1], 2, function(x) length(table(x)))
##n=500 
apply(sims[[3]][[1]]$c1[,-1], 2, function(x) length(table(x))) # p = 2
apply(sims[[3]][[1]]$c2[,-1], 2, function(x) length(table(x)))

apply(sims[[3]][[2]]$c1[,-1], 2, function(x) length(table(x))) # p = 10
apply(sims[[3]][[2]]$c2[,-1], 2, function(x) length(table(x)))

apply(sims[[3]][[3]]$c1[,-1], 2, function(x) length(table(x))) # p = 25
apply(sims[[3]][[3]]$c2[,-1], 2, function(x) length(table(x)))

# computation times
##n=100
round(sims[[1]][[1]]$t,3) # p = 2
round(sims[[1]][[2]]$t,3) # p = 10
round(sims[[1]][[3]]$t,3) # p = 25
##n=200
round(sims[[2]][[1]]$t,3) # p = 2
round(sims[[2]][[2]]$t,3) # p = 10
round(sims[[2]][[3]]$t,3) # p = 25
##n=500
round(sims[[3]][[1]]$t,3) # p = 2
round(sims[[3]][[2]]$t,3) # p = 10
round(sims[[3]][[3]]$t,3) # p = 25
