source("rfuncts/conditional_dependent.R")
source("rfuncts/conditional_independent.R")
source("rfuncts/conditional_trivial.R")
require(Rcpp)
require(RcppArmadillo)
require(mclust)
require(mcclust)
require(mcclust.ext)
# ari_mat function
ari_mat <- function(cmat) {
  c0 = cmat[,1]
  cmat = cmat[,-1]
  ari <- c()
  for (i in 1:ncol(cmat)) {
    ari[i] <- adjustedRandIndex(c0,cmat[,i])
  }
  return(round(ari,3))
} 

# running the three different settings
R <- 30000 # number of iterations
B <- 10000 # number of burn-in
Th <- 2 # thinning
c.dep <- conditional_dependent(R=R, B=B, Th=Th)
c.ind <- conditional_independent(R=R, B=B, Th=Th)
c.triv <- conditional_trivial(R=R, B=B, Th=Th)
conditional.results <- list(c.dep = c.dep,
                            c.ind = c.ind,
                            c.triv = c.triv)
saveRDS(conditional.results,"simulations/conditional.rds")

