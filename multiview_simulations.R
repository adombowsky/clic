source("rfuncts/multiview_dependent.R")
source("rfuncts/multiview_independent.R")
source("rfuncts/multiview_trivial.R")
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
R <- 30000 # number of iterations
B <-10000 # burn-in
Th <- 2 # thinning
c.dep <- multiview_dependent(R=R, B=B, Th=Th,eta=sqrt(0.2)) # set eta = sqrt(0.45) for overlap 2
c.ind <- multiview_independent(R=R, B=B, Th=Th,eta=sqrt(0.2)) # set eta = sqrt(0.45) for overlap 2
c.triv <- multiview_trivial(R=R, B=B, Th=Th,eta=sqrt(0.2)) # set eta = sqrt(0.45) for overlap 2
multiview.results <- list(c.dep = c.dep,
                            c.ind = c.ind,
                            c.triv = c.triv)
saveRDS(multiview.results,"simulations/multiview.rds")

