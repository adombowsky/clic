source("rfuncts/conditional_dependent.R")
source("rfuncts/conditional_independent.R")
source("rfuncts/conditional_trivial.R")
# ari_mat function
# calculating ARI
ari_mat <- function(cmat) {
  c0 = cmat[,1]
  cmat = cmat[,-1]
  ari <- c()
  for (i in 1:ncol(cmat)) {
    ari[i] <- adjustedRandIndex(c0,cmat[,i])
  }
  return(round(ari,3))
} 
R <- 30000
B <- 10000
Th <- 2
c.dep <- conditional_dependent(R=R, B=B, Th=Th)
c.ind <- conditional_independent(R=R, B=B, Th=Th)
c.triv <- conditional_trivial(R=R, B=B, Th=Th)
conditional.results <- list(c.dep = c.dep,
                            c.ind = c.ind,
                            c.triv = c.triv)
saveRDS(conditional.results,"simulations/results/conditional.rds")


