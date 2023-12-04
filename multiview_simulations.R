source("rfuncts/multiview_dependent.R")
source("rfuncts/multiview_independent.R")
source("rfuncts/multiview_trivial.R")
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
c.dep <- multiview_dependent(R=R, B=B, Th=Th)
c.ind <- multiview_independent(R=R, B=B, Th=Th)
c.triv <- multiview_trivial(R=R, B=B, Th=Th)
multiview.results <- list(c.dep = c.dep,
                            c.ind = c.ind,
                            c.triv = c.triv)
saveRDS(multiview.results,"simulations/results/multiview.rds")

