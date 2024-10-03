multiview_independent <- function(n=200, omega0=1/2,seed=1, R=1000, B=500, Th = 1, eta=sqrt(0.2)) {
  # n = sample size
  # omega0 = proportion in cluster 1 in C1
  # seed = random number seed
  # R = number of iterations
  # B = burn-in
  # Th = thin
  # eta = view 2 overlap
  set.seed(seed)
  require(mclust)
  require(mcclust)
  require(mcclust.ext)
  require(Rcpp)
  require(RcppArmadillo)
  sourceCpp("rcppfuncts/sampling.cpp")
  source("rfuncts/t_hdp.R")
  require(dirichletprocess)
  cmat <- matrix(sample(2, size = n*2, replace = T, prob = c(omega0, 1-omega0)), nrow = 2)
  c1.true <- cmat[1,]
  c2.true <- cmat[2,]
  theta1 <- c(1, -1)
  theta2 <- -theta1
  X1 <- X2 <- c()
  for (i in 1:n) {
    u <- runif(1)
    X1[i] <- rnorm(n=1, mean = theta1[c1.true[i]], sqrt(0.2))
    X2[i] <- rnorm(n=1, mean = theta2[c2.true[i]], eta)
  }
  X <- cbind(X1,X2)
  a.clic <- Sys.time()
  fit <- gibbs_markov(n=n,
               X=X,
               gamma1 = 1,
               gamma2 = 1,
               mu01=0,
               sigma01=1,
               mu02=0,
               sigma02=1,
               alpha1=1,
               beta1=1,
               alpha2=1,
               beta2=1,
               L1 = 5,
               L2 = 5,
               a_rho = 1,
               b_rho = 1,
               R = R)
  # burn-in and thinning
  rho <- fit$rho[-(1:B)]
  c1 <- fit$c1[-(1:B),]
  c2 <- fit$c2[-(1:B),]
  ind <- seq(1, nrow(c1), by = Th)
  rho <- rho[ind]
  c1 = c1[ind,]
  c2 = c2[ind,]
  nu <- 1/(rho+1)
  # point estimate of c1
  psm1 <- mcclust::comp.psm(c1)
  mv1 <- mcclust.ext::minVI(psm1,c1)
  c1.minVI <- mv1$cl
  #mclust::adjustedRandIndex(c1.true, c1.minVI)
  # point estimate of c2
  psm2 <- mcclust::comp.psm(c2)
  mv2 <- mcclust.ext::minVI(psm2,c2)
  c2.minVI <- mv2$cl
  b.clic <- Sys.time()
  t.clic <- round(as.numeric(difftime(b.clic,a.clic,units="secs")),3)
  rand <- rand_index_MCMC(c1,c2)
  #mclust::adjustedRandIndex(c2.true, c2.minVI)
  # distribution of nu and rand
  #hist(nu)
  #hist(rand)
  # competitor 1: t-HDP
  a.t_hdp <- Sys.time()
  t_hdp <- telescopic_HDP_NNIG_uni(data=X, hyper = c(1,1,1,1,1), totiter=R)
  ## c1
  c1_t_hdp <- t_hdp[-(1:B),1,]
  c1_t_hdp <- c1_t_hdp[ind,]
  psm1_t_hdp <- mcclust::comp.psm(c1_t_hdp)
  mv1_t_hdp <- mcclust.ext::minVI(psm1_t_hdp)
  c1.t_hdp <- mv1_t_hdp$cl
  #mclust::adjustedRandIndex(c1.true, c1.t_hdp)
  ## c2
  c2_t_hdp <- t_hdp[-(1:B),2,]
  c2_t_hdp <- c2_t_hdp[ind,]
  psm2_t_hdp <- mcclust::comp.psm(c2_t_hdp)
  mv2_t_hdp <- mcclust.ext::minVI(psm2_t_hdp)
  c2.t_hdp <- mv2_t_hdp$cl
  b.t_hdp <- Sys.time()
  t.t_hdp <- round(as.numeric(difftime(b.t_hdp,a.t_hdp,units="secs")),3)
  #mclust::adjustedRandIndex(c2.true, c2.t_hdp)
  rand_t_hdp <- mapply(fossil::rand.index, asplit(c1_t_hdp,1), asplit(c2_t_hdp,1))
  # competiror 2: mclust on marginals
  ## c1
  a.mcl <- Sys.time()
  mcl1 <- Mclust(X[,1],G=1:5)
  c1.mcl <- mcl1$classification
  #mclust::adjustedRandIndex(c1.true, c1.mcl)
  ## c2
  mcl2 <- Mclust(X[,2],G=1:5)
  c2.mcl <- mcl2$classification
  b.mcl <- Sys.time()
  t.mcl <- round(as.numeric(difftime(b.mcl,a.mcl,units="secs")),3)
  #mclust::adjustedRandIndex(c2.true, c2.mcl)
  # competitor 3: independent DPs
  ## c1
  a.dp <- Sys.time()
  c1_dp <- gibbs_vanilla(n=n,
                         X=X[,1],
                         gam=1,
                         mu0=0,
                         sigma=1,
                         alpha=1,
                         beta=1,
                         L=10,
                         R=R)[-(1:B),]
  psm1_dp <- mcclust::comp.psm(c1_dp)
  mv1_dp <- mcclust.ext::minVI(psm1_dp)
  c1.dp <- mv1_dp$cl
  #mclust::adjustedRandIndex(c1.true, c1.dp)
  ## c2
  c2_dp <- gibbs_vanilla(n=n,
                         X=X[,2],
                         gam=1,
                         mu0=0,
                         sigma=1,
                         alpha=1,
                         beta=1,
                         L=10,
                         R=R)[-(1:B),]
  psm2_dp <- mcclust::comp.psm(c2_dp)
  mv2_dp <- mcclust.ext::minVI(psm2_dp)
  c2.dp <- mv2_dp$cl
  b.dp <- Sys.time()
  t.dp <- round(as.numeric(difftime(b.dp,a.dp,units="secs")),3)
  #mclust::adjustedRandIndex(c2.true, c2.dp)
  c1mat <- cbind(c1.true, c1.minVI, c1.t_hdp, c1.mcl, c1.dp)
  c2mat <- cbind(c2.true, c2.minVI, c2.t_hdp, c2.mcl, c2.dp)
  tvec <- c(t.clic, t.t_hdp, t.mcl, t.dp)
  return(list(c1mat = c1mat, c2mat = c2mat, c1.clic = c1, c2.clic = c2, X=X, rand=rand, rand_t_hdp = rand_t_hdp, t = tvec))
  # c1mat = c1 point estimates
  # c2mat = c2 point estimates
  # c1.clic = c1 posterior samples from CLIC
  # c2.clic = c2 posterior samples from CLIC
  # X = simulated data
  # rand = posterior samples of the Rand index for CLIC
  # rand_t_hdp = posterior samples of the Rand index for the t-HDP
  # tvec = computation times
}