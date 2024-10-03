conditional_independent <- function(n=200, omega0=1/2, seed=1, R=1000, B=500, Th = 1) {
  # n = sample size
  # omega0 = proportion in cluster 1 in C1
  # seed = random number seed
  # R = number of iterations
  # B = burn-in
  # Th = thin
  set.seed(seed)
  require(mcclust)
  require(mcclust.ext)
  require(Rcpp)
  require(RcppArmadillo)
  sourceCpp("rcppfuncts/sampling.cpp")
  sourceCpp("rcppfuncts/postprocessing.cpp")
  source("rfuncts/t_hdp.R")
  cmat <- matrix(sample(2, size = n*2, replace = T, prob = c(omega0, 1-omega0)), nrow = 2)
  c1.true <- cmat[1,]
  c2.true <- cmat[2,]
  theta1 <- c(1, -1)
  theta2 <- -theta1
  X1 <- X2 <- c()
  for (i in 1:n) {
    X1[i] <- rnorm(n=1, mean = theta1[c1.true[i]], sqrt(0.2))
    X2[i] <- rnorm(n=1, mean = X1[i]*theta2[c2.true[i]], sqrt(0.2))
  }
  X <- cbind(X1,X2)
  X_df <- data.frame(X1=X1,X2=X2)
  fit <- gibbs_conditional_markov(n=n,
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
  rand <- rand_index_MCMC(c1,c2)
  ind <- seq(1, length(rho), by = Th)
  rho <- rho[ind]
  c1 = c1[ind,]
  c2 = c2[ind,]
  nu <- 1/(rho+1)
  rand <- rand[ind]
  # point estimate of c1
  psm1 <- mcclust::comp.psm(c1)
  mv1 <- mcclust.ext::minVI(psm1,c1)
  c1.minVI <- mv1$cl
  mclust::adjustedRandIndex(c1.true, c1.minVI)
  # point estimate of c2
  psm2 <- mcclust::comp.psm(c2)
  mv2 <- mcclust.ext::minVI(psm2,c2)
  c2.minVI <- mv2$cl
  mclust::adjustedRandIndex(c2.true, c2.minVI)
  # comparison: conditional t-HDP
  t_hdp <- telescopic_conditional_HDP_NNIG_uni(data=X, hyper = c(1,1,1,1,1), totiter = R)
  ## c1
  c1_t_hdp <- t_hdp[-(1:B),1,]
  c1_t_hdp <- c1_t_hdp[ind,]
  psm1_t_hdp <- mcclust::comp.psm(c1_t_hdp)
  mv1_t_hdp <- mcclust.ext::minVI(psm1_t_hdp)
  c1.t_hdp <- mv1_t_hdp$cl
  mclust::adjustedRandIndex(c1.true, c1.t_hdp)
  ## c2
  c2_t_hdp <- t_hdp[-(1:B),2,]
  c2_t_hdp <- c2_t_hdp[ind,]
  psm2_t_hdp <- mcclust::comp.psm(c2_t_hdp)
  mv2_t_hdp <- mcclust.ext::minVI(psm2_t_hdp)
  c2.t_hdp <- mv2_t_hdp$cl
  mclust::adjustedRandIndex(c2.true, c2.t_hdp)
  rand_t_hdp <- rand_index_MCMC(c1_t_hdp,c2_t_hdp)
  # comparison: separate clusterings
  ## c1
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
  mclust::adjustedRandIndex(c1.true, c1.dp)
  ## c2
  c2_dp <- gibbs_vanilla_profile(n=n,
                         X_out=X[,2],
                         X_covars = X[,1],
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
  mclust::adjustedRandIndex(c2.true, c2.dp)
  c1mat <- cbind(c1.true, c1.minVI, c1.dp)
  c2mat <- cbind(c2.true, c2.minVI, c2.dp)
  # comparison: EM algorithm
  mcl <- Mclust(X[,1],G=1:10)
  c1.mcl <- mcl$classification
  mclust::adjustedRandIndex(c1.true, c1.mcl)
  lmfit <- lm(X2 ~ X1, data = X_df)
  fmfit <- flexmix(X2~-1+X1, data = X_df, k = 2)
  c2.em <- fmfit@cluster
  mclust::adjustedRandIndex(c2.true, c2.em)
  c1mat <- cbind(c1.true, c1.minVI, c1.t_hdp, c1.dp, c1.mcl)
  c2mat <- cbind(c2.true, c2.minVI, c2.t_hdp, c2.dp, c2.em)
  return(list(c1mat = c1mat, c2mat = c2mat, c1.clic = c1, c2.clic = c2, X=X, rand = rand, rand_t_hdp = rand_t_hdp))
  # c1mat = c1 point estimates
  # c2mat = c2 point estimates
  # c1.clic = c1 posterior samples from CLIC
  # c2.clic = c2 posterior samples from CLIC
  # X = simulated data
  # rand = posterior samples of the Rand index for CLIC
  # rand_t_hdp = posterior samples of the Rand index for the t-HDP
}