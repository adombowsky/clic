samp_size_and_dim <- function(n=100, p1=2, p2=2, omega0=1/2, R=1000, B=500, Th = 1,eta=sqrt(0.2), seed=1){
  # n = sample size
  # p1 = view 1 dimension
  # p2 = view 2 dimension
  # omega0 = proportion in cluster 1 in C1
  # seed = random number seed
  # R = number of iterations
  # B = burn-in
  # Th = thin
  # eta = overlap
  require(Rcpp)
  require(RcppArmadillo)
  require(mclust)
  require(mcclust)
  require(mcclust.ext)
  sourceCpp("rcppfuncts/sampling.cpp")
  source("rfuncts/t_hdp_multivariate.R")
  p <- p1 + p2
  # sampling data
  set.seed(seed)
  c1.true <- c2.true <- c()
  theta1 <- matrix(0,nrow=2,ncol=p1)
  theta1[1,] <- c(1, -1)[rep(1:2,length=p1)]
  theta1[2,] <- c(-1, 1)[rep(1:2,length=p1)]
  theta2 <- matrix(0,nrow=2,ncol=p2)
  theta2[1,] <- c(-1, 1)[rep(1:2,length=p2)]
  theta2[2,] <- c(1, -1)[rep(1:2,length=p2)]
  X1 <- matrix(0,nrow=n,ncol=p1)
  X2 <- matrix(0,nrow=n,ncol=p2)
  for (i in 1:n) {
    u <- runif(1)
    if (u < 8/10) {
      c1.true[i] = c2.true[i] = sample(2,1,T,c(omega0,1-omega0))
    } else {
      c1.true[i] = sample(2,1,T,c(omega0,1-omega0))
      c2.true[i] = sample(2,1,T,c(omega0,1-omega0))
    }
    X1[i,] <-   theta1[c1.true[i],] + eta*rnorm(p1) 
    X2[i,] <- theta2[c2.true[i],] + eta*rnorm(p2) 
  }
  X <- cbind(X1,X2)
  #pairs(X)
  
  # fitting model
  a.clic <- Sys.time()
  fit <- gibbs_markov_multivariate(n=n,
                                   X1=scale(X1),
                                   X2=scale(X2),
                                   gamma1 = 1,
                                   gamma2 = 1,
                                   mu01=rep(0,p1),
                                   sigma01=diag(1,p1),
                                   mu02=rep(0,p2),
                                   sigma02=diag(1,p2),
                                   alpha1=1,
                                   beta1=1,
                                   alpha2=1,
                                   beta2=1,
                                   L1 = 10,
                                   L2 = 10,
                                   a_rho = 1,
                                   b_rho = 1,
                                   R = R)
  # burn-in and thinning
  c1 <- fit$c1[-(1:B),]
  c2 <- fit$c2[-(1:B),]
  ind <- seq(1, nrow(c1), by = Th)
  c1 = c1[ind,]
  c2 = c2[ind,]
  #rand <- mapply(fossil::rand.index, asplit(c1,1), asplit(c2,1))
  # point estimate of c1
  psm1 <- mcclust::comp.psm(c1)
  mv1 <- mcclust.ext::minVI(psm1,c1)
  c1.minVI <- mv1$cl
  # point estimate of c2
  psm2 <- mcclust::comp.psm(c2)
  mv2 <- mcclust.ext::minVI(psm2,c2)
  c2.minVI <- mv2$cl
  b.clic <- Sys.time()
  t.clic <- round(as.numeric(difftime(b.clic,a.clic,units="secs")),3)
  
  # using t-HDP
  a.t_hdp <- Sys.time()
  t_hdp <- telescopic_HDP_NNIW_multi(data=X,
                                     colayer = c(rep(1,p1), rep(2,p2)),
                                     totiter = R)
  ##c1
  c1_t_hdp <- t_hdp[,1,]
  c1_t_hdp <- c1_t_hdp[ind,]
  psm1_t_hdp <- mcclust::comp.psm(c1_t_hdp)
  mv1_t_hdp <- mcclust.ext::minVI(psm1_t_hdp)
  c1.t_hdp <- mv1_t_hdp$cl
  ## c2
  c2_t_hdp <- t_hdp[,2,]
  c2_t_hdp <- c2_t_hdp[ind,]
  psm2_t_hdp <- mcclust::comp.psm(c2_t_hdp)
  mv2_t_hdp <- mcclust.ext::minVI(psm2_t_hdp)
  c2.t_hdp <- mv2_t_hdp$cl
  b.t_hdp <- Sys.time()
  t.t_hdp <- round(as.numeric(difftime(b.t_hdp,a.t_hdp,units="secs")),3)
  # competiror: mclust on marginals
  ## c1
  a.mcl <- Sys.time()
  mcl1 <- Mclust(scale(X1),G=1:10)
  c1.mcl <- mcl1$classification
  #mclust::adjustedRandIndex(c1.true, c1.mcl)
  ## c2
  mcl2 <- Mclust(scale(X2),G=1:10)
  c2.mcl <- mcl2$classification
  b.mcl <- Sys.time()
  t.mcl <- round(as.numeric(difftime(b.mcl,a.mcl,units="secs")),3)
  #mclust::adjustedRandIndex(c2.true, c2.mcl)
  # competitor 3: independent DPs
  ## c1
  a.dp <- Sys.time()
  c1_dp <- gibbs_vanilla_multivariate(n=n,
                                      X=scale(X1),
                                      gam=1,
                                      mu0 = rep(0,p1),
                                      sigma0 = diag(1,p1),
                                      alpha=1,
                                      beta=1,
                                      L=10,
                                      R=R)
  c1_dp <- c1_dp[-(1:B),]
  c1_dp <- c1_dp[ind,]
  psm1_dp <- mcclust::comp.psm(c1_dp)
  mv1_dp <- mcclust.ext::minVI(psm1_dp)
  c1.dp <- mv1_dp$cl
  #mclust::adjustedRandIndex(c1.true, c1.dp)
  ## c2
  c2_dp <- gibbs_vanilla_multivariate(n=n,
                                      X=scale(X2),
                                      gam=1,
                                      mu0 = rep(0,p2),
                                      sigma0 = diag(1,p2),
                                      alpha=1,
                                      beta=1,
                                      L=10,
                                      R=R)
  c2_dp <- c2_dp[-(1:B),]
  c2_dp <- c2_dp[ind,]
  psm2_dp <- mcclust::comp.psm(c2_dp)
  mv2_dp <- mcclust.ext::minVI(psm2_dp)
  c2.dp <- mv2_dp$cl
  b.dp <- Sys.time()
  t.dp <- round(as.numeric(difftime(b.dp,a.dp,units="secs")),3)
  # vanilla joint model
  a.jdp <- Sys.time()
  c_vjm <- gibbs_vanilla_multivariate(n=n,
                                      X=X,
                                      gam=1,
                                      mu0 = rep(0,p),
                                      sigma0 = diag(1,p),
                                      alpha=1,
                                      beta=1,
                                      L=10,
                                      R=R)
  psm_vjm <- mcclust::comp.psm(c_vjm)
  mv_vjm <- mcclust.ext::minVI(psm_vjm,c_vjm)
  c.vjm <- mv_vjm$cl
  b.jdp <- Sys.time()
  t.jdp <- round(as.numeric(difftime(b.jdp,a.jdp,units="secs")),3)
  
  # collecting all info
  c1mat <- cbind(c1.true, c1.minVI, c1.t_hdp, c1.mcl, c1.dp,c.vjm)
  c2mat <- cbind(c2.true, c2.minVI, c2.t_hdp, c2.mcl, c2.dp,c.vjm)
  tvec <- c(t.clic, t.t_hdp, t.mcl, t.dp,t.jdp)
  return(list(c1 = c1mat, c2 = c2mat, t = tvec))
  # c1mat = c1 point estimates
  # c2mat = c2 point estimates
  # tvec = computation times
}