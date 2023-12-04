sample_CLIC <- function(n, L, gamma1, gamma2, rho) {
  # sample dependent partitions from the CLIC model
  sourceCpp("rcppfuncts/sampling.cpp")
  q1 <- rdirichlet_arma(rep(gamma1/L,L))
  q2 <- rdirichlet_arma(rep(gamma2/L, L))
  # simulate a
  q <- rho * tcrossprod(q1, q2)
  qvec <- as.vector(q) # goes by column
  pvec <- rdirichlet_arma(qvec)
  p <- matrix(pvec, byrow = FALSE, nrow = L, ncol = L)
  p1 <- rowSums(p)
  # simulate s and sprime
  c1 <- c2 <- c()
  for (i in 1:n) {
    c1[i] <- sample(L, size = 1, replace = TRUE, prob = p1)
    c2[i] <- sample(L, size = 1, replace = TRUE, prob = p[c1[i],]/p1[c1[i]])
  }
  return(list(c1 = c1, c2 = c2))
}