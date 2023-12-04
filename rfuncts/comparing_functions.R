cramer <- function(s_and_sprime) {
  n <- length(s_and_sprime)/2
  s <- s_and_sprime[1:n]
  sprime <- s_and_sprime[-(1:n)]
  cV <- cramerV(x=s,y=sprime)
  return(cV)
}

rand <- function(s_and_sprime) {
  n <- length(s_and_sprime)/2
  s <- s_and_sprime[1:n]
  sprime <- s_and_sprime[-(1:n)]
  rand <- rand.index(group1 = s, group2 = sprime)
  return(rand)
}

expectedrand <- function(rho, gamma1, gamma2) {
  nu <- 1/(rho+1)
  eri <- nu + (1-nu)*(1 + gamma1*gamma2)/((gamma1+1) * (gamma2+1))
  return(eri)
}