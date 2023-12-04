#t-HDP _ conditional algorithm - NNIG model - univariate
library(progress)#to draw the progress bar

#markovian dependence for L layers, univariate layers
#normal kernel and normal-inverse-gamma base
#the mean of the base is fixed to 0
#hyper are the 2 concentration parameters, the 3 param k0, alpha, beta of the NIG base
#if alpharandom=TRUE the two first entries of hyper are only the initialization
#data is the dataset nXL 
#H0 and H are the number of mixture components used in the approximation 
#of the infinite mixtures

telescopic_conditional_HDP_NNIG_uni = function(data, hyper=c(0.1, 0.1, 0.1, 0.1, 0.1),
                                   alpharandom = FALSE,
                                   H0 = 5, H = 5, totiter = 1000){
  #has to work for many X and one param
  kernel_eval_joint = function(x, hyper1, hyper2){
    return(sum(dnorm(x, hyper1, sqrt(hyper2), log = TRUE)))
  }
  
  #has to work for many X and one param
  kernel_eval_joint_outcome = function(y, x, hyper1, hyper2) {
    return(sum(dnorm(y, x*hyper1, sqrt(hyper2), log = TRUE)))
  }
  
  #has to work for one X and many param
  kernel_eval = function(x, hyper1, hyper2){
    return(dnorm(x, hyper1, sqrt(hyper2), log = TRUE))
  }
  
  #has to work for one X and many param
  kernel_eval_outcome = function(y, x, hyper1, hyper2){
    return(dnorm(y, x*hyper1, sqrt(hyper2), log = TRUE))
  }
  
  posterior_sampler = function(x, k0 = 0.1, a = 0.1, b = 0.1){
    n_temp = length(x)
    if(n_temp>0){
      xbar = mean(x)
      sum_squared = sum((x - xbar)^2)
      sigma = 1 / rgamma(1, a + n_temp / 2, b + sum_squared/2 +
                           xbar^2/2 * n*k0 / (k0+n_temp))
      mu = rnorm(1, mean = n_temp / (n_temp + k0) * xbar,
                 sd = sqrt(sigma / (n_temp+k0)))
      return( list(mu = mu, sigma = sigma) )
    }else{
      sigma = 1 / rgamma(1, a , b )
      mu = rnorm(1, mean = 0,
                 sd = sqrt(sigma / k0))
      return( list(mu = mu, sigma = sigma) )
    }
  }
  
  posterior_sampler_outcome = function(y, x, k0 = 0.1, a = 0.1, b = 0.1){
    n_temp = length(y)
    if(n_temp>0){
      beta_hat = crossprod(x,y)/crossprod(x,x)
      k_n = (crossprod(x,x) + k0)
      beta_n = crossprod(x,x)*beta_hat/k_n
      a_n = a + n_temp/2
      b_n = b + 0.5*(crossprod(y,y) - beta_n^2*k_n)
      sigma = 1 / rgamma(1, a_n, b_n)
      mu = rnorm(1, mean = beta_n,
                 sd = sqrt(sigma/k_n))
      return( list(mu = mu, sigma = sigma) )
    }else{
      sigma = 1 / rgamma(1, a , b )
      mu = rnorm(1, mean = 0,
                 sd = sqrt(sigma / k0))
      return( list(mu = mu, sigma = sigma) )
    }
  }
  
  
  
  L = dim(data)[2] #number of layers
  n_tot = dim(data)[1] #number of items
  
  #the data in a matrix X, which is nxL
  X = data 
  
  #RANDOM PROB.S
  pim = array(NA, dim = c(H0, H, L)) #weights (in log scale)
  b = array(NA, dim = c(H0, H, L)) #sticks for pi
  
  pi0 = matrix(NA, nrow = H0, ncol = L) #weights of the common overall in HDP
  b0 = matrix(NA, nrow = H0, ncol = L) #sticks for pi0 (in log scale)
  theta0 = matrix(NA, nrow = H0, ncol = L) #atoms
  sigma0 = matrix(NA, nrow = H0, ncol = L) #atoms #var
  
  #PARTITION 
  m = matrix(NA, nrow = L, ncol = n_tot) #clustering configuration
  c = array(NA, dim = c(L, n_tot)) #auxiliary: tables per item
  k = array(0, dim = c(L, H0*H)) #auxiliary: dishes per table
  q = matrix(NA, nrow = L, ncol = H0) #dish freq
  n = matrix(NA, nrow = L, ncol = H*H0)#table freq
  nbar = array(NA, dim = c(L, H0, H))#table freq remapped
  
  m_saved = array(NA,c(totiter, L, n_tot))
  
  #initialize
  for(l in 1:L){
    m[l, ] = stats::kmeans(as.matrix(X[,l]),min(10,H0))$cluster
  }
  c[1, ] = m[1, ]
  for(l in 2:L){
    c[l, ] = (m[l-1,]-1)*H + m[l, ]
  }
  for(l in 1:L){
    for (cc in 1:(H*H0)){
      if (cc %in% c[l,]){
        k[l,cc] = m[l,c[l,]==cc][1]
      }else{
        k[l,cc] = sample(1:H0, 1)
      }
    }
  }
  
  
  q[,] = 0
  for(l in 1:L){
    
    temp = rep(0, H*H0)
    temp[unique(c[l,])] = 1
    for(h in 1:H0){
      q[l,h] = sum((k[l,]==h)*temp)
    }
    
    for(cc in 1:(H*H0)){
      n[l,cc] = sum(c[l,]==cc)
    }
    
    for(h1 in 1:H0){
      for(h2 in 1:H){
        nbar[l, h1, h2] = n[l, (h1-1)*H + h2]
      }
    }
  }
  
  if(!alpharandom){
    alpha0 = hyper[1]; alpha = hyper[2]
    k0 = hyper[3]; a_sigma = hyper[4]; b_sigma = hyper[5]
    pb <- progress_bar$new(
      format = " MCMC [:bar] :percent Estimated completion time: :eta",
      total = totiter, clear = FALSE, width= 100)
    #mcmc
    for (iter in 1:totiter){
      pb$tick()
      for (l in 1:L){
        ql_sum = c(rev(cumsum(rev(q[l,])))[2:H0],0)
        for (h1 in 1:H0){
          b0[h1, l] = rbeta(1, 1 + q[l,h1], alpha0 + ql_sum[h1])
          if(h1>1){
            pi0[h1, l] = log(b0[h1, l]) + sum(log(1 - b0[1:(h1-1), l]))
          }else{
            pi0[h1, l] = log(b0[h1, l])
          }
          if (l==1) { # covariates
            out = posterior_sampler(X[m[l,]==h1,l], k0, a_sigma, b_sigma)
          } else { # outcome
            out = posterior_sampler_outcome(X[m[2,]==h1, 2],
                                            X[m[2,]==h1,1],
                                            k0, a_sigma, b_sigma)
          }
          
          theta0[h1, l] = out$mu
          sigma0[h1, l] = out$sigma
          
          nbar_sum = c(rev(cumsum(rev(nbar[l,h1,])))[2:H],0)
          for(h2 in 1:H){
            b[h1, h2, l] = rbeta(1, 1 + nbar[l,h1,h2], alpha + nbar_sum[h2])
            if(h2>1){
              pim[h1, h2, l] = log(b[h1, h2, l]) + sum(log(1 - b[h1, 1:(h2-1),l]))
            }else{
              pim[h1, h2, l] = log(b[h1, h2, l])
            }
          }
        }
      }
      for (l in 1:L){
        for (i in 1:n_tot){
          if(l>1){past = m[l-1, i]}else{past = 1}
          fut = rep(0,H)
          if(l<L){
            for (cc in 1:H){
              for (d in 1:H0){
                fut[cc] = fut[cc] + 
                  (k[l+1,(k[l,cc]-1)*H + d] == m[l+1,i])*exp(pim[k[l,cc],d,l+1])
              } 
            }
            fut = log(fut)
          }
          if (l==1) {
            prob = (pim[past,,l]) + kernel_eval(X[i,l], theta0[ ,l], sigma0[, l]) + fut
          } else {
            prob = (pim[past,,l]) + kernel_eval_outcome(X[i,2], X[i,1], theta0[ ,l], sigma0[, l]) + fut
          }
          if(max(prob)==-Inf){prob[]=1}
          prob = prob - max(prob)
          prob = exp(prob) 
          if(sum(exp(fut))!=0){ #if fut is zero it cannot move
            c[l, i] = sample(((past-1)*H+1):(past*H), 1, prob = prob)
          }
        }
        for (cc in 1:(H*H0)){
          prob = NULL
          for(kk in 1:H0){
            if (l==1) {
              prob[kk] = pi0[kk,l] + kernel_eval_joint(X[c[l,]==cc,l], theta0[kk,l], sigma0[kk,l])   
            } else {
              prob[kk] = pi0[kk,l] + kernel_eval_joint_outcome(X[c[2,]==cc,2], X[c[2,]==cc,1], theta0[kk,l], sigma0[kk,l])   
            }
            
          }
          prob = prob - max(prob)
          prob = exp(prob) 
          k[l, cc] = sample(1:H0, 1, prob = prob)
        }
        for (i in 1:n_tot){
          m[l,i] = k[l,c[l,i]]
        }
      } 
      q[,] = 0
      for(l in 1:L){
        
        temp = rep(0, H*H0)
        temp[unique(c[l,])] = 1
        for(h in 1:H0){
          q[l,h] = sum((k[l,]==h)*temp)
        }
        
        for(cc in 1:(H*H0)){
          n[l,cc] = sum(c[l,]==cc)
        }
        
        for(h1 in 1:H0){
          for(h2 in 1:H){
            nbar[l, h1, h2] = n[l, (h1-1)*H + h2]
          }
        }
      }
      m_saved[iter,,] = m
    }
    return(m_saved)
  }else{
    alpha0 = rep(hyper[1],L); alpha = rep(hyper[2],L)
    k0 = hyper[3]; a_sigma = hyper[4]; b_sigma = hyper[5]
    pb <- progress_bar$new(
      format = " MCMC [:bar] :percent Estimated completion time: :eta",
      total = totiter, clear = FALSE, width= 100)
    #mcmc
    for (iter in 1:totiter){
      pb$tick()
      for (l in 1:L){
        ql_sum = c(rev(cumsum(rev(q[l,])))[2:H0],0)
        for (h1 in 1:H0){
          b0[h1, l] = rbeta(1, 1 + q[l,h1], alpha0[l] + ql_sum[h1])
          if(h1>1){
            pi0[h1, l] = log(b0[h1, l]) + sum(log(1 - b0[1:(h1-1), l]))
          }else{
            pi0[h1, l] = log(b0[h1, l])
          }
          if (l==1) { # covariates
            out = posterior_sampler(X[m[l,]==h1,l], k0, a_sigma, b_sigma)
          } else { # outcome
            out = posterior_sampler_outcome(X[m[2,]==h1, 2],
                                            X[m[2,]==h1,1],
                                            k0, a_sigma, b_sigma)
          }
          theta0[h1, l] = out$mu
          sigma0[h1, l] = out$sigma
          
          nbar_sum = c(rev(cumsum(rev(nbar[l,h1,])))[2:H],0)
          for(h2 in 1:H){
            b[h1, h2, l] = rbeta(1, 1 + nbar[l,h1,h2], alpha[l] + nbar_sum[h2])
            if(h2>1){
              pim[h1, h2, l] = log(b[h1, h2, l]) + sum(log(1 - b[h1, 1:(h2-1),l]))
            }else{
              pim[h1, h2, l] = log(b[h1, h2, l])
            }
          }
        }
      }
      for (l in 1:L){
        for (i in 1:n_tot){
          if(l>1){past = m[l-1, i]}else{past = 1}
          fut = rep(0,H)
          if(l<L){
            for (cc in 1:H){
              for (d in 1:H0){
                fut[cc] = fut[cc] + 
                  (k[l+1,(k[l,cc]-1)*H + d] == m[l+1,i])*exp(pim[k[l,cc],d,l+1])
              } 
            }
            fut = log(fut)
          }
          prob = (pim[past,,l]) + kernel_eval(X[i,l], theta0[ ,l], sigma0[, l]) + fut
          if(max(prob)==-Inf){prob[]=1}
          prob = prob - max(prob)
          prob = exp(prob) 
          if(sum(exp(fut))!=0){ #if fut is zero it cannot move
            c[l, i] = sample(((past-1)*H+1):(past*H), 1, prob = prob)
          }
        }
        for (cc in 1:(H*H0)){
          prob = NULL
          for(kk in 1:H0){
            prob[kk] = pi0[kk,l] + kernel_eval_joint(X[c[l,]==cc,l], theta0[kk,l], sigma0[kk,l])   #check what happens if there is no one in cc
          }
          prob = prob - max(prob)
          prob = exp(prob) 
          k[l, cc] = sample(1:H0, 1, prob = prob)
        }
        for (i in 1:n_tot){
          m[l,i] = k[l,c[l,i]]
        }
      } 
      q[,] = 0
      for(l in 1:L){
        
        temp = rep(0, H*H0)
        temp[unique(c[l,])] = 1
        for(h in 1:H0){
          q[l,h] = sum((k[l,]==h)*temp)
        }
        
        for(cc in 1:(H*H0)){
          n[l,cc] = sum(c[l,]==cc)
        }
        
        for(h1 in 1:H0){
          for(h2 in 1:H){
            nbar[l, h1, h2] = n[l, (h1-1)*H + h2]
          }
        }
        alpha[l] = rgamma(1, 3 + H*H0, rate = 3 - sum(log(1-b[,,l])))
        if(alpha[l]<0.1){alpha[l]=0.1}
        alpha0[l] = rgamma(1, 3 + H0, rate = 3 - sum(log(1-b0[,l])))
        if(alpha0[l]<0.1){alpha0[l]=0.1}
      }
      #m_saved[iter,,] = m
      for (layer in (1:L)){
        write.table(t(m[layer,]), file = paste("layer",layer,".csv",sep=""), append = TRUE, row.names = FALSE,
                    col.names= FALSE)
      }
      #print(rand.index(m[1,],true_layer[,1]))
      #print(rand.index(m[100,],true_layer[,100]))
    }
    return(m_saved)
  }
}