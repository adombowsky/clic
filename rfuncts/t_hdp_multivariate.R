#t-HDP _ conditional algorithm - NNIW model - multivariate
#from: https://github.com/beatricefranzolini/CPE/blob/main/telescopic_HDP_NNIW_multi.R
library(progress)#to draw the progress bar
library(mvtnorm) #for multivariate normal density
library(LaplacesDemon) #for inverse Wishart 

#markovian dependence for L layers, multivariate layers
#normal kernel and normal-inverse-wishart base (correlation within cluster allowed)
#the mean of the base is fixed to 0
#hyper are the two concentration parameters, the scale parameter k0 for the mean
#if alpharandom=TRUE the two first entries of hyper are only the initialization
#data is the dataset nX(p_1+...p_L) where p_l is the dimension of layer l 
#H0 and H are the number of mixture components used in the approximation 
#of the infinite mixtures

#LAYERS CAN HAVE DIFFERENT DIMENSIONS
#WARNING: but MAY NOT WORK IF THERE ARE UNIVARIATE LAYERS #to be tested
telescopic_HDP_NNIW_multi <- function(data, colayer, hyper=c(0.1, 0.1, 0.1),
                                      alpharandom = FALSE,
                                      H0 = 10, H = 10, totiter = 1000){
  #has to work for many X and one param
  #dmvnorm(x, mean = rep(0, p), sigma = diag(p), log = FALSE, checkSymmetry = TRUE)
  kernel_eval_joint = function(x, mean, sigma){
    if(!is.matrix(sigma)){sigma = matrix(sigma)}
    return(sum(dmvnorm(x, mean, sigma, log = TRUE)))
  }
  
  #has to work for one X and many param
  kernel_eval = function(x, mean, sigma){
    temp = NULL
    for (jj in (1:dim(mean)[1])){
      sigma_temp = sigma[jj,,]
      if(!is.matrix(sigma_temp)){sigma_temp = matrix(sigma_temp)}
      temp = c(temp, dmvnorm(x, mean[jj,], sigma_temp, log = TRUE))
    }
    return(temp)
  }
  
  posterior_sampler = function(x, mu0 = NULL, k0 = 0.1, Lambda0 = NULL, v0 = NULL){
    if(!is.matrix(x)){x = matrix(x, ncol=length(x))}
    n_temp = dim(x)[1]; p = dim(x)[2] 
    if(is.null(mu0)){mu0 = rep(0,p)}
    if(is.null(Lambda0)){Lambda0 = diag(p)}
    if(is.null(v0)){v0 = p}
    vn = v0 + n_temp; kn = k0 + n_temp
    if(n_temp == 1){
      S = matrix(0, nrow = p, ncol = p)
      xbar = colMeans(x)
    }else if(n_temp==0){
      S = matrix(0, nrow = p, ncol = p)
      xbar = mu0
    }else{
      S = cov(x) * (n_temp - 1)
      xbar = colMeans(x)
    }
    Lambdan = Lambda0 + S + 
      k0*n_temp / (k0 + n_temp) * (xbar - mu0)%*%t((xbar - mu0))
    mun = (k0 * mu0 + n_temp * xbar) / (kn)
    Sigma_post = MCMCpack::riwish(vn, Lambdan)
    mu_post = rmvnorm(1, mun, kn^(-1)*Sigma_post)
    return(list("mu" = mu_post, 
                "Sigma" = Sigma_post))
  }
  
  
  
  k0 = hyper[3]
  L = length(unique(colayer)) #number of layers
  Layerdim = table(colayer) #number of column for each layer
  n_tot = dim(data)[1] #number of items
  
  #the data in a matrix X, which is nxL
  X = data 
  
  #RANDOM PROB.S
  pim = array(NA, dim = c(H0, H, L)) #weights (in log scale)
  b = array(NA, dim = c(H0, H, L)) #sticks for pi
  
  pi0 = matrix(NA, nrow = H0, ncol = L) #weights of the common overall in HDP
  b0 = matrix(NA, nrow = H0, ncol = L) #sticks for pi0 (in log scale)
  mu0 = array(NA, dim = c(H0, dim(data)[2])) #atoms
  Sigma0 = array(NA, dim = c(H0, dim(data)[2], dim(data)[2]))
  
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
    col_ends = cumsum(Layerdim)[l]
    col_init = col_ends - Layerdim[l] + 1
    m[l, ] = stats::kmeans(as.matrix(X[,col_init:col_ends]),min(10,H0))$cluster
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
    pb <- progress_bar$new(
      format = " MCMC [:bar] :percent Estimated completion time: :eta",
      total = totiter, clear = FALSE, width= 100)
    #mcmc
    for (iter in 1:totiter){
      pb$tick()
      for (l in 1:L){
        col_ends = cumsum(Layerdim)[l]
        col_init = col_ends - Layerdim[l] + 1
        ql_sum = c(rev(cumsum(rev(q[l,])))[2:H0],0)
        for (h1 in 1:H0){
          b0[h1, l] = rbeta(1, 1 + q[l,h1], alpha0 + ql_sum[h1])
          if(h1>1){
            pi0[h1, l] = log(b0[h1, l]) + sum(log(1 - b0[1:(h1-1), l]))
          }else{
            pi0[h1, l] = log(b0[h1, l])
          }
          temp = posterior_sampler(X[m[l,]==h1, col_init:col_ends])
          mu0[h1, col_init:col_ends ] = temp$mu
          Sigma0[h1, col_init:col_ends, col_init:col_ends] = temp$Sigma
          
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
        col_ends = cumsum(Layerdim)[l]
        col_init = col_ends - Layerdim[l] + 1
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
          prob = (pim[past,,l]) + kernel_eval(X[i,col_init:col_ends ], 
                                              mu0[ ,col_init:col_ends], 
                                              Sigma0[,col_init:col_ends, col_init:col_ends]) + fut
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
            prob[kk] = pi0[kk,l] + kernel_eval_joint(X[c[l,]==cc,col_init:col_ends], 
                                                     mu0[kk,col_init:col_ends],
                                                     Sigma0[kk,col_init:col_ends, col_init:col_ends])
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
      # print(rand_index(true_layer[,1],m[1,]))
      # print(rand_index(true_layer[,2],m[2,]))
    }
    return(m_saved)
  }else{
    alpha0 = rep(hyper[1],L); alpha = rep(hyper[2],L)
    pb <- progress_bar$new(
      format = " MCMC [:bar] :percent Estimated completion time: :eta",
      total = totiter, clear = FALSE, width= 100)
    #mcmc
    for (iter in 1:totiter){
      pb$tick()
      for (l in 1:L){
        col_ends = cumsum(Layerdim)[l]
        col_init = col_ends - Layerdim[l] + 1
        ql_sum = c(rev(cumsum(rev(q[l,])))[2:H0],0)
        for (h1 in 1:H0){
          b0[h1, l] = rbeta(1, 1 + q[l,h1], alpha0[l] + ql_sum[h1])
          if(h1>1){
            pi0[h1, l] = log(b0[h1, l]) + sum(log(1 - b0[1:(h1-1), l]))
          }else{
            pi0[h1, l] = log(b0[h1, l])
          }
          temp = posterior_sampler(X[m[l,]==h1, col_init:col_ends])
          mu0[h1, col_init:col_ends ] = temp$mu
          Sigma0[h1, col_init:col_ends, col_init:col_ends] = temp$Sigma
          
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
        col_ends = cumsum(Layerdim)[l]
        col_init = col_ends - Layerdim[l] + 1
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
          prob = (pim[past,,l]) + kernel_eval(X[i,col_init:col_ends ], 
                                              mu0[ ,col_init:col_ends], 
                                              Sigma0[,col_init:col_ends, col_init:col_ends]) + fut
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
            prob[kk] = pi0[kk,l] + kernel_eval_joint(X[c[l,]==cc,col_init:col_ends], 
                                                     mu0[kk,col_init:col_ends],
                                                     Sigma0[kk,col_init:col_ends, col_init:col_ends])
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
      m_saved[iter,,] = m
      # print(rand_index(true_layer[,1],m[1,]))
      # print(rand_index(true_layer[,2],m[2,]))
    }
    return(m_saved)
  }
}