#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double ldmvnorm_arma(arma::vec y, arma::vec mu, arma::mat Sigma, int d){
  double C = -d * 0.5 * log(2*arma::datum::pi);
  double D = - 0.5 * arma::log_det_sympd(Sigma);
  double E = - 0.5 * dot(y - mu, arma::inv_sympd(Sigma) * (y-mu));
  double x = C + D + E;
  return(x);
}

// [[Rcpp::export]]
double ldnorm_arma(double y, double mu, double sigma){
  double C = - 0.5 * log(2*arma::datum::pi);
  double D = - log(sigma);
  double E = - 0.5 * pow((y-mu)/sigma,2);
  double x = C + D + E;
  return(x);
}

// [[Rcpp::export]]
arma::vec rdirichlet_arma(arma::vec a){
  int K = a.n_rows;
  arma::vec Y = arma::zeros(K,1);
  for (int k=0; k<K; k++) {
    arma::vec samp = Rcpp::rgamma(1, a(k), 1);
    Y(k) = as_scalar(samp);
  }
  arma::vec X = Y/sum(Y);
  return(X);
}

//[[Rcpp::export]]
double rbeta_arma(double a, double b) {
  double x = Rcpp::rgamma(1,a,1)(0);
  double y = Rcpp::rgamma(1,b,1)(0);
  double z = x/(x+y);
  return(z);
}

// [[Rcpp::export]]
arma::vec rmvnorm_arma(arma::vec mu, arma::mat Sigma){
  int p = mu.n_rows;
  arma::vec x = mu + chol(Sigma) * arma::randn(p);
  return(x);
}

// [[Rcpp::export]]
double rnorm_arma(double mu, double sigma){
  double x = mu + sigma * as_scalar(arma::randn(1));
  return(x);
}

// [[Rcpp::export]]
arma::mat rinvwishart_arma(int nu, arma::mat Psi){
  int p = Psi.n_rows;
  int nu_inv = nu + p - 1;
  arma::mat Psi_inv = arma::inv_sympd(Psi);
  arma::mat Y = arma::zeros(p, nu_inv);
  for (int i = 0; i<nu_inv; i++){
    Y.col(i) = rmvnorm_arma(arma::zeros(p,1), Psi_inv);
  }
  arma::mat S_Y = Y * Y.t();
  return(inv_sympd(S_Y));
}

// [[Rcpp::export]]
int sample_arma(arma::vec probs) {
  int K = probs.n_rows;
  IntegerVector clusts = Range(1,K);
  IntegerVector samp = RcppArmadillo::sample(clusts, 1, TRUE, probs);
  int s = samp(0);
  return(s);
}

// [[Rcpp::export]]
// for computing the unsigned Stirling numbers of the first kind
// note: rows and columns go from 0:N-1, so best to set N=N+1
arma::mat stirling(int N) {
  arma::mat stir = arma::eye(N,N);
  for (int n = 2; n < N; n++) {
    for (int m = 1; m < n; m++ ) {
      stir(n, m) = stir(n-1, m-1) + (n-1)*stir(n-1, m);
    }
  }
  return(stir);
}

// [[Rcpp::export]]
// from https://github.com/danieledurante/ESBM/blob/master/Source/stirling.cpp
// computes the log-stirling numbers
arma::mat lstirling(int n){
  arma::mat LogS(n+1,n+1);  LogS.fill(-arma::datum::inf);
  
  // Fill the starting values
  LogS(0,0) = 0;
  LogS(1,1) = 0;
  
  for(int i = 2; i <= n; i++){
    for(int j = 1; j < i; j++){
      LogS(i,j) = LogS(i-1,j) + std::log(i-1 + std::exp(LogS(i-1,j-1) - LogS(i-1,j))); 
    }
    LogS(i,i)  = 0;
  }
  return(LogS);
}

// [[Rcpp::export]]
int gibbs_m_kk(int n_kk, double rho, double psi_k, double psi_k_prime) {
  arma::mat ls = lstirling(n_kk);
  arma::vec lprbs = arma::zeros(n_kk+1);
  for (int m = 0; m <= n_kk; m++) {
    lprbs(m) = ls(n_kk, m) + m * (log(rho) + log(psi_k) + log(psi_k_prime));
  }
  arma::vec prbs = arma::exp(lprbs - arma::max(lprbs));
  int m_kk = sample_arma(prbs/arma::sum(prbs)) - 1;
  return(m_kk);
}


// [[Rcpp::export]]
arma::mat gibbs_M(arma::mat N, double rho, arma::vec Psi, arma::vec Psi_prime) {
  int K = Psi.n_rows;
  int K_prime = Psi_prime.n_rows;
  arma::mat M = arma::zeros(K, K_prime);
  for (int k=0; k<K; k++) {
    for (int k_prime=0; k_prime<K_prime; k_prime++) {
      if (N(k,k_prime)==0) {
        continue;
      } else {
        M(k, k_prime) = gibbs_m_kk(N(k,k_prime), rho, Psi(k), Psi_prime(k_prime));
      }
    }
  }
  return(M);
}
// 
// arma::mat metropolis_M(int n, arma::mat N, arma::mat M, double rho, double gamma1, double gamma_2) {
//   int M_sum = arma::accu(M);
//   // step 1, simulate M_sum
//   arma::vec lprobs = arma::zeros(n,1);
//   for (int i = 0; i < n; i++) {
//     lprobs(i) = i * log(rho);
//   }
//   arma::vec probs = arma::exp(lprobs - arma::max(lprobs));
//   int M_sum_star = sample_arma(probs/arma::sum(probs));
//   // step 2, simulate M
// }

// [[Rcpp::export]]
double gibbs_rho(double rho, int n, double M, double a_rho, double b_rho) {
  // first, sample eta
  double eta = rbeta_arma(rho+1.0,n*1.0);
  // next sample from mixture of gammas
  double o_eta = (a_rho + M - 1.0)/( n*1.0 *(b_rho - log(eta)) );
  double pi_eta = o_eta/(1.0 + o_eta);
  int t = Rcpp::rbinom(1, 1, pi_eta)(0);
  if (t==0) {
    //rho = Rcpp::rgamma(1, a_rho+ M - 1.0, 1.0/(b_rho - log(eta)))(0);
    rho = arma::randg(1, arma::distr_param(a_rho+ M - 1.0, 1.0/(b_rho - log(eta))))(0);
  } else {
    //rho = Rcpp::rgamma(1, a_rho + M, 1.0/(b_rho - log(eta)))(0);
    rho = arma::randg(1, arma::distr_param(a_rho + M, 1.0/(b_rho - log(eta))))(0);
  }
  return(rho);
}

// [[Rcpp::export]]
double gibbs_rho_grid(arma::vec rho_grid, int n, arma::mat N, arma::vec q1, arma::vec q2) {
  int W = rho_grid.n_rows;
  int K1 = N.n_rows;
  int K2 = N.n_cols;
  arma::vec lp = arma::zeros(W,1);
  for (int w=0; w<W; w++) {
    lp(w) = lgamma(rho_grid(w))-lgamma(rho_grid(w) + n);
    for (int k1 = 0; k1<K1; k1++) {
      for (int k2 = 0; k2<K2; k2++) {
        lp(w) = lp(w) +
          lgamma(N(k1,k2) + rho_grid(w) * q1(k1) * q2(k2))-
          lgamma(rho_grid(w) * q1(k1) * q2(k2));
      }
    }
  }
  arma::vec p = exp(lp - arma::max(lp));
  int rho_index = sample_arma(p/arma::sum(p)) - 1;
  return(rho_grid(rho_index));
}

// // [[Rcpp::export]]
// double gibbs_rho_grid(arma::vec rho_grid, int n, int M) {
//   int W = rho_grid.n_rows;
//   arma::vec lp = arma::zeros(W,1);
//   for (int w=0; w<W; w++) {
//     lp(w) = lgamma(rho_grid(w))-lgamma(rho_grid(w) + n) + M * log(rho_grid(w));
//   }
//   arma::vec p = exp(lp - arma::max(lp));
//   int rho_index = sample_arma(p/arma::sum(p)) - 1;
//   return(rho_grid(rho_index));
// }


// [[Rcpp::export]]
List gibbs_simple(arma::mat N, int n, 
           double gamma1, double gamma2,
           double a_rho, double b_rho,
           int R) {
  // initialization
  int L1 = N.n_rows;
  int L2 = N.n_cols;
  arma::vec q1 = rdirichlet_arma((gamma1/(L1*1.0)) * arma::ones(L1,1));
  arma::vec q2 = rdirichlet_arma( (gamma2/(L2*1.0)) * arma::ones(L2,1));
  arma::mat M = N;
  double M_sum = arma::accu(M);
  double rho = 1.0;
  arma::vec rho_out = arma::ones(R);
  arma::vec M_sum_out = n * arma::ones(R);
  // sampling
  for (int r=1; r<R; r++) {
    // sample M
    M = gibbs_M(N,rho,q1,q2);
    //M = gibbs_M(N=N, rho=rho, Psi=q1, Psi_prime=q2);
    // sample psi and psi_prime
    arma::vec M_r = arma::sum(M,1);
    q1 = rdirichlet_arma(M_r + (gamma1/(L1*1.0)) * arma::ones(L1,1));
    arma::rowvec M_c = arma::sum(M,0);
    q2 = rdirichlet_arma(M_c.t() + (gamma2/(L2*1.0)) * arma::ones(L2,1));
    // sample rho
    M_sum = arma::accu(M);
    rho = gibbs_rho(rho, n, M_sum, a_rho, b_rho);
    // saving
    rho_out(r) = rho;
    M_sum_out(r) = M_sum;
  }
  return(List::create(Named("rho") = rho_out, Named("M") = M_sum_out));
}

// [[Rcpp::export]]
arma::mat gibbs_C_multiview(arma::vec c1, arma::vec c2, double rho,
                  int L1, int L2,
                  arma::vec theta1, arma::vec theta2,
                  arma::vec X1, arma::vec X2,
                  arma::vec q1, arma::vec q2,
                  double sigma1, double sigma2){
  int n = X1.n_rows;
  int n_k1k2 = n;
  arma::mat lp = arma::zeros(L1,L2);
  arma::mat p = arma::zeros(L1,L2);
  arma::vec p1 = arma::zeros(L1,1);
  for (int i=0; i<n; i++){
    for (int k1=0; k1<L1;k1++) {
      for (int k2=0; k2<L2;k2++) {
        if ((c1(i)==k1+1 & c2(i)==k2+1)) {
          n_k1k2 = arma::sum( (c1==k1+1) % (c2 == k2+1) ) -1;
        } else {
          n_k1k2 = arma::sum( (c1==k1+1) % (c2 == k2+1) );
        }
        lp(k1,k2) = log(rho*q1(k1)*q2(k2) + n_k1k2) + 
          arma::log_normpdf(X1(i), theta1(k1), sigma1) + 
          arma::log_normpdf(X2(i), theta2(k2), sigma2);
      }
    }
    p = arma::exp(lp - arma::max(max(lp,0)))/arma::accu( exp(lp - max(max(lp,0))) );
    p1 = sum(p, 1);
    c1(i) = sample_arma(p1);
    c2(i) = sample_arma( p.row(c1(i)-1).t()/p1(c1(i)-1)  );
  }
  return(arma::join_horiz(c1,c2));
}

// [[Rcpp::export]]
arma::mat gibbs_C_multiview_markov(arma::vec c1, arma::vec c2, double rho,
                            int L1, int L2,
                            arma::vec theta1, arma::vec theta2,
                            arma::vec X1, arma::vec X2,
                            arma::vec q1, arma::vec q2,
                            double sigma1, double sigma2){
  int n = X1.n_rows;
  int n_k1 = n;
  int n_k1k2 = n;
  arma::vec lp1 = arma::zeros(L1,1);
  arma::vec p1 = arma::zeros(L1,1);
  arma::vec lp2 = arma::zeros(L2,1);
  arma::vec p2 = arma::zeros(L2,1);
  for (int i=0; i<n; i++){
    // first sample c1
    for (int k1=0; k1<L1;k1++) {
      if (c1(i)==k1+1) {
        n_k1 = arma::sum(c1==k1+1) - 1;
      } else {
        n_k1 = arma::sum(c1==k1+1);
      }
      lp1(k1) = log(rho * q1(k1) + n_k1) + arma::log_normpdf(X1(i), theta1(k1), sigma1);
    }
    p1 = arma::exp(lp1 - arma::max(lp1));
    c1(i) = sample_arma(p1/arma::sum(p1));
    // next sample c2
    for (int k2=0; k2<L2; k2++) {
      if (c2(i)==k2+1) {
        n_k1k2 =  arma::sum( (c1 == c1(i)) % (c2== k2+1) ) -1;
      } else {
        n_k1k2 =  arma::sum( (c1 == c1(i)) % (c2== k2+1) );
      }
      lp2(k2) = log(rho*q1(c1(i)-1)*q2(k2) + n_k1k2) + arma::log_normpdf(X2(i), theta2(k2), sigma2);
    }
    p2 = arma::exp(lp2 - arma::max(lp2));
    c2(i) = sample_arma(p2/arma::sum(p2));
  }
  return(arma::join_horiz(c1,c2));
}

// [[Rcpp::export]]
arma::mat gibbs_C_multiview_markov_jitter(arma::vec c1, arma::vec c2, double rho,
                                   int L1, int L2,
                                   arma::vec theta1, arma::vec theta2,
                                   arma::vec X1, arma::vec X2,
                                   arma::vec q1, arma::vec q2,
                                   double sigma1, double sigma2){
  int n = X1.n_rows;
  int n_k1 = n;
  int n_k1k2 = n;
  arma::vec lp1 = arma::zeros(L1,1);
  arma::vec p1 = arma::zeros(L1,1);
  arma::vec lp2 = arma::zeros(L2,1);
  arma::vec p2 = arma::zeros(L2,1);
  int b = rbinom(1,1,0.5)(0);
  if (b==1) {
    for (int i=0; i<n; i++){
      // first sample c1
      for (int k1=0; k1<L1;k1++) {
        if (c1(i)==k1+1) {
          n_k1 = arma::sum(c1==k1+1) - 1;
        } else {
          n_k1 = arma::sum(c1==k1+1);
        }
        lp1(k1) = log(rho * q1(k1) + n_k1) + arma::log_normpdf(X1(i), theta1(k1), sigma1);
      }
      p1 = arma::exp(lp1 - arma::max(lp1));
      c1(i) = sample_arma(p1/arma::sum(p1));
      // next sample c2
      for (int k2=0; k2<L2; k2++) {
        if (c2(i)==k2+1) {
          n_k1k2 =  arma::sum( (c1 == c1(i)) % (c2== k2+1) ) -1;
        } else {
          n_k1k2 =  arma::sum( (c1 == c1(i)) % (c2== k2+1) );
        }
        lp2(k2) = log(rho*q1(c1(i)-1)*q2(k2) + n_k1k2) + arma::log_normpdf(X2(i), theta2(k2), sigma2);
      }
      p2 = arma::exp(lp2 - arma::max(lp2));
      c2(i) = sample_arma(p2/arma::sum(p2));
    }
  } else {
    for (int i=0; i<n; i++){
      // first sample c1
      for (int k2=0; k2<L2;k2++) {
        if (c2(i)==k2+1) {
          n_k1 = arma::sum(c2==k2+1) - 1;
        } else {
          n_k1 = arma::sum(c2==k2+1);
        }
        lp1(k2) = log(rho * q2(k2) + n_k1) + arma::log_normpdf(X2(i), theta2(k2), sigma2);
      }
      p1 = arma::exp(lp1 - arma::max(lp1));
      c2(i) = sample_arma(p1/arma::sum(p1));
      // next sample c2
      for (int k1=0; k1<L1; k1++) {
        if (c1(i)==k1+1) {
          n_k1k2 =  arma::sum( (c2 == c2(i)) % (c1== k1+1) ) -1;
        } else {
          n_k1k2 =  arma::sum( (c2 == c2(i)) % (c1== k1+1) );
        }
        lp2(k1) = log(rho*q2(c2(i)-1)*q1(k1) + n_k1k2) + arma::log_normpdf(X1(i), theta1(k1), sigma1);
      }
      p2 = arma::exp(lp2 - arma::max(lp2));
      c1(i) = sample_arma(p2/arma::sum(p2));
    }
  }
  return(arma::join_horiz(c1,c2));
}

// [[Rcpp::export]]
arma::mat gibbs_C_multiview_markov_simple(arma::vec c1, arma::vec c2, double rho,
                                   int L1, int L2,
                                   arma::vec theta1, arma::vec theta2,
                                   arma::vec X1, arma::vec X2,
                                   double gamma1, double gamma2,
                                   double sigma1, double sigma2){
  int n = X1.n_rows;
  int n_k1 = n;
  int n_k1k2 = n;
  arma::vec lp1 = arma::zeros(L1,1);
  arma::vec p1 = arma::zeros(L1,1);
  arma::vec lp2 = arma::zeros(L2,1);
  arma::vec p2 = arma::zeros(L2,1);
  for (int i=0; i<n; i++){
    // first sample c1
    for (int k1=0; k1<L1;k1++) {
      if (c1(i)==k1+1) {
        n_k1 = arma::sum(c1==k1+1) - 1;
      } else {
        n_k1 = arma::sum(c1==k1+1);
      }
      lp1(k1) = log(rho * (gamma1 * gamma2/L1) + n_k1) + arma::log_normpdf(X1(i), theta1(k1), sigma1);
    }
    p1 = arma::exp(lp1 - arma::max(lp1));
    c1(i) = sample_arma(p1/arma::sum(p1));
    // next sample c2
    for (int k2=0; k2<L2; k2++) {
      if (c2(i)==k2+1) {
        n_k1k2 =  arma::sum( (c1 == c1(i)) % (c2== k2+1) ) -1;
      } else {
        n_k1k2 =  arma::sum( (c1 == c1(i)) % (c2== k2+1) );
      }
      lp2(k2) = log(rho*(gamma1/L1)*(gamma2/L2) + n_k1k2) + arma::log_normpdf(X2(i), theta2(k2), sigma2);
    }
    p2 = arma::exp(lp2 - arma::max(lp2));
    c2(i) = sample_arma(p2/arma::sum(p2));
  }
  return(arma::join_horiz(c1,c2));
}

// [[Rcpp::export]]
arma::mat gibbs_C_multiview_markov_fixed_nu(arma::vec c1, arma::vec c2, double rho,
                                   int L1, int L2,
                                   arma::vec theta1, arma::vec theta2,
                                   double gamma1, double gamma2,
                                   arma::vec X1, arma::vec X2,
                                   double sigma1, double sigma2){
  int n = X1.n_rows;
  int n_k1 = n;
  int n_k1k2 = n;
  arma::vec lp1 = arma::zeros(L1,1);
  arma::vec p1 = arma::zeros(L1,1);
  arma::vec lp2 = arma::zeros(L2,1);
  arma::vec p2 = arma::zeros(L2,1);
  for (int i=0; i<n; i++){
    // first sample c1
    for (int k1=0; k1<L1;k1++) {
      if (c1(i)==k1+1) {
        n_k1 = arma::sum(c1==k1+1) - 1;
      } else {
        n_k1 = arma::sum(c1==k1+1);
      }
      lp1(k1) = log(rho * (gamma1 * gamma2/L1) + n_k1) + arma::log_normpdf(X1(i), theta1(k1), sigma1);
    }
    p1 = arma::exp(lp1 - arma::max(lp1));
    c1(i) = sample_arma(p1/arma::sum(p1));
    // next sample c2
    for (int k2=0; k2<L2; k2++) {
      if (c2(i)==k2+1) {
        n_k1k2 =  arma::sum( (c1 == c1(i)) % (c2== k2+1) ) -1;
      } else {
        n_k1k2 =  arma::sum( (c1 == c1(i)) % (c2== k2+1) );
      }
      lp2(k2) = log(rho*(gamma1 * gamma2/(L1*L2)) + n_k1k2) + arma::log_normpdf(X2(i), theta2(k2), sigma2);
    }
    p2 = arma::exp(lp2 - arma::max(lp2));
    c2(i) = sample_arma(p2/arma::sum(p2));
  }
  return(arma::join_horiz(c1,c2));
}

// [[Rcpp::export]]
arma::mat gibbs_C_multiview_markov_fixed_nu_jitter(arma::vec c1, arma::vec c2, double rho,
                                            int L1, int L2,
                                            arma::vec theta1, arma::vec theta2,
                                            double gamma1, double gamma2,
                                            arma::vec X1, arma::vec X2,
                                            double sigma1, double sigma2){
  int n = X1.n_rows;
  int n_k1 = n;
  int n_k1k2 = n;
  arma::vec lp1 = arma::zeros(L1,1);
  arma::vec p1 = arma::zeros(L1,1);
  arma::vec lp2 = arma::zeros(L2,1);
  arma::vec p2 = arma::zeros(L2,1);
  int b = rbinom(1,1,0.5)(0);
  if (b==1) {
    for (int i=0; i<n; i++){
      // first sample c1
      for (int k1=0; k1<L1;k1++) {
        if (c1(i)==k1+1) {
          n_k1 = arma::sum(c1==k1+1) - 1;
        } else {
          n_k1 = arma::sum(c1==k1+1);
        }
        lp1(k1) = log(rho * (gamma1 * gamma2/L1) + n_k1) + arma::log_normpdf(X1(i), theta1(k1), sigma1);
      }
      p1 = arma::exp(lp1 - arma::max(lp1));
      c1(i) = sample_arma(p1/arma::sum(p1));
      // next sample c2
      for (int k2=0; k2<L2; k2++) {
        if (c2(i)==k2+1) {
          n_k1k2 =  arma::sum( (c1 == c1(i)) % (c2== k2+1) ) -1;
        } else {
          n_k1k2 =  arma::sum( (c1 == c1(i)) % (c2== k2+1) );
        }
        lp2(k2) = log(rho*(gamma1 * gamma2/(L1*L2)) + n_k1k2) + arma::log_normpdf(X2(i), theta2(k2), sigma2);
      }
      p2 = arma::exp(lp2 - arma::max(lp2));
      c2(i) = sample_arma(p2/arma::sum(p2));
    }
  } else {
    for (int i=0; i<n; i++){
      // first sample c1
      for (int k2=0; k2<L2;k2++) {
        if (c2(i)==k2+1) {
          n_k1 = arma::sum(c2==k2+1) - 1;
        } else {
          n_k1 = arma::sum(c2==k2+1);
        }
        lp1(k2) = log(rho * (gamma1 * gamma2/L2) + n_k1) + arma::log_normpdf(X2(i), theta2(k2), sigma2);
      }
      p1 = arma::exp(lp1 - arma::max(lp1));
      c2(i) = sample_arma(p1/arma::sum(p1));
      // next sample c2
      for (int k1=0; k1<L1; k1++) {
        if (c1(i)==k1+1) {
          n_k1k2 =  arma::sum( (c2 == c2(i)) % (c1== k1+1) ) -1;
        } else {
          n_k1k2 =  arma::sum( (c2 == c2(i)) % (c1== k1+1) );
        }
        lp2(k1) = log(rho*(gamma1 * gamma2/(L1*L2)) + n_k1k2) + arma::log_normpdf(X1(i), theta1(k1), sigma1);
      }
      p2 = arma::exp(lp2 - arma::max(lp2));
      c1(i) = sample_arma(p2/arma::sum(p2));
    }
  }
  return(arma::join_horiz(c1,c2));
}


// [[Rcpp::export]]
arma::mat arma_table(arma::vec c1, arma::vec c2, int L1, int L2){
  arma::mat N = arma::zeros(L1,L2);
  for (int k1 = 0; k1<L1; k1++) {
    for (int k2 = 0; k2<L2; k2++) {
      N(k1,k2) = sum((c1==k1+1) % (c2==k2+1));
    }
  }
  return(N);
}

// [[Rcpp::export]]
int arma_which_index(arma::vec v, int l) {
  int K = v.n_rows;
  int ind = 0;
  for (int k=0; k<K; k++) {
    if (v(k)==l+1) {
      ind = k;
      break;
    } else {
      continue;
    }
  }
  return(ind);
}

// [[Rcpp::export]]
arma::vec expand(arma::vec v, arma::vec c, int L) {
  arma::vec u = unique(c);
  int K = u.n_rows;
  arma::vec w = arma::zeros(L,1);
  for (int l=0; l<L; l++) {
    if (sum(u==l+1)>0) {
      w(l) = v(arma_which_index(u,l));
    }
    else {
      continue;
    }
  }
  return(w);
}

// [[Rcpp::export]]
arma::vec gibbs_theta_multiview(double mu0, double sigma0,
                                double sigma,
                                arma::vec X,
                                arma::vec c,
                                int L) {
  arma::vec theta = arma::zeros(L,1);
  double mu_l = 0.0;
  double sigma_l = 1.0;
  int n_l = X.n_rows;
  arma::vec X_l = X;
  for (int l=0; l<L; l++) {
    n_l = sum(c == l+1);
    if (n_l==0) {
      theta(l) = mu0 + sigma0*arma::randn();
    } else {
      X_l = X.rows(find(c==l+1));
      sigma_l = pow( (1/pow(sigma0,2.0)) + (n_l/pow(sigma,2.0)), -0.5);
      mu_l = pow(sigma_l,2.0) * ( (mu0/pow(sigma0,2.0)) + (sum(X_l)/pow(sigma,2.0))  );
      theta(l) = mu_l + sigma_l * arma::randn();
    }
  }
  //for (int l=0; l<L2; l++) {
  //  n_l = sum(c2 == l+1);
  //  X_2l = X2.rows(c2==l+1);
  //  sigma_l = pow( (1/pow(sigma02,2.0)) + (n_l/pow(sigma2,2.0)), 0.5);
  //  mu_l = pow(sigma_l,2.0) * ( (mu02/pow(sigma02,2.0)) + (sum(X_2l)/pow(sigma2,2.0))  );
  ///  theta2(l) = mu_l + sigma_l * arma::randn();
  //}
  return(theta);
}

// [[Rcpp::export]]
arma::mat gibbs_theta_multiview_bivariate(arma::vec mu0, 
                                arma::mat sigma0,
                                arma::mat sigma,
                                arma::mat X,
                                arma::vec c,
                                int L) {
  int p = X.n_cols;
  arma::mat theta = arma::zeros(L,p);
  arma::mat mu_l = arma::zeros(L,p);
  arma::mat sigma_l = arma::eye(p,p);
  arma::mat sigma0_inv = arma::inv_sympd(sigma0);
  arma::mat sigma_inv = arma::inv_sympd(sigma);
  int n_l = X.n_rows;
  arma::mat X_l = X;
  for (int l=0; l<L; l++) {
    n_l = sum(c == l+1);
    if (n_l==0) {
      theta.row(l) = (mu0 + chol(sigma0)*arma::randn(p,1)).t();
    } else {
      X_l = X.rows(find(c==l+1));
      sigma_l = arma::inv_sympd( sigma0_inv + n_l *  sigma_inv );
      mu_l = sigma_l * ( sigma0_inv * mu0 + n_l * sigma_inv * arma::mean(X_l,0).t() );
      theta.row(l) = (mu_l + chol(sigma_l) * arma::randn(p,1)).t();
    }
  }
  return(theta);
}

// [[Rcpp::export]]
double gibbs_sigma_multiview(double alpha, double beta,
                             arma::vec theta,
                             arma::vec X,
                             arma::vec c,
                             int L) {
  int n = X.n_rows;
  arma::vec ssqvec = arma::zeros(L,1);
  int n_l = n;
  for (int l=0; l<L; l++) {
    n_l = sum(c==l+1);
    ssqvec(l) = arma::as_scalar(arma::sum(square( X.rows(find(c == l+1)) - theta(l) * arma::ones(n_l,1)) ));
  }
  double ssq = arma::sum(ssqvec);
  double alpha_star = alpha + 0.5*n;
  double beta_star = beta + 0.5 * ssq;
  return( pow(rgamma(1, alpha_star, 1.0/beta_star)(0), -0.5));
}

// [[Rcpp::export]]
double gibbs_sigma_multiview_bivariate(double alpha, double beta,
                             arma::mat theta,
                             arma::mat X,
                             arma::vec c,
                             int L) {
  int n = X.n_rows;
  arma::vec X0 = X.col(0);
  arma::vec X1 = X.col(1);
  arma::vec ssqvec = arma::zeros(L,1);
  int n_l = n;
  for (int l=0; l<L; l++) {
    n_l = sum(c==l+1);
    ssqvec(l) = arma::as_scalar(arma::sum(square( X0.rows(find(c == l+1)) - theta(l,0) * arma::ones(n_l,1)) )) +
      arma::as_scalar(arma::sum(square( X1.rows(find(c == l+1)) - theta(l,1) * arma::ones(n_l,1)) ));
  }
  double ssq = arma::sum(ssqvec);
  double alpha_star = alpha + 0.5*n;
  double beta_star = beta + 0.5 * ssq;
  return( pow(rgamma(1, alpha_star, 1.0/beta_star)(0), -0.5));
}

// [[Rcpp::export]]
arma::mat gibbs_theta_sigma_multiview(double mu0, double kappa0,
                                double alpha, double beta,
                                arma::vec X,
                                arma::vec c,
                                int L) {
  arma::vec theta = arma::zeros(L,1);
  arma::vec sigma = arma::zeros(L,1);
  double mu0_l = 0.0;
  double kappa0_l = 0.0;
  double alpha_l = 1.0;
  double beta_l = 1.0;
  double sigma_l = 1.0;
  int n_l = X.n_rows;
  arma::vec X_l = X;
  for (int l=0; l<L; l++) {
    n_l = sum(c == l+1);
    if (n_l==0) {
      sigma(l) = pow(arma::randg(arma::distr_param(alpha, pow(beta, -1.0))), -0.5);
      theta(l) = mu0 + sigma(l)*pow(kappa0, -0.5) * arma::randn();
    } else {
      X_l = X.rows(find(c==l+1));
      mu0_l = ((kappa0 * mu0) + arma::sum(X_l))/(kappa0 + n_l);
      kappa0_l = kappa0 + n_l;
      alpha_l = alpha + (n_l/2.0);
      beta_l = beta + 0.5 * arma::sum(arma::square(X_l - arma::mean(X_l))) +
        ((n_l * kappa0)/(n_l + kappa0)) * 0.5 * pow(arma::mean(X_l) - mu0,2.0);
      sigma(l) = pow(arma::randg(arma::distr_param(alpha_l, pow(beta_l, -1.0))), -0.5);
      theta(l) = mu0_l + sigma(l) * pow(kappa0_l, -0.5) * arma::randn();
    }
  }
  return(arma::join_horiz(theta,sigma));
}

// [[Rcpp::export]]
List gibbs(int n, arma::mat X,
                  double gamma1, double gamma2,
                  double mu01, double sigma01,
                  double mu02, double sigma02,
                  double alpha1, double beta1,
                  double alpha2, double beta2,
                  int L1, int L2,
                  double a_rho, double b_rho,
                  int R) {
  // initialization
  arma::vec q1 = rdirichlet_arma((gamma1/(L1*1.0)) * arma::ones(L1,1));
  arma::vec q2 = rdirichlet_arma( (gamma2/(L2*1.0)) * arma::ones(L2,1));
  arma::vec c1 = arma::zeros(n,1);
  arma::vec c2 = arma::zeros(n,1);
  for (int i=0; i<n; i++) {
    c1(i) = sample_arma(q1);
    c2(i) = sample_arma(q2);
  }
  arma::mat c1_out = arma::zeros(R,n);
  c1_out.row(0) = c1.t();
  arma::mat c2_out = arma::zeros(R,n);
  c2_out.row(0) = c2.t();
  arma::mat N = arma_table(c1,c2,L1,L2);
  arma::mat M = N;
  double M_sum = arma::accu(M);
  double rho = 1.0;
  arma::vec rho_out = arma::ones(R);
  arma::mat Cmat = arma::zeros(n,2);
  arma::vec theta1 = arma::zeros(L1,1);
  arma::vec theta2 = arma::zeros(L2,1);
  double sigma1 = 1.0;
  double sigma2 = 1.0;
  // sampling
  for (int r=1; r<R; r++) {
    // sample c
    Cmat = gibbs_C_multiview(c1,c2,rho,
                             L1,L2,theta1,theta2,
                             X.col(0),X.col(1),
                             q1,q2,
                             sigma1,sigma2);
    c1 = Cmat.col(0);
    c2 = Cmat.col(1);
    N = arma_table(c1,c2,L1,L2);
    // sample M
    M = gibbs_M(N,rho,q1,q2);
    // sample psi and psi_prime
    arma::vec M_r = arma::sum(M,1);
    q1 = rdirichlet_arma(M_r + (gamma1/(L1*1.0)) * arma::ones(L1,1));
    arma::rowvec M_c = arma::sum(M,0);
    q2 = rdirichlet_arma(M_c.t() + (gamma2/(L2*1.0)) * arma::ones(L2,1));
    // sample rho
    M_sum = arma::accu(M);
    rho = gibbs_rho(rho, n, M_sum, a_rho, b_rho);
    // sample theta
    theta1 = gibbs_theta_multiview(mu01, sigma01,
                                   sigma1,
                                   X.col(0),
                                   c1,
                                   L1);
    theta2 = gibbs_theta_multiview(mu02, sigma02,
                                   sigma2,
                                   X.col(1),
                                   c2,
                                   L2);
    // sample sigma
    sigma1 = gibbs_sigma_multiview(alpha1, beta1,
                                  theta1,
                                   X.col(0),
                                   c1,
                                   L1);
    sigma2 = gibbs_sigma_multiview(alpha2, beta2,
                                   theta2,
                                   X.col(1),
                                   c2,
                                   L2);
    // saving
    rho_out(r) = rho;
    c1_out.row(r) = c1.t();
    c2_out.row(r) = c2.t();
  }
  return(List::create(Named("rho") = rho_out,
                      Named("c1") = c1_out,
                      Named("c2") = c2_out));
}

// [[Rcpp::export]]
List gibbs_markov(int n, arma::mat X,
           double gamma1, double gamma2,
           double mu01, double sigma01,
           double mu02, double sigma02,
           double alpha1, double beta1,
           double alpha2, double beta2,
           int L1, int L2,
           double a_rho, double b_rho,
           int R) {
  // initialization
  arma::vec q1 = rdirichlet_arma((gamma1/(L1*1.0)) * arma::ones(L1,1));
  arma::vec q2 = rdirichlet_arma( (gamma2/(L2*1.0)) * arma::ones(L2,1));
  arma::vec c1 = arma::zeros(n,1);
  arma::vec c2 = arma::zeros(n,1);
  for (int i=0; i<n; i++) {
    c1(i) = sample_arma(q1);
    c2(i) = sample_arma(q2);
  }
  arma::mat c1_out = arma::zeros(R,n);
  c1_out.row(0) = c1.t();
  arma::mat c2_out = arma::zeros(R,n);
  c2_out.row(0) = c2.t();
  arma::mat N = arma_table(c1,c2,L1,L2);
  arma::mat M = N;
  double M_sum = arma::accu(M);
  double rho = 1.0;
  arma::vec rho_out = arma::ones(R);
  arma::mat Cmat = arma::zeros(n,2);
  arma::vec theta1 = arma::zeros(L1,1);
  arma::vec theta2 = arma::zeros(L2,1);
  double sigma1 = 1.0;
  double sigma2 = 1.0;
  // sampling
  for (int r=1; r<R; r++) {
    // sample c
    Cmat = gibbs_C_multiview_markov(c1,c2,rho,
                             L1,L2,theta1,theta2,
                             X.col(0),X.col(1),
                             q1,q2,
                             sigma1,sigma2);
    c1 = Cmat.col(0);
    c2 = Cmat.col(1);
    N = arma_table(c1,c2,L1,L2);
    // sample M
    M = gibbs_M(N,rho,q1,q2);
    // sample psi and psi_prime
    arma::vec M_r = arma::sum(M,1);
    q1 = rdirichlet_arma(M_r + (gamma1/(L1*1.0)) * arma::ones(L1,1));
    arma::rowvec M_c = arma::sum(M,0);
    q2 = rdirichlet_arma(M_c.t() + (gamma2/(L2*1.0)) * arma::ones(L2,1));
    // sample rho
    M_sum = arma::accu(M);
    rho = gibbs_rho(rho, n, M_sum, a_rho, b_rho);
    // sample theta
    theta1 = gibbs_theta_multiview(mu01, sigma01,
                                   sigma1,
                                   X.col(0),
                                   c1,
                                   L1);
    theta2 = gibbs_theta_multiview(mu02, sigma02,
                                   sigma2,
                                   X.col(1),
                                   c2,
                                   L2);
    // sample sigma
    sigma1 = gibbs_sigma_multiview(alpha1, beta1,
                                   theta1,
                                   X.col(0),
                                   c1,
                                   L1);
    sigma2 = gibbs_sigma_multiview(alpha2, beta2,
                                   theta2,
                                   X.col(1),
                                   c2,
                                   L2);
    // saving
    rho_out(r) = rho;
    c1_out.row(r) = c1.t();
    c2_out.row(r) = c2.t();
  }
  return(List::create(Named("rho") = rho_out,
                      Named("c1") = c1_out,
                      Named("c2") = c2_out));
}

// [[Rcpp::export]]
List grid_gibbs_markov(int n, arma::mat X,
                  double gamma1, double gamma2,
                  double mu01, double sigma01,
                  double mu02, double sigma02,
                  double alpha1, double beta1,
                  double alpha2, double beta2,
                  int L1, int L2,
                  arma::vec rho_grid,
                  int R) {
  // initialization
  arma::vec q1 = rdirichlet_arma((gamma1/(L1*1.0)) * arma::ones(L1,1));
  arma::vec q2 = rdirichlet_arma( (gamma2/(L2*1.0)) * arma::ones(L2,1));
  arma::vec c1 = arma::zeros(n,1);
  arma::vec c2 = arma::zeros(n,1);
  for (int i=0; i<n; i++) {
    c1(i) = sample_arma(q1);
    c2(i) = sample_arma(q2);
  }
  arma::mat c1_out = arma::zeros(R,n);
  c1_out.row(0) = c1.t();
  arma::mat c2_out = arma::zeros(R,n);
  c2_out.row(0) = c2.t();
  arma::mat N = arma_table(c1,c2,L1,L2);
  arma::mat M = N;
  double M_sum = arma::accu(M);
  int W = rho_grid.n_cols;
  int rho_index = sample_arma((1/W) * arma::ones(W,1))-1;
  double rho = rho_grid(rho_index);
  arma::vec rho_out = arma::ones(R);
  arma::mat Cmat = arma::zeros(n,2);
  arma::vec theta1 = arma::zeros(L1,1);
  arma::vec theta2 = arma::zeros(L2,1);
  double sigma1 = 1.0;
  double sigma2 = 1.0;
  // sampling
  for (int r=1; r<R; r++) {
    // sample c
    Cmat = gibbs_C_multiview_markov(c1,c2,rho,
                                    L1,L2,theta1,theta2,
                                    X.col(0),X.col(1),
                                    q1,q2,
                                    sigma1,sigma2);
    c1 = Cmat.col(0);
    c2 = Cmat.col(1);
    N = arma_table(c1,c2,L1,L2);
    // sample M
    M = gibbs_M(N,rho,q1,q2);
    // sample psi and psi_prime
    arma::vec M_r = arma::sum(M,1);
    q1 = rdirichlet_arma(M_r + (gamma1/(L1*1.0)) * arma::ones(L1,1));
    arma::rowvec M_c = arma::sum(M,0);
    q2 = rdirichlet_arma(M_c.t() + (gamma2/(L2*1.0)) * arma::ones(L2,1));
    // sample rho
    M_sum = arma::accu(M);
    rho = gibbs_rho_grid(rho_grid, n, N, q1, q2);
    // sample theta
    theta1 = gibbs_theta_multiview(mu01, sigma01,
                                   sigma1,
                                   X.col(0),
                                   c1,
                                   L1);
    theta2 = gibbs_theta_multiview(mu02, sigma02,
                                   sigma2,
                                   X.col(1),
                                   c2,
                                   L2);
    // sample sigma
    sigma1 = gibbs_sigma_multiview(alpha1, beta1,
                                   theta1,
                                   X.col(0),
                                   c1,
                                   L1);
    sigma2 = gibbs_sigma_multiview(alpha2, beta2,
                                   theta2,
                                   X.col(1),
                                   c2,
                                   L2);
    // saving
    rho_out(r) = rho;
    c1_out.row(r) = c1.t();
    c2_out.row(r) = c2.t();
  }
  return(List::create(Named("rho") = rho_out,
                      Named("c1") = c1_out,
                      Named("c2") = c2_out));
}

// [[Rcpp::export]]
List grid_gibbs_longnecker(int n, arma::mat X,
                        arma::vec c1, arma::vec c2,
                       double gamma1, double gamma2,
                       double mu01, double sigma01,
                       double mu02, double sigma02,
                       double alpha1, double beta1,
                       double alpha2, double beta2,
                       int L1, int L2,
                       arma::vec rho_grid,
                       int R) {
  // initialization
  arma::vec q1 = rdirichlet_arma((gamma1/(L1*1.0)) * arma::ones(L1,1));
  arma::vec q2 = rdirichlet_arma( (gamma2/(L2*1.0)) * arma::ones(L2,1));
  arma::mat c1_out = arma::zeros(R,n);
  c1_out.row(0) = c1.t();
  arma::mat c2_out = arma::zeros(R,n);
  c2_out.row(0) = c2.t();
  arma::mat N = arma_table(c1,c2,L1,L2);
  arma::mat M = N;
  double M_sum = arma::accu(M);
  int W = rho_grid.n_cols;
  int rho_index = sample_arma((1/W) * arma::ones(W,1))-1;
  double rho = rho_grid(rho_index);
  arma::vec rho_out = arma::ones(R);
  arma::mat Cmat = arma::zeros(n,2);
  arma::vec theta1 = arma::zeros(L1,1);
  arma::vec theta2 = arma::zeros(L2,1);
  double sigma1 = 1.0;
  double sigma2 = 1.0;
  // sampling
  for (int r=1; r<R; r++) {
    // sample M
    M = gibbs_M(N,rho,q1,q2);
    // sample psi and psi_prime
    arma::vec M_r = arma::sum(M,1);
    q1 = rdirichlet_arma(M_r + (gamma1/(L1*1.0)) * arma::ones(L1,1));
    arma::rowvec M_c = arma::sum(M,0);
    q2 = rdirichlet_arma(M_c.t() + (gamma2/(L2*1.0)) * arma::ones(L2,1));
    // sample rho
    M_sum = arma::accu(M);
    rho = gibbs_rho_grid(rho_grid, n, N, q1, q2);
    // sample theta
    theta1 = gibbs_theta_multiview(mu01, sigma01,
                                   sigma1,
                                   X.col(0),
                                   c1,
                                   L1);
    theta2 = gibbs_theta_multiview(mu02, sigma02,
                                   sigma2,
                                   X.col(1),
                                   c2,
                                   L2);
    // sample sigma
    sigma1 = gibbs_sigma_multiview(alpha1, beta1,
                                   theta1,
                                   X.col(0),
                                   c1,
                                   L1);
    sigma2 = gibbs_sigma_multiview(alpha2, beta2,
                                   theta2,
                                   X.col(1),
                                   c2,
                                   L2);
    // sample c
    Cmat = gibbs_C_multiview_markov_jitter(c1,c2,rho,
                                    L1,L2,theta1,theta2,
                                    X.col(0),X.col(1),
                                    q1,q2,
                                    sigma1,sigma2);
    c1 = Cmat.col(0);
    c2 = Cmat.col(1);
    N = arma_table(c1,c2,L1,L2);
    // saving
    rho_out(r) = rho;
    c1_out.row(r) = c1.t();
    c2_out.row(r) = c2.t();
  }
  return(List::create(Named("rho") = rho_out,
                      Named("c1") = c1_out,
                      Named("c2") = c2_out));
}



// [[Rcpp::export]]
List gibbs_markov_fixed_nu(int n, arma::mat X,
                  double gamma1, double gamma2,
                  double mu01, double sigma01,
                  double mu02, double sigma02,
                  double alpha1, double beta1,
                  double alpha2, double beta2,
                  int L1, int L2,
                  double nu,
                  int R) {
  // initialization
  arma::vec q1 = rdirichlet_arma((gamma1/(L1*1.0)) * arma::ones(L1,1));
  arma::vec q2 = rdirichlet_arma( (gamma2/(L2*1.0)) * arma::ones(L2,1));
  arma::vec c1 = arma::zeros(n,1);
  arma::vec c2 = arma::zeros(n,1);
  for (int i=0; i<n; i++) {
    c1(i) = sample_arma(q1);
    c2(i) = sample_arma(q2);
  }
  arma::mat c1_out = arma::zeros(R,n);
  c1_out.row(0) = c1.t();
  arma::mat c2_out = arma::zeros(R,n);
  c2_out.row(0) = c2.t();
  double rho = (1.0-nu)/(nu);
  arma::mat N = arma_table(c1,c2,L1,L2);
  arma::mat Cmat = arma::zeros(n,2);
  arma::vec theta1 = arma::zeros(L1,1);
  arma::vec theta2 = arma::zeros(L2,1);
  double sigma1 = 1.0;
  double sigma2 = 1.0;
  // sampling
  for (int r=1; r<R; r++) {
    // sample c
    Cmat = gibbs_C_multiview_markov_fixed_nu_jitter(c1,c2,rho,
                                    L1,L2,theta1,theta2,
                                    gamma1, gamma2,
                                    X.col(0),X.col(1),
                                    sigma1,sigma2);
    c1 = Cmat.col(0);
    c2 = Cmat.col(1);
    // sample theta
    theta1 = gibbs_theta_multiview(mu01, sigma01,
                                   sigma1,
                                   X.col(0),
                                   c1,
                                   L1);
    theta2 = gibbs_theta_multiview(mu02, sigma02,
                                   sigma2,
                                   X.col(1),
                                   c2,
                                   L2);
    // sample sigma
    sigma1 = gibbs_sigma_multiview(alpha1, beta1,
                                   theta1,
                                   X.col(0),
                                   c1,
                                   L1);
    sigma2 = gibbs_sigma_multiview(alpha2, beta2,
                                   theta2,
                                   X.col(1),
                                   c2,
                                   L2);
    // saving
    c1_out.row(r) = c1.t();
    c2_out.row(r) = c2.t();
  }
  return(List::create(Named("c1") = c1_out,
                      Named("c2") = c2_out));
}


// [[Rcpp::export]]
List gibbs_fixed_nu(int n, arma::mat X,
           double gamma1, double gamma2,
           double mu01, double sigma01,
           double mu02, double sigma02,
           double alpha1, double beta1,
           double alpha2, double beta2,
           int L1, int L2,
           double nu,
           int R) {
  // initialization
  arma::vec q1 = rdirichlet_arma((gamma1/(L1*1.0)) * arma::ones(L1,1));
  arma::vec q2 = rdirichlet_arma( (gamma2/(L2*1.0)) * arma::ones(L2,1));
  arma::vec c1 = arma::zeros(n,1);
  arma::vec c2 = arma::zeros(n,1);
  for (int i=0; i<n; i++) {
    c1(i) = sample_arma(q1);
    c2(i) = sample_arma(q2);
  }
  arma::mat c1_out = arma::zeros(R,n);
  c1_out.row(0) = c1.t();
  arma::mat c2_out = arma::zeros(R,n);
  c2_out.row(0) = c2.t();
  arma::mat N = arma_table(c1,c2,L1,L2);
  arma::mat M = N;
  double rho = (1.0-nu)/(nu);
  arma::vec rho_out = arma::ones(R);
  arma::mat Cmat = arma::zeros(n,2);
  arma::vec theta1 = arma::zeros(L1,1);
  arma::vec theta2 = arma::zeros(L2,1);
  double sigma1 = 1.0;
  double sigma2 = 1.0;
  // sampling
  for (int r=1; r<R; r++) {
    // sample c
    Cmat = gibbs_C_multiview(c1,c2,rho,
                             L1,L2,theta1,theta2,
                             X.col(0),X.col(1),
                             q1,q2,
                             sigma1,sigma2);
    c1 = Cmat.col(0);
    c2 = Cmat.col(1);
    N = arma_table(c1,c2,L1,L2);
    // sample M
    M = gibbs_M(N,rho,q1,q2);
    // sample psi and psi_prime
    arma::vec M_r = arma::sum(M,1);
    q1 = rdirichlet_arma(M_r + (gamma1/(L1*1.0)) * arma::ones(L1,1));
    arma::rowvec M_c = arma::sum(M,0);
    q2 = rdirichlet_arma(M_c.t() + (gamma2/(L2*1.0)) * arma::ones(L2,1));
    // sample theta
    theta1 = gibbs_theta_multiview(mu01, sigma01,
                                   sigma1,
                                   X.col(0),
                                   c1,
                                   L1);
    theta2 = gibbs_theta_multiview(mu02, sigma02,
                                   sigma2,
                                   X.col(1),
                                   c2,
                                   L2);
    // sample sigma
    sigma1 = gibbs_sigma_multiview(alpha1, beta1,
                                   theta1,
                                   X.col(0),
                                   c1,
                                   L1);
    sigma2 = gibbs_sigma_multiview(alpha2, beta2,
                                   theta2,
                                   X.col(1),
                                   c2,
                                   L2);
    // saving
    rho_out(r) = rho;
    c1_out.row(r) = c1.t();
    c2_out.row(r) = c2.t();
  }
  return(List::create(Named("c1") = c1_out,
                      Named("c2") = c2_out));
}


// [[Rcpp::export]]
arma::vec gibbs_C_vanilla(arma::vec c,
                          int L,
                          arma::vec theta,
                          arma::vec X,
                          double gam,
                          double sigma){
  int n = X.n_rows;
  int n_k = n;
  arma::vec lp = arma::zeros(L,1);
  arma::vec p = arma::zeros(L,1);
  for (int i=0; i<n; i++) {
    for (int k=0; k<L; k++) {
      if (c(i)==k+1) {
        n_k = arma::sum(c==k+1)-1;
      } else {
        n_k = arma::sum(c==k+1);
      }
      lp(k) = log(n_k + (gam/(L*1.0))) + arma::log_normpdf(X(i), theta(k), sigma);
    }
    p = arma::exp(lp - arma::max(lp));
    c(i) = sample_arma(p/sum(p));
  }
  return(c);
}

// [[Rcpp::export]]
arma::mat gibbs_vanilla(int n, arma::vec X,
           double gam,
           double mu0, double sigma0,
           double alpha, double beta,
           int L,
           int R) {
  // initialization
  arma::vec q = rdirichlet_arma((gam/(L*1.0)) * arma::ones(L,1));
  arma::vec c = arma::zeros(n,1);
  for (int i=0; i<n; i++) {
    c(i) = sample_arma(q);
  }
  arma::mat c_out = arma::zeros(R,n);
  c_out.row(0) = c.t();
  arma::vec theta = arma::zeros(L,1);
  double sigma = 1.0;
  // sampling
  for (int r=1; r<R; r++) {
    // sample c
    c = gibbs_C_vanilla(c,
                        L,
                        theta,
                        X,
                        gam,
                        sigma);
    // sample theta
    theta = gibbs_theta_multiview(mu0, 
                                  sigma0,
                                  sigma,
                                  X,
                                  c,
                                  L);
    // sample sigma
    sigma = gibbs_sigma_multiview(alpha, 
                                  beta,
                                  theta,
                                  X,
                                  c,
                                  L);
    // saving
    c_out.row(r) = c.t();
  }
  return(c_out);
}

// [[Rcpp::export]]
arma::mat gibbs_C_conditional(arma::vec c1, arma::vec c2, double rho,
                            int L1, int L2,
                            arma::vec theta1, arma::vec theta2,
                            arma::vec X1, arma::vec X2,
                            arma::vec q1, arma::vec q2,
                            double sigma1, double sigma2){
  int n = X1.n_rows;
  int n_k1k2 = n;
  arma::mat lp = arma::zeros(L1,L2);
  arma::mat p = arma::zeros(L1,L2);
  arma::vec p1 = arma::zeros(L1,1);
  for (int i=0; i<n; i++){
    for (int k1=0; k1<L1;k1++) {
      for (int k2=0; k2<L2;k2++) {
        if ((c1(i)==k1+1 & c2(i)==k2+1)) {
          n_k1k2 = arma::sum( (c1==k1+1) % (c2 == k2+1) ) -1;
        } else {
          n_k1k2 = arma::sum( (c1==k1+1) % (c2 == k2+1) );
        }
        lp(k1,k2) = log(rho*q1(k1)*q2(k2) + n_k1k2) + 
          arma::log_normpdf(X1(i), theta1(k1), sigma1) + 
          arma::log_normpdf(X2(i), X1(i)*theta2(k2), sigma2);
      }
    }
    p = arma::exp(lp - arma::max(max(lp,0)))/arma::accu( exp(lp - max(max(lp,0))) );
    p1 = sum(p, 1);
    c1(i) = sample_arma(p1);
    c2(i) = sample_arma( p.row(c1(i)-1).t()/p1(c1(i)-1)  );
  }
  return(arma::join_horiz(c1,c2));
}

// [[Rcpp::export]]
arma::mat gibbs_C_conditional_markov(arma::vec c1, arma::vec c2, double rho,
                              int L1, int L2,
                              arma::vec theta1, arma::vec theta2,
                              arma::vec X1, arma::vec X2,
                              arma::vec q1, arma::vec q2,
                              double sigma1, double sigma2){
  int n = X1.n_rows;
  int n_k1 = n;
  int n_k1k2 = n;
  arma::vec lp1 = arma::zeros(L1,1);
  arma::vec p1 = arma::zeros(L1,1);
  arma::vec lp2 = arma::zeros(L2,1);
  arma::vec p2 = arma::zeros(L2,1);
  for (int i=0; i<n; i++){
    // first sample c1
    for (int k1=0; k1<L1;k1++) {
      if (c1(i)==k1+1) {
        n_k1 = arma::sum(c1==k1+1) - 1;
      } else {
        n_k1 = arma::sum(c1==k1+1);
      }
      lp1(k1) = log(rho * q1(k1) + n_k1) + arma::log_normpdf(X1(i), theta1(k1), sigma1);
    }
    p1 = arma::exp(lp1 - arma::max(lp1));
    c1(i) = sample_arma(p1/arma::sum(p1));
    // next sample c2
    for (int k2=0; k2<L2; k2++) {
      if (c2(i)==k2+1) {
        n_k1k2 =  arma::sum( (c1 == c1(i)) % (c2== k2+1) ) -1;
      } else {
        n_k1k2 =  arma::sum( (c1 == c1(i)) % (c2== k2+1) );
      }
      lp2(k2) = log(rho*q1(c1(i)-1)*q2(k2) + n_k1k2) + arma::log_normpdf(X2(i), X1(i)*theta2(k2), sigma2);
    }
    p2 = arma::exp(lp2 - arma::max(lp2));
    c2(i) = sample_arma(p2/arma::sum(p2));
  }
  return(arma::join_horiz(c1,c2));
}

// [[Rcpp::export]]
arma::mat gibbs_C_conditional_markov_ls(arma::vec c1, arma::vec c2, double rho,
                                     int L1, int L2,
                                     arma::mat theta1, arma::mat theta2,
                                     arma::vec mu2,
                                     arma::mat X1, arma::vec X2,
                                     arma::vec q1, arma::vec q2,
                                     arma::vec sigma1, arma::vec sigma2){
  int n = X1.n_rows;
  int p = X2.n_cols;
  int n_k1 = n;
  int n_k1k2 = n;
  arma::vec lp1 = arma::zeros(L1,1);
  arma::vec p1 = arma::zeros(L1,1);
  arma::vec lp2 = arma::zeros(L2,1);
  arma::vec p2 = arma::zeros(L2,1);
  for (int i=0; i<n; i++){
    // first sample c1
    for (int k1=0; k1<L1;k1++) {
      if (c1(i)==k1+1) {
        n_k1 = arma::sum(c1==k1+1) - 1;
      } else {
        n_k1 = arma::sum(c1==k1+1);
      }
      lp1(k1) = log(rho * q1(k1) + n_k1) + 
        arma::log_normpdf(X1(i,0), theta1(k1,0), sigma1(0)) + 
        arma::log_normpdf(X1(i,1), theta1(k1,1), sigma1(1));
    }
    p1 = arma::exp(lp1 - arma::max(lp1));
    c1(i) = sample_arma(p1/arma::sum(p1));
    // next sample c2
    for (int k2=0; k2<L2; k2++) {
      if (c2(i)==k2+1) {
        n_k1k2 =  arma::sum( (c1 == c1(i)) % (c2== k2+1) ) -1;
      } else {
        n_k1k2 =  arma::sum( (c1 == c1(i)) % (c2== k2+1) );
      }
      lp2(k2) = log(rho*q1(c1(i)-1)*q2(k2) + n_k1k2) + arma::log_normpdf(X2(i), mu2(k2) + dot(X1.row(i).t(), theta2.row(k2).t()), sigma2(k2));
    }
    p2 = arma::exp(lp2 - arma::max(lp2));
    c2(i) = sample_arma(p2/arma::sum(p2));
  }
  return(arma::join_horiz(c1,c2));
}

// [[Rcpp::export]]
arma::vec gibbs_theta_conditional(double mu0, double sigma0,
                                double sigma,
                                arma::vec X_out,
                                arma::vec X_covars,
                                arma::vec c,
                                int L) {
  arma::vec theta = arma::zeros(L,1);
  double mu_l = 0.0;
  double sigma_l = 1.0;
  int n_l = X_out.n_rows;
  arma::vec X_out_l = X_out;
  arma::vec X_covars_l = X_covars;
  for (int l=0; l<L; l++) {
    n_l = sum(c == l+1);
    if (n_l==0) {
      theta(l) = mu0 + sigma0*arma::randn();
    } else {
      X_out_l = X_out.rows(find(c==l+1));
      X_covars_l = X_covars.rows(find(c==l+1));
      sigma_l = pow( dot(X_covars_l.t(), X_covars_l.t()) + pow(sigma0,-2.0), -0.5);
      mu_l = pow(sigma_l,2.0) * ( dot(X_covars_l.t(), X_out_l.t()) + mu0*pow(sigma0,-2.0) );
      theta(l) = mu_l + sigma * sigma_l * arma::randn();
    }
  }
  return(theta);
}

// [[Rcpp::export]]
double gibbs_sigma_conditional(double alpha, double beta,
                             arma::vec theta,
                             arma::vec X_out,
                             arma::vec X_covars,
                             arma::vec c,
                             int L) {
  int n = X_out.n_rows;
  arma::vec ssqvec = arma::zeros(L,1);
  int n_l = n;
  for (int l=0; l<L; l++) {
    n_l = sum(c==l+1);
    ssqvec(l) = arma::as_scalar(arma::sum(square( X_out.rows(find(c == l+1)) - theta(l) * X_covars(find(c==l+1)) ) ));
  }
  double ssq = arma::sum(ssqvec);
  double alpha_star = alpha + 0.5*n;
  double beta_star = beta + 0.5 * ssq;
  return( pow(rgamma(1, alpha_star, 1.0/beta_star)(0), -0.5));
}

// [[Rcpp::export]]
arma::mat outer_prod(arma::mat A, arma::mat B) {
  return(A.t() * B);
}


// [[Rcpp::export]]
arma::mat gibbs_theta_sigma_conditional(arma::vec mu0, arma::mat kappa0,
                                      double alpha, double beta,
                                      arma::vec X_out,
                                      arma::mat X_in,
                                      arma::vec c,
                                      int L) {
  int n = X_out.n_rows;
  arma::mat X_covars = arma::join_horiz(arma::ones(n,1), X_in);
  int p = X_covars.n_cols;
  arma::mat theta = arma::zeros(L,p);
  arma::vec sigma = arma::zeros(L,1);
  arma::vec mu0_l = arma::zeros(p,1);
  arma::mat kappa0_l = arma::eye(p,p);
  double alpha_l = 1.0;
  double beta_l = 1.0;
  int n_l = X_out.n_rows;
  arma::vec X_out_l = X_out;
  arma::mat X_covars_l = X_covars;
  for (int l=0; l<L; l++) {
    n_l = sum(c == l+1);
    if (n_l==0) {
      sigma(l) = pow(arma::randg(arma::distr_param(alpha, pow(beta, -1.0))), -0.5);
      theta.row(l) = (mu0 + sigma(l)*arma::chol(arma::inv(kappa0)) * arma::randn(p,1)).t();
    } else {
      X_out_l = X_out.rows(find(c==l+1));
      X_covars_l = X_covars.rows(find(c==l+1));
      kappa0_l = kappa0 + outer_prod(X_covars_l, X_covars_l);
      mu0_l = inv(kappa0_l) * (outer_prod(X_covars_l, X_out_l) + kappa0*mu0);
      alpha_l = alpha + (n_l/2.0);
      beta_l = beta + 0.5 * (dot(X_out_l, X_out_l) + dot(mu0, kappa0*mu0) - dot(mu0_l, kappa0_l * mu0_l));
      sigma(l) = pow(arma::randg(arma::distr_param(alpha_l, pow(beta_l, -1.0))), -0.5);
      theta.row(l) = (mu0_l + sigma(l) * arma::chol(arma::inv(kappa0)) * arma::randn(p,1)).t();
    }
  }
  return(arma::join_horiz(theta,sigma));
}

// [[Rcpp::export]]
List gibbs_conditional(int n, arma::mat X,
           double gamma1, double gamma2,
           double mu01, double sigma01,
           double mu02, double sigma02,
           double alpha1, double beta1,
           double alpha2, double beta2,
           int L1, int L2,
           double a_rho, double b_rho,
           int R) {
  // initialization
  arma::vec q1 = rdirichlet_arma((gamma1/(L1*1.0)) * arma::ones(L1,1));
  arma::vec q2 = rdirichlet_arma( (gamma2/(L2*1.0)) * arma::ones(L2,1));
  arma::vec c1 = arma::zeros(n,1);
  arma::vec c2 = arma::zeros(n,1);
  for (int i=0; i<n; i++) {
    c1(i) = sample_arma(q1);
    c2(i) = sample_arma(q2);
  }
  arma::mat c1_out = arma::zeros(R,n);
  c1_out.row(0) = c1.t();
  arma::mat c2_out = arma::zeros(R,n);
  c2_out.row(0) = c2.t();
  arma::mat N = arma_table(c1,c2,L1,L2);
  arma::mat M = N;
  double M_sum = arma::accu(M);
  double rho = 1.0;
  arma::vec rho_out = arma::ones(R);
  arma::mat Cmat = arma::zeros(n,2);
  arma::vec theta1 = arma::zeros(L1,1);
  arma::vec theta2 = arma::zeros(L2,1);
  double sigma1 = 1.0;
  double sigma2 = 1.0;
  // sampling
  for (int r=1; r<R; r++) {
    // sample c
    Cmat = gibbs_C_conditional(c1,c2,rho,
                             L1,L2,theta1,theta2,
                             X.col(0),X.col(1),
                             q1,q2,
                             sigma1,sigma2);
    c1 = Cmat.col(0);
    c2 = Cmat.col(1);
    N = arma_table(c1,c2,L1,L2);
    // sample M
    M = gibbs_M(N,rho,q1,q2);
    // sample psi and psi_prime
    arma::vec M_r = arma::sum(M,1);
    q1 = rdirichlet_arma(M_r + (gamma1/(L1*1.0)) * arma::ones(L1,1));
    arma::rowvec M_c = arma::sum(M,0);
    q2 = rdirichlet_arma(M_c.t() + (gamma2/(L2*1.0)) * arma::ones(L2,1));
    // sample rho
    M_sum = arma::accu(M);
    rho = gibbs_rho(rho, n, M_sum, a_rho, b_rho);
    // sample theta
    theta1 = gibbs_theta_multiview(mu01, sigma01,
                                   sigma1,
                                   X.col(0),
                                   c1,
                                   L1);
    theta2 = gibbs_theta_conditional(mu02, sigma02,
                                   sigma2,
                                   X.col(1),
                                   X.col(0),
                                   c2,
                                   L2);
    // sample sigma
    sigma1 = gibbs_sigma_multiview(alpha1, beta1,
                                   theta1,
                                   X.col(0),
                                   c1,
                                   L1);
    sigma2 = gibbs_sigma_conditional(alpha2, beta2,
                                   theta2,
                                   X.col(1),
                                   X.col(0),
                                   c2,
                                   L2);
    // saving
    rho_out(r) = rho;
    c1_out.row(r) = c1.t();
    c2_out.row(r) = c2.t();
  }
  return(List::create(Named("rho") = rho_out,
                      Named("c1") = c1_out,
                      Named("c2") = c2_out));
}

// [[Rcpp::export]]
List gibbs_conditional_markov(int n, arma::mat X,
                       double gamma1, double gamma2,
                       double mu01, double sigma01,
                       double mu02, double sigma02,
                       double alpha1, double beta1,
                       double alpha2, double beta2,
                       int L1, int L2,
                       double a_rho, double b_rho,
                       int R) {
  // initialization
  arma::vec q1 = rdirichlet_arma((gamma1/(L1*1.0)) * arma::ones(L1,1));
  arma::vec q2 = rdirichlet_arma( (gamma2/(L2*1.0)) * arma::ones(L2,1));
  arma::vec c1 = arma::zeros(n,1);
  arma::vec c2 = arma::zeros(n,1);
  for (int i=0; i<n; i++) {
    c1(i) = sample_arma(q1);
    c2(i) = sample_arma(q2);
  }
  arma::mat c1_out = arma::zeros(R,n);
  c1_out.row(0) = c1.t();
  arma::mat c2_out = arma::zeros(R,n);
  c2_out.row(0) = c2.t();
  arma::mat N = arma_table(c1,c2,L1,L2);
  arma::mat M = N;
  double M_sum = arma::accu(M);
  double rho = 1.0;
  arma::vec rho_out = arma::ones(R);
  arma::mat Cmat = arma::zeros(n,2);
  arma::vec theta1 = arma::zeros(L1,1);
  arma::vec theta2 = arma::zeros(L2,1);
  double sigma1 = 1.0;
  double sigma2 = 1.0;
  // sampling
  for (int r=1; r<R; r++) {
    // sample c
    Cmat = gibbs_C_conditional_markov(c1,c2,rho,
                               L1,L2,theta1,theta2,
                               X.col(0),X.col(1),
                               q1,q2,
                               sigma1,sigma2);
    c1 = Cmat.col(0);
    c2 = Cmat.col(1);
    N = arma_table(c1,c2,L1,L2);
    // sample M
    M = gibbs_M(N,rho,q1,q2);
    // sample psi and psi_prime
    arma::vec M_r = arma::sum(M,1);
    q1 = rdirichlet_arma(M_r + (gamma1/(L1*1.0)) * arma::ones(L1,1));
    arma::rowvec M_c = arma::sum(M,0);
    q2 = rdirichlet_arma(M_c.t() + (gamma2/(L2*1.0)) * arma::ones(L2,1));
    // sample rho
    M_sum = arma::accu(M);
    rho = gibbs_rho(rho, n, M_sum, a_rho, b_rho);
    // sample theta
    theta1 = gibbs_theta_multiview(mu01, sigma01,
                                   sigma1,
                                   X.col(0),
                                   c1,
                                   L1);
    theta2 = gibbs_theta_conditional(mu02, sigma02,
                                     sigma2,
                                     X.col(1),
                                     X.col(0),
                                     c2,
                                     L2);
    // sample sigma
    sigma1 = gibbs_sigma_multiview(alpha1, beta1,
                                   theta1,
                                   X.col(0),
                                   c1,
                                   L1);
    sigma2 = gibbs_sigma_conditional(alpha2, beta2,
                                     theta2,
                                     X.col(1),
                                     X.col(0),
                                     c2,
                                     L2);
    // saving
    rho_out(r) = rho;
    c1_out.row(r) = c1.t();
    c2_out.row(r) = c2.t();
  }
  return(List::create(Named("rho") = rho_out,
                      Named("c1") = c1_out,
                      Named("c2") = c2_out));
}

// [[Rcpp::export]]
List gibbs_conditional_markov_ls(int n, int p, arma::mat X,
                              double gamma1, double gamma2,
                              arma::vec mu01, arma::mat sigma01,
                              arma::vec mu02, arma::mat kappa02,
                              double alpha1, double beta1,
                              double alpha2, double beta2,
                              int L1, int L2,
                              double a_rho, double b_rho,
                              int R) {
  // initialization
  arma::mat X_covars = arma::join_horiz(X.col(0), X.col(1));
  arma::vec q1 = rdirichlet_arma((gamma1/(L1*1.0)) * arma::ones(L1,1));
  arma::vec q2 = rdirichlet_arma( (gamma2/(L2*1.0)) * arma::ones(L2,1));
  arma::vec c1 = arma::zeros(n,1);
  arma::vec c2 = arma::zeros(n,1);
  for (int i=0; i<n; i++) {
    c1(i) = sample_arma(q1);
    c2(i) = sample_arma(q2);
  }
  arma::mat c1_out = arma::zeros(R,n);
  c1_out.row(0) = c1.t();
  arma::mat c2_out = arma::zeros(R,n);
  c2_out.row(0) = c2.t();
  arma::mat N = arma_table(c1,c2,L1,L2);
  arma::mat M = N;
  double M_sum = arma::accu(M);
  double rho = 1.0;
  arma::vec rho_out = arma::ones(R);
  arma::mat Cmat = arma::zeros(n,2);
  arma::vec mu2 = arma::zeros(L2,1);
  arma::mat theta1 = arma::zeros(L1,p);
  arma::mat theta2 = arma::zeros(L2,p);
  arma::vec sigma1 = arma::ones(p,1);
  arma::vec sigma2 = arma::ones(L2,1);
  arma::mat ts2 = arma::zeros(L2,p+1);
  // sampling
  for (int r=1; r<R; r++) {
    // sample c
    Cmat = gibbs_C_conditional_markov_ls(c1, c2, rho,
                                          L1, L2,
                                          theta1, theta2,
                                          mu2,
                                          X_covars, X.col(2),
                                          q1, q2,
                                          sigma1, sigma2);
    c1 = Cmat.col(0);
    c2 = Cmat.col(1);
    N = arma_table(c1,c2,L1,L2);
    // sample M
    M = gibbs_M(N,rho,q1,q2);
    // sample psi and psi_prime
    arma::vec M_r = arma::sum(M,1);
    q1 = rdirichlet_arma(M_r + (gamma1/(L1*1.0)) * arma::ones(L1,1));
    arma::rowvec M_c = arma::sum(M,0);
    q2 = rdirichlet_arma(M_c.t() + (gamma2/(L2*1.0)) * arma::ones(L2,1));
    // sample rho
    M_sum = arma::accu(M);
    rho = gibbs_rho(rho, n, M_sum, a_rho, b_rho);
    // sample theta and sigma
    theta1 = gibbs_theta_multiview_bivariate(mu01, sigma01,
                                             arma::diagmat(sigma1),
                                             X_covars,
                                             c1,
                                             L1);
    sigma1(0) = gibbs_sigma_multiview(alpha1, beta1,
                                      theta1.col(0),
                                      X_covars.col(0),
                                      c1,
                                      L1);
    sigma1(1) = gibbs_sigma_multiview(alpha1, beta1,
                                      theta1.col(1),
                                      X_covars.col(1),
                                      c1,
                                      L1);
    ts2 = gibbs_theta_sigma_conditional(mu02, kappa02,
                                        alpha2, beta2,
                                        X.col(2),
                                        X_covars,
                                        c2,
                                        L2);
    mu2 = ts2.col(0);
    theta2 = arma::join_horiz(ts2.col(1), ts2.col(2));
    sigma2 = ts2.col(3);
    // saving
    rho_out(r) = rho;
    c1_out.row(r) = c1.t();
    c2_out.row(r) = c2.t();
  }
  return(List::create(Named("rho") = rho_out,
                      Named("c1") = c1_out,
                      Named("c2") = c2_out));
}


// [[Rcpp::export]]
arma::vec gibbs_C_profile(arma::vec c,
                          int L,
                          arma::vec theta,
                          arma::vec X_out,
                          arma::vec X_covars,
                          double gam,
                          double sigma){
  int n = X_out.n_rows;
  int n_k = n;
  arma::vec lp = arma::zeros(L,1);
  arma::vec p = arma::zeros(L,1);
  for (int i=0; i<n; i++) {
    for (int k=0; k<L; k++) {
      if (c(i)==k+1) {
        n_k = arma::sum(c==k+1)-1;
      } else {
        n_k = arma::sum(c==k+1);
      }
      lp(k) = log(n_k + (gam/(L*1.0))) + arma::log_normpdf(X_out(i), X_covars(i)*theta(k), sigma);
    }
    p = arma::exp(lp - arma::max(lp));
    c(i) = sample_arma(p/sum(p));
  }
  return(c);
}

// [[Rcpp::export]]
arma::mat gibbs_vanilla_profile(int n, arma::vec X_out, arma::vec X_covars,
                        double gam,
                        double mu0, double sigma0,
                        double alpha, double beta,
                        int L,
                        int R) {
  // initialization
  arma::vec q = rdirichlet_arma((gam/(L*1.0)) * arma::ones(L,1));
  arma::vec c = arma::zeros(n,1);
  for (int i=0; i<n; i++) {
    c(i) = sample_arma(q);
  }
  arma::mat c_out = arma::zeros(R,n);
  c_out.row(0) = c.t();
  arma::vec theta = arma::zeros(L,1);
  double sigma = 1.0;
  // sampling
  for (int r=1; r<R; r++) {
    // sample c
    c = gibbs_C_profile(c,
                        L,
                        theta,
                        X_out,
                        X_covars,
                        gam,
                        sigma);
    // sample theta
    theta = gibbs_theta_conditional(mu0, sigma0,
                                    sigma,
                                    X_out,
                                    X_covars,
                                    c,
                                    L);
    // sample sigma
    sigma = gibbs_sigma_conditional(alpha, beta,
                                    theta,
                                    X_out,
                                    X_covars,
                                    c,
                                    L);
    // saving
    c_out.row(r) = c.t();
  }
  return(c_out);
}



// [[Rcpp::export]]
double debug(arma::vec x, arma::vec mu, arma::vec sigma) {
  return(arma::as_scalar(arma::sum(arma::log_normpdf(x, mu, sigma))));
}

// [[Rcpp::export]]
arma::mat gibbs_C_multivariate(arma::vec c1, arma::vec c2, double rho,
                               int L1, int L2,
                               arma::mat theta1, arma::mat theta2,
                               arma::mat X1, arma::mat X2,
                               arma::vec q1, arma::vec q2,
                               arma::vec sigma1, arma::vec sigma2){
  int n = X1.n_rows;
  int n_k1 = n;
  int n_k1k2 = n;
  arma::vec lp1 = arma::zeros(L1,1);
  arma::vec p1 = arma::zeros(L1,1);
  arma::vec lp2 = arma::zeros(L2,1);
  arma::vec p2 = arma::zeros(L2,1);
  for (int i=0; i<n; i++){
    // first sample c1
    for (int k1=0; k1<L1;k1++) {
      if (c1(i)==k1+1) {
        n_k1 = arma::sum(c1==k1+1) - 1;
      } else {
        n_k1 = arma::sum(c1==k1+1);
      }
      lp1(k1) = log(rho * q1(k1) + n_k1) +
        arma::as_scalar(arma::sum(arma::log_normpdf(X1.row(i).t(), theta1.row(k1).t(), sigma1)));
    }
    p1 = arma::exp(lp1 - arma::max(lp1));
    c1(i) = sample_arma(p1/arma::sum(p1));
    // next sample c2
    for (int k2=0; k2<L2; k2++) {
      if (c2(i)==k2+1) {
        n_k1k2 =  arma::sum( (c1 == c1(i)) % (c2== k2+1) ) -1;
      } else {
        n_k1k2 =  arma::sum( (c1 == c1(i)) % (c2== k2+1) );
      }
      lp2(k2) = log(rho*q1(c1(i)-1)*q2(k2) + n_k1k2) +
        arma::as_scalar(arma::sum(arma::log_normpdf(X2.row(i).t(), theta2.row(k2).t(), sigma2)));
    }
    p2 = arma::exp(lp2 - arma::max(lp2));
    c2(i) = sample_arma(p2/arma::sum(p2));
  }
  return(arma::join_horiz(c1,c2));
}

// [[Rcpp::export]]
arma::mat gibbs_theta_multivariate(arma::vec mu0, 
                                  arma::mat sigma0,
                                  arma::mat sigma,
                                  arma::mat X,
                                  arma::vec c,
                                  int L) {
  int p = X.n_cols;
  arma::mat theta = arma::zeros(L,p);
  arma::vec mu_l = arma::zeros(p,1);
  arma::mat sigma_l = arma::eye(p,p);
  arma::mat sigma0_inv = arma::inv_sympd(sigma0);
  arma::mat sigma_inv = arma::inv_sympd(sigma);
  int n_l = X.n_rows;
  arma::mat X_l = X;
  for (int l=0; l<L; l++) {
    n_l = sum(c == l+1);
    if (n_l==0) {
      theta.row(l) = (mu0 + chol(sigma0)*arma::randn(p,1)).t();
    } else {
      X_l = X.rows(find(c==l+1));
      sigma_l = arma::inv_sympd( sigma0_inv + n_l *  sigma_inv );
      mu_l = sigma_l * ( sigma0_inv * mu0 + n_l * sigma_inv * arma::mean(X_l,0).t() );
      theta.row(l) = (mu_l + chol(sigma_l) * arma::randn(p,1)).t();
    }
  }
  return(theta);
}

// [[Rcpp::export]]
double gibbs_sigma_multivariate(double alpha, double beta,
                                arma::mat theta,
                                arma::mat X,
                                arma::vec c,
                                int L) {
  int n = X.n_rows;
  int p = X.n_cols;
  arma::vec ssqvec = arma::zeros(n,1);
  for (int i=0; i<n; i++) {
    ssqvec(i) = arma::sum(arma::square(X.row(i) - theta.row(c(i)-1)));
  }
  double ssq = arma::sum(ssqvec);
  double alpha_star = alpha + 0.5*n*p;
  double beta_star = beta + 0.5 * ssq;
  return( pow(rgamma(1, alpha_star, 1.0/beta_star)(0), -0.5));
}


// [[Rcpp::export]]
List gibbs_markov_multivariate(int n,
                               arma::mat X1, arma::mat X2,
                               double gamma1, double gamma2,
                               arma::vec mu01, arma::mat sigma01,
                               arma::vec mu02, arma::mat sigma02,
                               double alpha1, double beta1,
                               double alpha2, double beta2,
                               int L1, int L2,
                               double a_rho, double b_rho,
                               int R) {
  // initialization
  int p1 = X1.n_cols;
  int p2 = X2.n_cols;
  arma::vec q1 = rdirichlet_arma((gamma1/(L1*1.0)) * arma::ones(L1,1));
  arma::vec q2 = rdirichlet_arma( (gamma2/(L2*1.0)) * arma::ones(L2,1));
  arma::vec c1 = arma::zeros(n,1);
  arma::vec c2 = arma::zeros(n,1);
  for (int i=0; i<n; i++) {
    c1(i) = sample_arma(q1);
    c2(i) = sample_arma(q2);
  }
  arma::mat c1_out = arma::zeros(R,n);
  c1_out.row(0) = c1.t();
  arma::mat c2_out = arma::zeros(R,n);
  c2_out.row(0) = c2.t();
  arma::mat N = arma_table(c1,c2,L1,L2);
  arma::mat M = N;
  double M_sum = arma::accu(M);
  double rho = 1.0;
  arma::vec rho_out = arma::ones(R);
  arma::mat Cmat = arma::zeros(n,2);
  arma::mat theta1 = arma::zeros(L1,p1);
  arma::mat theta2 = arma::zeros(L2,p2);
  double sigma1 = 1.0;
  double sigma2 = 1.0;
  // sampling
  for (int r=1; r<R; r++) {
    // sample c
    Cmat = gibbs_C_multivariate(c1, c2, rho,
                                L1, L2,
                                theta1, theta2,
                                X1, X2,
                                q1, q2,
                                sigma1*arma::ones(p1), sigma2*arma::ones(p2));
    c1 = Cmat.col(0);
    c2 = Cmat.col(1);
    N = arma_table(c1,c2,L1,L2);
    // sample M
    M = gibbs_M(N,rho,q1,q2);
    // sample psi and psi_prime
    arma::vec M_r = arma::sum(M,1);
    q1 = rdirichlet_arma(M_r + (gamma1/(L1*1.0)) * arma::ones(L1,1));
    arma::rowvec M_c = arma::sum(M,0);
    q2 = rdirichlet_arma(M_c.t() + (gamma2/(L2*1.0)) * arma::ones(L2,1));
    // sample rho
    M_sum = arma::accu(M);
    rho = gibbs_rho(rho, n, M_sum, a_rho, b_rho);
    // sample theta
    theta1 = gibbs_theta_multivariate(mu01, 
                                      sigma01,
                                      pow(sigma1,2.0)*arma::eye(p1,p1),
                                      X1,
                                      c1,
                                      L1);
    theta2 = gibbs_theta_multivariate(mu02, 
                                      sigma02,
                                      pow(sigma2,2.0)*arma::eye(p2,p2),
                                      X2,
                                      c2,
                                      L2);
    // sample sigma
    sigma1 = gibbs_sigma_multivariate(alpha1, beta1,
                                      theta1,
                                      X1,
                                      c1,
                                      L1);
    sigma2 = gibbs_sigma_multivariate(alpha2, beta2,
                                      theta2,
                                      X2,
                                      c2,
                                      L2);
    // saving
    rho_out(r) = rho;
    c1_out.row(r) = c1.t();
    c2_out.row(r) = c2.t();
  }
  return(List::create(Named("rho") = rho_out,
                      Named("c1") = c1_out,
                      Named("c2") = c2_out));
}

// [[Rcpp::export]]
arma::vec gibbs_C_vanilla_multivariate(arma::vec c,
                          int L,
                          arma::mat theta,
                          arma::mat X,
                          double gam,
                          arma::vec sigma){
  int n = X.n_rows;
  int n_k = n;
  arma::vec lp = arma::zeros(L,1);
  arma::vec p = arma::zeros(L,1);
  for (int i=0; i<n; i++) {
    for (int k=0; k<L; k++) {
      if (c(i)==k+1) {
        n_k = arma::sum(c==k+1)-1;
      } else {
        n_k = arma::sum(c==k+1);
      }
      lp(k) = log(n_k + (gam/(L*1.0))) +
        arma::as_scalar(arma::sum(arma::log_normpdf(X.row(i).t(), theta.row(k).t(), sigma)));
    }
    p = arma::exp(lp - arma::max(lp));
    c(i) = sample_arma(p/sum(p));
  }
  return(c);
}

// [[Rcpp::export]]
arma::mat gibbs_vanilla_multivariate(int n, 
                                     arma::mat X,
                                     double gam,
                                     arma::vec mu0, arma::mat sigma0,
                                     double alpha, double beta,
                                     int L, int R) {
  // initialization
  int p = X.n_cols;
  arma::vec q = rdirichlet_arma((gam/(L*1.0)) * arma::ones(L,1));
  arma::vec c = arma::zeros(n,1);
  for (int i=0; i<n; i++) {
    c(i) = sample_arma(q);
  }
  arma::mat c_out = arma::zeros(R,n);
  c_out.row(0) = c.t();
  arma::mat theta = arma::zeros(L,p);
  double sigma = 1.0;
  // sampling
  for (int r=1; r<R; r++) {
    // sample c
    c = gibbs_C_vanilla_multivariate(c,
                                    L,
                                    theta,
                                    X,
                                    gam,
                                    sigma*arma::ones(p));
    // sample theta
    theta = gibbs_theta_multivariate(mu0, 
                                    sigma0,
                                    pow(sigma,2.0)*arma::eye(p,p),
                                    X,
                                    c,
                                    L);
    // sample sigma
    sigma = gibbs_sigma_multivariate(alpha, beta,
                                     theta,
                                     X,
                                     c,
                                     L);
    // saving
    c_out.row(r) = c.t();
  }
  return(c_out);
}
