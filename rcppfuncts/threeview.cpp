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
int gibbs_m_kkk(int n_kkk, double rho, double psi_k, double psi_k_prime, double psi_k_tilde) {
  arma::mat ls = lstirling(n_kkk);
  arma::vec lprbs = arma::zeros(n_kkk+1);
  for (int m = 0; m <= n_kkk; m++) {
    lprbs(m) = ls(n_kkk, m) + m * (log(rho) + log(psi_k) + log(psi_k_prime) + log(psi_k_tilde));
  }
  arma::vec prbs = arma::exp(lprbs - arma::max(lprbs));
  int m_kkk = sample_arma(prbs/arma::sum(prbs)) - 1;
  return(m_kkk);
}


// [[Rcpp::export]]
arma::cube gibbs_M_tv(arma::cube N, double rho, arma::vec Psi, arma::vec Psi_prime, arma::vec Psi_tilde) {
  int K = Psi.n_rows;
  int K_prime = Psi_prime.n_rows;
  int K_tilde = Psi_tilde.n_rows;
  arma::cube M = arma::zeros(K, K_prime, K_tilde);
  for (int k=0; k<K; k++) {
    for (int k_prime=0; k_prime<K_prime; k_prime++) {
      for (int k_tilde=0; k_tilde<K_tilde;k_tilde++) {
        if (N(k,k_prime,k_tilde)==0) {
          continue;
        } else {
          M(k, k_prime, k_tilde) = gibbs_m_kkk(N(k,k_prime,k_tilde), rho, Psi(k), Psi_prime(k_prime), Psi_tilde(k_tilde));
        }
      }
    }
  }
  return(M);
}

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
double gibbs_rho_grid(arma::vec rho_grid, int n, arma::cube N, arma::vec q1, arma::vec q2, arma::vec q3) {
  int W = rho_grid.n_rows;
  int K1 = N.n_rows;
  int K2 = N.n_cols;
  int K3 = N.n_slices;
  arma::vec lp = arma::zeros(W,1);
  for (int w=0; w<W; w++) {
    lp(w) = lgamma(rho_grid(w))-lgamma(rho_grid(w) + n);
    for (int k1 = 0; k1<K1; k1++) {
      for (int k2 = 0; k2<K2; k2++) {
        for (int k3=0; k3<K3; k3++) {
          lp(w) = lp(w) +
            lgamma(N(k1,k2,k3) + rho_grid(w) * q1(k1) * q2(k2) * q3(k3))-
            lgamma(rho_grid(w) * q1(k1) * q2(k2) * q3(k3));
        }
      }
    }
  }
  arma::vec p = exp(lp - arma::max(lp));
  int rho_index = sample_arma(p/arma::sum(p)) - 1;
  return(rho_grid(rho_index));
}


// [[Rcpp::export]]
arma::mat gibbs_C_multiview_markov_tv(arma::vec c1, arma::vec c2, arma::vec c3, double rho,
                                   int L1, int L2, int L3,
                                   arma::vec theta1, arma::vec theta2, arma::vec theta3,
                                   arma::vec X1, arma::vec X2, arma::vec X3,
                                   arma::vec q1, arma::vec q2, arma::vec q3,
                                   double sigma1, double sigma2, double sigma3){
  int n = X1.n_rows;
  int n_k1 = n;
  int n_k1k2 = n;
  int n_k1k2k3 = n;
  arma::vec lp1 = arma::zeros(L1,1);
  arma::vec p1 = arma::zeros(L1,1);
  arma::vec lp2 = arma::zeros(L2,1);
  arma::vec p2 = arma::zeros(L2,1);
  arma::vec lp3 = arma::zeros(L3,1);
  arma::vec p3 = arma::zeros(L3,1);
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
    // next sample c3
    for (int k3=0; k3<L3; k3++) {
      if (c3(i)==k3+1) {
        n_k1k2k3 =  arma::sum( (c1 == c1(i)) % (c2 == c2(i)) % (c3== k3+1) ) -1;
      } else {
        n_k1k2k3 =  arma::sum( (c1 == c1(i)) % (c2 == c2(i)) % (c3== k3+1) );
      }
      lp3(k3) = log(rho*q1(c1(i)-1)*q2(c2(i)-1)*q3(k3) + n_k1k2k3) + arma::log_normpdf(X3(i), theta3(k3), sigma3);
    }
    p3 = arma::exp(lp3 - arma::max(lp3));
    c3(i) = sample_arma(p3/arma::sum(p3));
  }
  return(arma::join_horiz(c1,c2,c3));
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
arma::cube arma_array(arma::vec c1, arma::vec c2, arma::vec c3, int L1, int L2, int L3){
  arma::cube N = arma::zeros(L1,L2,L3);
  for (int k1 = 0; k1<L1; k1++) {
    for (int k2 = 0; k2<L2; k2++) {
      for (int k3=0; k3<L3;k3++) {
        N(k1,k2,k3) = sum((c1==k1+1) % (c2==k2+1) % (c3==k3+1));
      }
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
arma::vec array_rowsum(arma::cube N) {
  arma::cube sl = arma::sum(arma::sum(N,2),1);
  return(sl.slice(0));
}

// [[Rcpp::export]]
arma::vec array_colsum(arma::cube N) {
  arma::cube sl = arma::sum(arma::sum(N,2),0);
  return(sl.slice(0).t());
}

// [[Rcpp::export]]
arma::vec array_slicesum(arma::cube N) {
  int L = N.n_slices;
  arma::vec sl = arma::zeros(L,1);
  for (int l=0; l<L; l++) {
    sl(l) = arma::accu(N.slice(l));
  }
  return(sl);
}

// [[Rcpp::export]]
List gibbs_markov_tv(int n, arma::mat X,
                  double gamma1, double gamma2, double gamma3,
                  double mu01, double sigma01,
                  double mu02, double sigma02,
                  double mu03, double sigma03,
                  double alpha1, double beta1,
                  double alpha2, double beta2,
                  double alpha3, double beta3,
                  int L1, int L2, int L3,
                  double a_rho, double b_rho,
                  int R) {
  // initialization
  arma::vec q1 = rdirichlet_arma((gamma1/(L1*1.0)) * arma::ones(L1,1));
  arma::vec q2 = rdirichlet_arma( (gamma2/(L2*1.0)) * arma::ones(L2,1));
  arma::vec q3 = rdirichlet_arma( (gamma3/(L3*1.0)) * arma::ones(L3,1));
  arma::vec c1 = arma::zeros(n,1);
  arma::vec c2 = arma::zeros(n,1);
  arma::vec c3 = arma::zeros(n,1);
  for (int i=0; i<n; i++) {
    c1(i) = sample_arma(q1);
    c2(i) = sample_arma(q2);
    c3(i) = sample_arma(q3);
  }
  arma::mat c1_out = arma::zeros(R,n);
  c1_out.row(0) = c1.t();
  arma::mat c2_out = arma::zeros(R,n);
  c2_out.row(0) = c2.t();
  arma::mat c3_out = arma::zeros(R,n);
  c3_out.row(0) = c3.t();
  arma::cube N = arma_array(c1,c2,c3,L1,L2,L3);
  arma::cube M = N;
  double M_sum = arma::accu(M);
  double rho = 1.0;
  arma::vec rho_out = arma::ones(R);
  arma::mat Cmat = arma::zeros(n,3);
  arma::vec theta1 = arma::zeros(L1,1);
  arma::vec theta2 = arma::zeros(L2,1);
  arma::vec theta3 = arma::zeros(L3,1);
  double sigma1 = 1.0;
  double sigma2 = 1.0;
  double sigma3 = 1.0;
  // sampling
  for (int r=1; r<R; r++) {
    // sample c
    Cmat = gibbs_C_multiview_markov_tv(c1, c2, c3, rho,
            L1, L2, L3,
            theta1, theta2, theta3,
            X.col(0), X.col(1), X.col(2),
            q1, q2, q3,
            sigma1, sigma2, sigma3);
    c1 = Cmat.col(0);
    c2 = Cmat.col(1);
    c3 = Cmat.col(2);
    N = arma_array(c1,c2,c3,L1,L2,L3);
    // sample M
    M = gibbs_M_tv(N, rho, q1, q2, q3);
    // sample psi and psi_prime
    arma::vec M_r = array_rowsum(M);
    q1 = rdirichlet_arma(M_r + (gamma1/(L1*1.0)) * arma::ones(L1,1));
    arma::vec M_c = array_colsum(M);
    q2 = rdirichlet_arma(M_c + (gamma2/(L2*1.0)) * arma::ones(L2,1));
    arma::vec M_s = array_slicesum(M);
    q3 = rdirichlet_arma(M_s + (gamma3/(L3*1.0)) * arma::ones(L3,1));
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
    theta3 = gibbs_theta_multiview(mu03, sigma03,
                                   sigma3,
                                   X.col(2),
                                   c3,
                                   L3);
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
    sigma3 = gibbs_sigma_multiview(alpha3, beta3,
                                   theta3,
                                   X.col(2),
                                   c3,
                                   L3);
    // saving
    rho_out(r) = rho;
    c1_out.row(r) = c1.t();
    c2_out.row(r) = c2.t();
    c3_out.row(r) = c3.t();
  }
  return(List::create(Named("rho") = rho_out,
                      Named("c1") = c1_out,
                      Named("c2") = c2_out,
                      Named("c3") = c3_out));
}

// [[Rcpp::export]]
List grid_gibbs_markov_tv(int n, arma::mat X,
                     double gamma1, double gamma2, double gamma3,
                     double mu01, double sigma01,
                     double mu02, double sigma02,
                     double mu03, double sigma03,
                     double alpha1, double beta1,
                     double alpha2, double beta2,
                     double alpha3, double beta3,
                     int L1, int L2, int L3,
                     arma::vec rho_grid,
                     int R) {
  // initialization
  arma::vec q1 = rdirichlet_arma((gamma1/(L1*1.0)) * arma::ones(L1,1));
  arma::vec q2 = rdirichlet_arma( (gamma2/(L2*1.0)) * arma::ones(L2,1));
  arma::vec q3 = rdirichlet_arma( (gamma3/(L3*1.0)) * arma::ones(L3,1));
  arma::vec c1 = arma::zeros(n,1);
  arma::vec c2 = arma::zeros(n,1);
  arma::vec c3 = arma::zeros(n,1);
  for (int i=0; i<n; i++) {
    c1(i) = sample_arma(q1);
    c2(i) = sample_arma(q2);
    c3(i) = sample_arma(q3);
  }
  arma::mat c1_out = arma::zeros(R,n);
  c1_out.row(0) = c1.t();
  arma::mat c2_out = arma::zeros(R,n);
  c2_out.row(0) = c2.t();
  arma::mat c3_out = arma::zeros(R,n);
  c3_out.row(0) = c3.t();
  arma::cube N = arma_array(c1,c2,c3,L1,L2,L3);
  arma::cube M = N;
  double M_sum = arma::accu(M);
  double rho = 1.0;
  arma::vec rho_out = arma::ones(R);
  arma::mat Cmat = arma::zeros(n,3);
  arma::vec theta1 = arma::zeros(L1,1);
  arma::vec theta2 = arma::zeros(L2,1);
  arma::vec theta3 = arma::zeros(L3,1);
  double sigma1 = 1.0;
  double sigma2 = 1.0;
  double sigma3 = 1.0;
  // sampling
  for (int r=1; r<R; r++) {
    // sample c
    Cmat = gibbs_C_multiview_markov_tv(c1, c2, c3, rho,
                                       L1, L2, L3,
                                       theta1, theta2, theta3,
                                       X.col(0), X.col(1), X.col(2),
                                       q1, q2, q3,
                                       sigma1, sigma2, sigma3);
    c1 = Cmat.col(0);
    c2 = Cmat.col(1);
    c3 = Cmat.col(2);
    N = arma_array(c1,c2,c3,L1,L2,L3);
    // sample M
    M = gibbs_M_tv(N, rho, q1, q2, q3);
    // sample psi and psi_prime
    arma::vec M_r = array_rowsum(M);
    q1 = rdirichlet_arma(M_r + (gamma1/(L1*1.0)) * arma::ones(L1,1));
    arma::vec M_c = array_colsum(M);
    q2 = rdirichlet_arma(M_c + (gamma2/(L2*1.0)) * arma::ones(L2,1));
    arma::vec M_s = array_slicesum(M);
    q3 = rdirichlet_arma(M_s + (gamma3/(L3*1.0)) * arma::ones(L3,1));
    // sample rho
    M_sum = arma::accu(M);
    rho = gibbs_rho_grid(rho_grid, n, N, q1, q2, q3);;
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
    theta3 = gibbs_theta_multiview(mu03, sigma03,
                                   sigma3,
                                   X.col(2),
                                   c3,
                                   L3);
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
    sigma3 = gibbs_sigma_multiview(alpha3, beta3,
                                   theta3,
                                   X.col(2),
                                   c3,
                                   L3);
    // saving
    rho_out(r) = rho;
    c1_out.row(r) = c1.t();
    c2_out.row(r) = c2.t();
    c3_out.row(r) = c3.t();
  }
  return(List::create(Named("rho") = rho_out,
                      Named("c1") = c1_out,
                      Named("c2") = c2_out,
                      Named("c3") = c3_out));
}


