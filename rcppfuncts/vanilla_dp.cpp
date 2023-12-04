#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec gibbs_C_vanilla(arma::vec c,
                            int L,
                            arma::vec theta,
                            arma::vec X,
                            double gam,
                            double sigma){
  int n = X.n_rows;
  int n_k = n;
  arma::mat lp = arma::zeros(L,1);
  arma::mat p = arma::zeros(L,1);
  for (int i=0; i<n; i++) {
    for (int k=0; k<L; k++) {
      if (c(i)==k+1) {
        n_k = arma::sum(c==k+1)-1;
      } else {
        n_k = arma::sum(c==k+1);
      }
      lp(k) = log(n_k + (gam/L)) + arma::log_normpdf(X(i), theta(k), sigma);
    }
    p = arma::exp(lp - arma::max(lp))/(arma::sum(arma::exp(lp - arma::max(lp))));
    c(i) = sample_arma(p);
  }
  return(c);
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
