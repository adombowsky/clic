#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double rand_index(arma::rowvec c1, arma::rowvec c2) {
  int n = c1.n_cols;
  arma::mat RI_mat = arma::zeros(n,n);
  for (int i=0; i<(n-1); i++) {
    for (int j=i+1; j<n; j++) {
      RI_mat(i,j) = ((c1(i) == c1(j)) & (c2(i) == c2(j))) + ((c1(i) != c1(j)) & (c2(i) != c2(j)));
    }
  }
  return(arma::accu(RI_mat)/(0.5 * n * (n-1)));
}

// [[Rcpp::export]]
arma::vec rand_index_MCMC(arma::mat c1, arma::mat c2) {
  int R = c1.n_rows;
  arma::vec RI = arma::zeros(R,1);
  for (int r=0; r<R; r++) {
    RI(r) = rand_index(c1.row(r), c2.row(r));
  }
  return(RI);
}