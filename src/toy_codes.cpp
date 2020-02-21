// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export]]
arma::vec cpp_colsums(arma::mat A){
  int N = A.n_rows;
  int P = A.n_cols;
  
  double tmp = 0.0;
  arma::vec output(P,fill::zeros);
  for (int p=0;p<P;p++){
    tmp = 0.0;
    for (int n=0;n<N;n++){
      tmp += A(n,p);
    }
    output(p) = tmp;
  }
  return(output);
}