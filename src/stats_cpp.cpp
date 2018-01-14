// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include "ace_kernelutilities.hpp"

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::rowvec stats_cpp(arma::colvec y, arma::mat& Kmat, arma::mat& invKmatn, arma::vec& eigenval, const double mu) {
  //Rcpp::List gradients = clone(parameters); //for the same list structure

  unsigned int n = y.size();
  //preallocate memory
  arma::rowvec stats(2);
  arma::colvec alpha(n);

  y = y - mu; //
  alpha = invKmatn * y;

  //RMSE
  stats(0) = pow(arma::norm(y - (Kmat * alpha)),2);
  //Evidence
  stats(1) = logevidence(y, alpha, eigenval, n);

  return stats;
}

