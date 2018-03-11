// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include "ace_kernel_utils.hpp"

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::rowvec stats_cpp(const arma::colvec& y,
                       const arma::mat& Kmat,
                       const arma::mat& invKmatn,
                       const arma::vec& eigenval,
                       const double mu,
                       double std_y = 1) {
  //Rcpp::List gradients = clone(parameters); //for the same list structure

  unsigned int n = y.n_rows;
  //preallocate memory
  arma::rowvec stats(2);
  arma::colvec alpha(n);

  arma::colvec ybar = y - mu; //
  alpha = invKmatn * ybar;

  //RMSE
  stats(0) = std_y * arma::norm(ybar - (Kmat * alpha)) / sqrt(n);

  //Evidence
  stats(1) = logevidence(y, alpha, eigenval, n);

  return stats;
}

