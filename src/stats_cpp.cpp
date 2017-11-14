// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <cmath>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::rowvec stats_SE_cpp(arma::colvec y, arma::mat& Kmat, arma::mat& invKmatn, arma::vec& eigenval, const double mu) {
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
  stats(1) = - 0.5 * (n * log( 2.0 * arma::datum::pi ) + sum(log(eigenval)) + arma::dot( y, alpha ) ) ;

  return stats;
}

