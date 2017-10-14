// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat normalize_train(arma::vec& y, arma::mat& X, arma::mat& Z) {
  unsigned int px = X.n_cols;
  unsigned int pz = Z.n_cols;

  arma::mat moments(1+px+pz,2);
  moments.col(0).zeros(); moments.col(1).ones();

  //create index of binary variables in X and Z
  arma::uvec isbinary(px+pz);
  for(unsigned int i = 0; i < px; i++){
    isbinary(i) = (arma::sum(arma::pow(arma::unique( X.col(i) ),0)) == 2); // ugly
  }
  for(unsigned int i = px; i < (px+pz); i++){
    isbinary(i) = (sum(pow(arma::unique( Z.col(i) ),0)) == 2); // ugly
  }

  //get mean and de-mean
  moments(0,0) = 0; // mean is a parameter
  y = y - moments(0,0);
  for(unsigned int i = 1; i < (px+1); i++){
    if(isbinary(i-1)==0){
      moments(i,0) = mean(X.col(i));
      X.col(i) = X.col(i) - moments(i,0);
    }
  }
  for(unsigned int i = px+1; i < (px+pz+1); i++){
    if(isbinary(i-1)==0){
      moments(i,0) = mean(Z.col(i));
      Z.col(i) = Z.col(i) - moments(i,0);
    }
  }

  //rescale
  moments(0,1) = arma::stddev(y); // mean is a parameter
  y = y / moments(0,1);
  for(unsigned int i = 1; i < (px+1); i++){
    if(isbinary(i-1)==0){
      moments(i,1) = arma::max(abs(X.col(i)));
      X.col(i) = X.col(i) / moments(i,1);
    }
  }
  for(unsigned int i = px+1; i < (px+pz+1); i++){
    if(isbinary(i-1)==0){
      moments(i,1) = arma::max(abs(Z.col(i)));
      Z.col(i) = Z.col(i) / moments(i,1);
    }
  }


  //return moments
  return moments;
}


// [[Rcpp::export]]
void normalize_test(arma::vec& y, arma::mat& X, arma::mat& Z, const arma::mat& moments) {
  unsigned int px = X.n_cols;
  unsigned int pz = Z.n_cols;

  y = (y - moments(0,0)) / moments(0,1);
  for(unsigned int i = 0 ; i < px; i++){
    X.col(i) = (X.col(i) - moments(i+1,0)) / moments(i+1,1);
  }
  for(unsigned int i = 0 ; i < pz; i++){
    Z.col(i) = (Z.col(i) - moments(i+1+px,0))/ moments(i+1+px,1);
  }
}
