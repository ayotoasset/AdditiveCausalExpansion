// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat normalize_train(arma::vec& y, arma::mat& X, arma::mat& Z) {
  unsigned int px = X.n_cols;
  unsigned int pz = Z.n_cols;

  arma::mat moments(1+px+pz,3);
  //Rcpp::colnames(moments) = CharacterVector::create("mean", "scale", "isbinary"); // for Rcpp objects
  moments.col(0).zeros(); moments.col(1).ones(); moments.col(2).zeros();

  //create index of binary variables in X and Z
  arma::uvec isbinary(px+pz);
  unsigned long long tmp;
  for(unsigned int i = 0; i < px; i++){
    tmp = sum(pow(arma::unique( X.col(i) ),0));
    isbinary(i) = moments(i+1,2) = (tmp <= 2); // ugly but works
    if(tmp==1){
      Rcout << "Column " << i+1 << " of X is constant." << std::endl;
      X.col(i).zeros();
    }
  }
  for(unsigned int i = px; i < (px+pz); i++){
    tmp = sum(pow(arma::unique( Z.col(i-px) ),0));
    isbinary(i) = moments(i+1,2) = (tmp <= 2); // ugly
    if(tmp==1){
      Rcout << "Column " << i << " of Z is constant." << std::endl;
      Z.col(i).zeros();
    }
  }

  //get mean and de-mean
  moments(0,0) = mean(y);//0; // mean is a parameter
  y = y - moments(0,0);
  for(unsigned int i = 1; i < (px+1); i++){
    if(isbinary(i-1)==0){
      moments(i,0) = mean(X.col(i-1));
      X.col(i-1) = X.col(i-1) - moments(i,0);
    }
  }
  for(unsigned int i = px+1; i < (px+pz+1); i++){
    if(isbinary(i-1)==0){
      moments(i,0) = mean(Z.col(i-px-1));
      Z.col(i-px-1) = Z.col(i-px-1) - moments(i,0);
    }
  }
  //rescale
  moments(0,1) = arma::stddev(y); // mean is a parameter
  y = y / moments(0,1);
  for(unsigned int i = 1; i < (px+1); i++){
    if(isbinary(i-1)==0){
      moments(i,1) = arma::max(abs(X.col(i-1)));
      X.col(i-1) = X.col(i-1) / moments(i,1);
    } else {
      moments(i,1) = 1;//arma::stddev(X.col(i-1));
      X.col(i-1) = X.col(i-1) / moments(i,1);
    }
  }
  //important for Z to be within unit circle for splines
  for(unsigned int i = px+1; i < (px+pz+1); i++){
    if(isbinary(i-px-1)==0){
      moments(i,1) = arma::max(abs(Z.col(i-px-1)));
      Z.col(i-px-1) = Z.col(i-px-1) / moments(i,1);
    }
  }
  //return moments
  return moments;
}


// [[Rcpp::export]]
void normalize_test(arma::mat& X, arma::mat& Z, const arma::mat& moments) {
  unsigned int px = X.n_cols;
  unsigned int pz = Z.n_cols;

  //y = (y - moments(0,0)) / moments(0,1);
  for(unsigned int i = 0 ; i < px; i++){
    X.col(i) = (X.col(i) - moments(i+1,0)) / moments(i+1,1);
  }
  for(unsigned int i = 0 ; i < pz; i++){
    Z.col(i) = (Z.col(i) - moments(i+1+px,0))/ moments(i+1+px,1);
  }
}
