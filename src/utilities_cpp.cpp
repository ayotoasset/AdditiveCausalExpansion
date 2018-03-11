// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
double mu_solution_cpp(arma::colvec& y,
                       arma::mat& invKmat) {
  // calculates the exact solutions to the maximization problem
  return  0.5 * arma::sum(invKmat * y) / arma::accu(invKmat);
}

// [[Rcpp::export]]
arma::mat normalize_train(arma::vec& y, arma::mat& X, arma::mat& Z) {
  unsigned int px = X.n_cols;
  unsigned int pz = Z.n_cols;

  arma::mat moments(1 + px + pz, 3);
  moments.col(0).zeros();
  moments.col(1).ones();
  moments.col(2).zeros();

  //create index of binary variables in X and Z
  arma::uvec isbinary(px + pz); isbinary.zeros();
  arma::vec tmp;
  // X
  for(unsigned int i = 0; i < px; i++){
    tmp = arma::unique(X.col(i));
    if(tmp.n_elem == 2) {
      isbinary(i) = moments(i + 1, 2) = 1;
      if(arma::min(tmp) != 0) {
        //set smaller value to zero
        moments(i, 0) = arma::min(tmp);
      }
      if(arma::max(tmp) != 1) {
        //set larger value to one
        moments(i, 1) = (arma::max(tmp)- arma::min(tmp));
      }
      // locate and scale:
      X.col(i) -= moments(i, 0);
      X.col(i) /= moments(i, 1);

    } else if(tmp.n_elem == 1) {
      Rcout << "Column " << i + 1 << " of X is constant." << std::endl;
      X.col(i).zeros();
    }
  }
  // Z
  for(unsigned int i = px; i < (px + pz); i++){
    tmp = arma::unique(Z.col(i - px));
    if(tmp.n_elem == 2) {
      isbinary(i) = moments(i + 1, 2) = 1;
      if(arma::min(tmp) != 0) {
        //set smaller value to zero
        moments(i, 0) = arma::min(tmp);
      }
      if(arma::max(tmp) != 1) {
        //set larger value to one
        moments(i, 1) = (arma::max(tmp)- arma::min(tmp));
      }
      // locate and scale:
      Z.col(i - px) -= moments(i, 0);
      Z.col(i - px) /= moments(i, 1);
    } else if(tmp.n_elem == 1){
      Rcout << "Column " << i << " of Z is constant." << std::endl;
      Z.col(i).zeros();
    }
  }

  //get mean and de-mean
  moments(0, 0) = mean(y);
  y = y - moments(0, 0);
  for(unsigned int i = 1; i < (px + 1); i++){
    if(isbinary(i - 1) == 0){
      moments(i, 0) = median(X.col(i - 1));
      X.col(i - 1) -= moments(i, 0);
    }
  }
  for(unsigned int i = px + 1; i < (px + pz + 1); i++){
    if(isbinary(i - 1) == 0){
      moments(i, 0) = median(Z.col(i - px - 1));
      Z.col(i - px - 1) -= moments(i, 0);
    }
  }
  //rescale
  moments(0, 1) = arma::stddev(y); // mean is a parameter
  y = y / moments(0, 1);
  for(unsigned int i = 1; i < (px + 1); i++){
    if(isbinary(i - 1) == 0){
      moments(i, 1) = arma::max(abs(X.col(i - 1)));
      X.col(i - 1) /= moments(i, 1);
    } /*else {
      moments(i, 1) = 1;//arma::stddev(X.col(i-1));
      X.col(i - 1) /=  moments(i, 1);
    }*/
  }
  //important for Z to be within unit circle for splines
  for(unsigned int i = px + 1; i < (px + pz + 1); i++){
    if(isbinary(i - px - 1) == 0){
      moments(i, 1) = arma::max(abs(Z.col(i - px - 1)));
      Z.col(i - px - 1) = Z.col(i - px - 1) / moments(i, 1);
    }
  }
  return moments;
}


// [[Rcpp::export]]
void normalize_test(arma::mat& X, arma::mat& Z, const arma::mat& moments) {
  unsigned int px = X.n_cols;
  unsigned int pz = Z.n_cols;

    for(unsigned int i = 0; i < px; i++){
    X.col(i) = (X.col(i) - moments(i + 1, 0)) / moments(i + 1, 1);
  }
  for(unsigned int i = 0; i < pz; i++){
    Z.col(i) = (Z.col(i) - moments(i + 1 + px, 0)) / moments(i + 1 + px, 1);
  }
}

// [[Rcpp::export]]
void norm_clip_cpp(bool flag, arma::vec& grads, double max_length){
  // cut gradients by reference if flag true and norm above max_length
  if(flag) {
    double L2_norm = arma::norm(grads);
    if ((L2_norm > max_length) &  arma::is_finite(L2_norm) & (L2_norm != 0)) {
      grads = grads / L2_norm;
    }
  }
}

