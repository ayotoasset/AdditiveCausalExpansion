// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include "ace_kernelutilities.hpp"

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List kernmat_SE_cpp(const arma::mat& X1,
                          const arma::mat& X2,
                          const arma::mat& Z1,
                          const arma::mat& Z2,
                          const arma::vec& parameters) {
  //input X1,X2,Z1,Z2,parameters
  unsigned int n1 = X1.n_rows;
  unsigned int n2 = X2.n_rows;
  unsigned int p = X2.n_cols;
  unsigned int B = Z1.n_cols + 1; //including nuisance term

  //(help) storage variables
  arma::mat Kfull(n1, n2);
  arma::rowvec tmprow(n2);
  arma::cube tmpX(n1, n2, B); tmpX.zeros(); //reuse for all additive elements
  // storage cost for training set: O(n^2*B)

  //create the Euclidian distance metric:
  for(unsigned int i = 0; i < p; i++){ // for every element of x
    for(unsigned int r = 0; r < n1; r++){ //for every output row
      tmprow = arma::pow(X1(r, i) - conv_to<arma::rowvec>::from(X2.col(i)), 2);

      for(unsigned int b = 0; b < B; b++){ //3rd dimension of array is called slice
        tmpX.slice(b).row(r) += tmprow * exp( - parameters[1+b+B*(i+1)]); //L(i,b)
      }
    }
  }

  //nuisance term does not interact with Z - simple copy
  tmpX.slice(0) = exp(parameters[2] - tmpX.slice(0)); //lambda[0]
  Kfull = tmpX.slice(0);
  //For all basis expanded terms, add interaction with Z and then add to Kfull
  for(unsigned int b = 1; b < B; b++){
    for(unsigned int r = 0; r < n1; r++){
      if (Z1(r,b-1) == 0 ) {
        tmpX.slice(b).row(r).zeros();
        continue;
      }
      for(unsigned int c = 0; c < n2; c++){
        if (Z2(c,b-1) == 0 ) {
          tmpX(r,c,b) = 0;
          //continue;
        } else {
          tmpX(r,c,b) = (sign(Z1(r, b - 1)) * sign(Z2(c, b - 1))) * exp(parameters[2 + b] - tmpX(r, c, b) + log(std::abs(Z1(r, b - 1))) + log(std::abs(Z2(c, b - 1))));
        }

        //tmpX(r,c,b) = exp(parameters[2+b] - tmpX(r,c,b)) * Z1(r,b-1) * Z2(c,b-1); //lambda[b]
        //tmpX(r,c,b) = exp(parameters[2+b] - tmpX(r,c,b) - pow(Z1(r,b-1) - Z2(c,b-1),2) ); // replace each element of tmpX
      }
    }
    Kfull += tmpX.slice(b);
  }
  return Rcpp::List::create(_("full") =  Kfull,
                            _("elements") = tmpX);
}

// [[Rcpp::export]]
Rcpp::List kernmat_SE_symmetric_cpp(const arma::mat& X,
                                    const arma::mat& Z,
                                    const arma::vec& parameters) {
  //input X1,X2,Z1,Z2,parameters
  unsigned int n = X.n_rows;
  unsigned int p = X.n_cols;
  unsigned int B = Z.n_cols + 1; //including nuisance term
  unsigned int cnt;

  //(help) storage variables
  double tmp = 0;
  arma::mat tmpX(n*(n+1)/2, B); tmpX.zeros(); //elements of upper/lower triangle
  arma::mat Kfull(n, n);
  arma::cube Ks(n, n, B); Ks.zeros();

  for(unsigned int i = 0; i < p; i++){ // for every element of x
    cnt = 0;
    for(unsigned int r = 0; r < n; r++){ //for every output row
      for(unsigned int c = r; c < n; c++){
        tmp = pow(X(r, i) - X(c, i),2);
        for(unsigned int b = 0; b < B; b++){
          tmpX(cnt, b) += tmp * exp( - parameters[1 + b + B * (i+1)]); //L(i,b)
        }
        cnt++;
      }
    }
  }

  tmpX.col(0) = exp(parameters[2] - tmpX.col(0));
  Ks.slice(0) = uppertri2symmat(tmpX.col(0), n);
  for(unsigned int b = 1; b < B; b++){
    cnt=0;
    for(unsigned int r = 0; r < n; r++){
      //skip row if Z is zero, important for binary variables
      if (Z(r, b - 1) == 0) {
        tmpX.submat(cnt, b, cnt + (n-r) - 1, b).zeros();
        // need to increase counter by the number of upper-triangle elements we skip
        cnt += (n - r);
        continue;
      }
      for(unsigned int c = r; c < n; c++){
        // replace each element of tmpX
        if (Z(c, b - 1) == 0) {
          tmpX(cnt, b) = 0;
          cnt++;
          continue;
        }
        if (Z(r, b - 1) == 0){
          tmpX(cnt, b) = 0;
        } else {
          tmpX(cnt, b) = (sign(Z(r, b - 1)) * sign(Z(c, b - 1))) * exp(parameters[2 + b] - tmpX(cnt, b) + log(std::abs(Z(r, b - 1))) + log(std::abs(Z(c, b - 1))));
        }
        cnt++;
      }
    }
    // being done with the elements in tmpX we use it to construct the kernel matrices in an efficient way
    Ks.slice(b) = uppertri2symmat( tmpX.col(b), n );
    for(unsigned int j = 0; j < (n * (n + 1) / 2); j++) {
      tmpX(j, 0) += tmpX(j, b);
      }
  }
  Kfull = uppertri2symmat(tmpX.col(0), n); //sum sparse vectorization instead of matrices

  return Rcpp::List::create(_("full") =  Kfull,
                            _("elements") = Ks);
}

// [[Rcpp::export]]
Rcpp::List invkernel_cpp(arma::mat pdmat,
                         const double& sigma){
  unsigned int n = pdmat.n_cols;
  arma::vec eigval(n); eigval.ones();
  arma::mat eigvec(n, n); eigvec.eye();
  pdmat.diag() += exp(sigma);

  if(!arma::eig_sym(eigval, eigvec, pdmat)){
    Rcout << "Eigenvalue decomp. not completed." << std::endl;
  }

  //get inverse and eigenvalues
  arma::mat invKmat(n, n);

  invKmat = eigvec;
  for(unsigned int i = 0; i < n; i++){
    invKmat.col(i) = invKmat.col(i) / eigval[i];
  }
  invKmat = invKmat * eigvec.t();

  return Rcpp::List::create(_("eigenval") =  eigval,
                            _("inv") = invKmat);
}

///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// gradients

inline arma::mat evid_scale_gradients(const arma::mat& X,
                               const arma::mat& Kaa,
                               const arma::cube& K,
                               arma::vec L,
                               unsigned int B){//,
                               //double length_scale_hyperparameter = 0){
  // produce a (p x B) matrix "grad" with the gradients of L
  unsigned int n = X.n_rows;
  unsigned int p = X.n_cols;
  arma::mat tmpX(n, n);

  for(unsigned int i=0; i < p; i++){
    for(unsigned int r=0; r < n; r++){
      tmpX.col(r) = pow(X(r, i) - X.col(i), 2);
    }
    for(unsigned int b=0; b < B; b++){
      // replace parameter with its gradient
      L[b + B*i] = evid_grad(Kaa, K.slice(b) % tmpX * exp( - L[b + B * i]));
                   // - length_scale_hyperparameter * sqrt(exp(L[b + B * i]));
                   // - 2 * exp(L[b + B * i]) * length_scale_hyperparameter
                   // / (1 + L[b + B*i] * length_scale_hyperparameter);
        //- 0.5 * arma::trace( Kaa * (K.slice(b) % tmpX)) * exp( - L[b + B*i] );
    }
  }
  return L;
}

// reduced to only an output list due to specification of the optimizers
// [[Rcpp::export]]
arma::vec grad_SE_cpp(const arma::vec& y,
                      const arma::mat& X,
                      const arma::mat& Z,
                      arma::mat& Kfull,
                      arma::cube& K,
                      arma::mat& invKmatn,
                      arma::vec& eigenval,
                      const arma::vec& parameters,
                      arma::vec& stats,
                      const unsigned int& B) {
  // set constants
  unsigned int n = X.n_rows;
  unsigned int px = X.n_cols;

  // initialize variables
  arma::vec gradients(parameters.size());
  gradients.zeros();

  arma::colvec ybar(n);
  ybar = y - parameters[1];

  arma::colvec alpha(n);
  alpha = invKmatn * ybar;

  arma::mat tmpK(n, n);
  tmpK = invKmatn - alpha * alpha.t();

  // sigma
  gradients[0] = sigma_gradient(tmpK, parameters[0]);

  // kernel-scale "lambda"
  for(unsigned int b=0; b<B; b++){
    gradients[2+b] = 0.1 * evid_grad(tmpK, K.slice(b));
  }

  // L
  gradients.rows(2 + B, 1 + B * (px + 1)) = evid_scale_gradients(X,
                 tmpK, K, parameters.rows(2 + B, 1 + B * (px + 1)), B);

  // mu - gradient approach
  gradients[1] = arma::sum(invKmatn * ybar);

  // stats is returned as a reference
  //RMSE
  stats(0) = pow(arma::norm(ybar - (Kfull * alpha)), 2);
  //Evidence
  stats(1) = logevidence(y, alpha, eigenval, n);

  return gradients;
}
