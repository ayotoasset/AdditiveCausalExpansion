// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <cmath>

using namespace arma;
using namespace Rcpp;

//more efficient with pointers and a class but well

arma::mat uppertri2symmat_noadd(arma::vec matvec,unsigned int dim){
  arma::mat out(dim,dim);
  unsigned int cnt = 0;
  for(unsigned int r = 0; r < dim; r++ ){
    for(unsigned int c = r; c < dim; c++ ){
      out(r,c) = out(c,r) = matvec(cnt);
      cnt++;
    }
  }
  return out;
}


// [[Rcpp::export]]
Rcpp::List kernmat_SE_noadd_cpp(const arma::mat& X1,const arma::mat& X2,
                          const arma::vec& parameters) {
  //input X1,X2,Z1,Z2,parameters
  unsigned int n1 = X1.n_rows;
  unsigned int n2 = X2.n_rows;
  unsigned int p = X2.n_cols;
  unsigned int B = 1;

  //(help) storage variables
  //arma::mat Kfull(n1, n2);
  arma::rowvec tmprow(n2);
  arma::cube tmpX(n1,n2,B); tmpX.zeros(); //reuse for all additive elements
  // storage cost for training set: O(n^2*B)

  for(unsigned int i = 0; i < p; i++){ // for every element of x
    for(unsigned int r = 0; r < n1; r++){ //for every output row
      tmprow = arma::pow(X1(r,i) - conv_to<rowvec>::from(X2.col(i)), 2);

      for(unsigned int b = 0; b < B; b++){ //3rd dimension of array is called slice
        tmpX.slice(b).row(r) += tmprow * exp(- parameters[1+b+B*(i+1)]  );//L(i,b)
      }
    }
  }

  //outest-most iteration is basis dimension
  tmpX.slice(0) = exp(parameters[2] - tmpX.slice(0)); //lambda[0]

  return Rcpp::List::create(_("full") =  tmpX.slice(0),
                            _("elements") = tmpX);
}


// [[Rcpp::export]]
Rcpp::List kernmat_SE_symmetric_noadd_cpp(const arma::mat& X, const arma::vec& parameters) {
  //input X1,X2,Z1,Z2,parameters
  unsigned int n = X.n_rows;
  unsigned int p = X.n_cols;
  unsigned int B = 1; //including nuisance term
  unsigned int cnt;

  //(help) storage variables
  double tmp = 0;
  arma::mat tmpX(n*(n+1)/2,B); tmpX.zeros();
  arma::mat Kfull(n,n);// Kfull.zeros();
  arma::cube Ks(n,n,B); Ks.zeros();

  for(unsigned int i = 0; i < p; i++){ // for every element of x
    cnt = 0;
    for(unsigned int r = 0; r < n; r++){ //for every output row
      for(unsigned int c = r; c < n; c++){
        tmp = pow(X(r,i) - X(c,i),2);
        for(unsigned int b = 0; b < B; b++){ //3rd dimension of array is called slice
          tmpX(cnt,b) += tmp * exp(- parameters[1+b+B*(i+1)] ); //L(i,b)
        } cnt++;
      }
    }
  }

  tmpX.col(0) = exp(parameters[2] - tmpX.col(0) ); Ks.slice(0) = uppertri2symmat_noadd( tmpX.col(0), n);
  return Rcpp::List::create(_("full") =  Ks.slice(0),
                            _("elements") = Ks);
}

// [[Rcpp::export]]
arma::mat kernmat_marginal_noadd_cpp(const arma::mat& Kmat,const arma::mat& Z1,const arma::mat& Z2,
                                const arma::vec& parameters) {

  unsigned int n1 = Z1.n_rows;
  unsigned int n2 = Z2.n_rows;
  arma::mat tmpX(n1,n2); tmpX.zeros();
  arma::rowvec tmprow(n2);

  for(unsigned int r = 0; r < n1; r++){ //for every output row
      tmprow = - 2 * (Z1(r,0) - conv_to<rowvec>::from(Z2.col(0))) * exp(- parameters[parameters.size()]) ;
      tmpX.row(r) = tmprow;
  }

  return Kmat % tmpX;
}


