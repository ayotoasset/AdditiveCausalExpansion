// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include "ace_kernelutilities.h"

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List kernmat_Matern52_cpp(const arma::mat& X1,const arma::mat& X2,const arma::mat& Z1,const arma::mat& Z2,
                          const arma::vec& parameters) {
  //input X1,X2,Z1,Z2,parameters
  unsigned int n1 = X1.n_rows;
  unsigned int n2 = X2.n_rows;
  unsigned int p = X2.n_cols;
  unsigned int B = Z1.n_cols + 1; //including nuisance term

  //(help) storage variables
  arma::mat Kfull(n1, n2);
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
  tmpX.slice(0) = (1 + arma::sqrt(5*tmpX.slice(0)) + 5*tmpX.slice(0)/3 ) % exp(parameters[2] - sqrt(5*tmpX.slice(0)) ); //lambda[0]
  Kfull = tmpX.slice(0);
  for(unsigned int b = 1; b < B; b++){
    for(unsigned int r = 0; r < n1; r++){
      if (Z1(r,b-1)==0 ) {tmpX.slice(b).row(r).fill(0); continue;}
      for(unsigned int c = 0; c < n2; c++){
        if (Z2(c,b-1)==0 ) {tmpX(r,c,b)=0; continue;}
        tmpX(r,c,b) = (1 + sqrt(5*tmpX(r,c,b)) + 5*tmpX(r,c,b)/3 ) * exp(parameters[2+b] - sqrt(5*tmpX(r,c,b))) * Z1(r,b-1) * Z2(c,b-1); //lambda[b]
      }
    }
    Kfull += tmpX.slice(b);
  }
  return Rcpp::List::create(_("full") =  Kfull,
                            _("elements") = tmpX);
}

// [[Rcpp::export]]
Rcpp::List kernmat_Matern32_cpp(const arma::mat& X1,const arma::mat& X2,const arma::mat& Z1,const arma::mat& Z2,
                                const arma::vec& parameters) {
  //input X1,X2,Z1,Z2,parameters
  unsigned int n1 = X1.n_rows;
  unsigned int n2 = X2.n_rows;
  unsigned int p = X2.n_cols;
  unsigned int B = Z1.n_cols + 1; //including nuisance term

  //(help) storage variables
  arma::mat Kfull(n1, n2);
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
  tmpX = arma::sqrt(tmpX);

  //outest-most iteration is basis dimension
  tmpX.slice(0) = (1 + sqrt(3)*tmpX.slice(0)) % exp(parameters[2] - sqrt(3)*tmpX.slice(0) ); //lambda[0]
  Kfull = tmpX.slice(0);
  for(unsigned int b = 1; b < B; b++){
    for(unsigned int r = 0; r < n1; r++){
      if (Z1(r,b-1)==0 ) {tmpX.slice(b).row(r).fill(0); continue;}
      for(unsigned int c = 0; c < n2; c++){
        if (Z2(c,b-1)==0 ) {tmpX(r,c,b)=0; continue;}
        tmpX(r,c,b) = (1 + sqrt(3)*tmpX(r,c,b) ) * exp(parameters[2+b] - sqrt(3)*tmpX(r,c,b)) * Z1(r,b-1) * Z2(c,b-1); //lambda[b]
      }
    }
    Kfull += tmpX.slice(b);
  }
  return Rcpp::List::create(_("full") =  Kfull,
                            _("elements") = tmpX);
}

// [[Rcpp::export]]
Rcpp::List kernmat_Matern12_cpp(const arma::mat& X1,const arma::mat& X2,const arma::mat& Z1,const arma::mat& Z2,
                                const arma::vec& parameters) {
  //input X1,X2,Z1,Z2,parameters
  unsigned int n1 = X1.n_rows;
  unsigned int n2 = X2.n_rows;
  unsigned int p = X2.n_cols;
  unsigned int B = Z1.n_cols + 1; //including nuisance term

  //(help) storage variables
  arma::mat Kfull(n1, n2);
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
  tmpX = arma::sqrt(tmpX);

  //outest-most iteration is basis dimension
  tmpX.slice(0) = exp(parameters[2] - tmpX.slice(0) ); //lambda[0]
  Kfull = tmpX.slice(0);
  for(unsigned int b = 1; b < B; b++){
    for(unsigned int r = 0; r < n1; r++){
      if (Z1(r,b-1)==0 ) {tmpX.slice(b).row(r).fill(0); continue;}
      for(unsigned int c = 0; c < n2; c++){
        if (Z2(c,b-1)==0 ) {tmpX(r,c,b)=0; continue;}
        tmpX(r,c,b) = exp(parameters[2+b] - tmpX(r,c,b)) * Z1(r,b-1) * Z2(c,b-1); //lambda[b]
      }
    }
    Kfull += tmpX.slice(b);
  }
  return Rcpp::List::create(_("full") =  Kfull,
                            _("elements") = tmpX);
}


// [[Rcpp::export]]
Rcpp::List kernmat_Matern52_symmetric_cpp(const arma::mat& X, const arma::mat& Z, const arma::vec& parameters) {
  //input X1,X2,Z1,Z2,parameters
  unsigned int n = X.n_rows;
  unsigned int p = X.n_cols;
  unsigned int B = Z.n_cols + 1; //including nuisance term
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

  tmpX.col(0) = (1 + sqrt(5*tmpX.col(0)) + 5*tmpX.col(0)/3 ) % exp(parameters[2] - sqrt(5*tmpX.col(0))); //lambda[b]
  Ks.slice(0) = uppertri2symmat( tmpX.col(0), n);
  for(unsigned int b = 1; b < B; b++){
    cnt=0;
    for(unsigned int r = 0; r < n; r++){
      if (Z(r,b-1)==0 ) { // need to increase counter by the number of upper-triangle elements we skip
        tmpX.submat(cnt,b,cnt+(n-r)-1,b).fill(0); cnt += n-r;
        continue; }
      for(unsigned int c = r; c < n; c++){
        tmpX(cnt,b) = (1 + sqrt(5*tmpX(cnt,b)) + 5*tmpX(cnt,b)/3 ) * exp(parameters[2+b] - sqrt(5*tmpX(cnt,b))) * Z(r,b-1) * Z(c,b-1); //lambda[b]
        cnt++;
    }
    }
    // being done with the elements in tmpX we use it to construct the kernel matrices in an efficient way
    Ks.slice(b) = uppertri2symmat(tmpX.col(b), n);
    for(unsigned int j = 0; j < n*(n+1)/2; j++){ tmpX(j,0) += tmpX(j,b); }
  }
  Kfull = uppertri2symmat(tmpX.col(0), n); //sum sparse vectorization instead of matrices

  //Rcpp::Rcout << "Kele(0,0,1)" << Ks.slice(1)(0,0) << std::endl;
  return Rcpp::List::create(_("full") =  Kfull,
                            _("elements") = Ks);
}

// [[Rcpp::export]]
Rcpp::List kernmat_Matern32_symmetric_cpp(const arma::mat& X, const arma::mat& Z, const arma::vec& parameters) {
  //input X1,X2,Z1,Z2,parameters
  unsigned int n = X.n_rows;
  unsigned int p = X.n_cols;
  unsigned int B = Z.n_cols + 1; //including nuisance term
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
  tmpX = arma::sqrt(tmpX);

  tmpX.col(0) = (1 + sqrt(3)*tmpX.col(0) ) % exp(parameters[2] - sqrt(3)*tmpX.col(0) ); //lambda[b]
  Ks.slice(0) = uppertri2symmat(tmpX.col(0), n);
  for(unsigned int b = 1; b < B; b++){
    cnt=0;
    for(unsigned int r = 0; r < n; r++){
      if (Z(r,b-1)==0 ) { // need to increase counter by the number of upper-triangle elements we skip
        tmpX.submat(cnt,b,cnt+(n-r)-1,b).fill(0); cnt += n-r;
        continue; }
      for(unsigned int c = r; c < n; c++){
        tmpX(cnt,b) = (1 + sqrt(3)*tmpX(cnt,b)) * exp(parameters[2+b] - sqrt(3)*tmpX(cnt,b)) * Z(r,b-1) * Z(c,b-1); //lambda[b]
        cnt++;
      }
    }
    // being done with the elements in tmpX we use it to construct the kernel matrices in an efficient way
    Ks.slice(b) = uppertri2symmat( tmpX.col(b), n );
    for(unsigned int j = 0; j < n*(n+1)/2; j++){ tmpX(j,0) += tmpX(j,b); }
  }
  Kfull = uppertri2symmat(tmpX.col(0), n); //sum sparse vectorization instead of matrices

  //Rcpp::Rcout << "Kele(0,0,1)" << Ks.slice(1)(0,0) << std::endl;
  return Rcpp::List::create(_("full") =  Kfull,
                            _("elements") = Ks);
}

// [[Rcpp::export]]
Rcpp::List kernmat_Matern12_symmetric_cpp(const arma::mat& X, const arma::mat& Z, const arma::vec& parameters) {
  //input X1,X2,Z1,Z2,parameters
  unsigned int n = X.n_rows;
  unsigned int p = X.n_cols;
  unsigned int B = Z.n_cols + 1; //including nuisance term
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
  tmpX = arma::sqrt(tmpX);

  tmpX.col(0) = exp(parameters[2] - tmpX.col(0) ); //lambda[b]
  Ks.slice(0) = uppertri2symmat( tmpX.col(0), n);
  for(unsigned int b = 1; b < B; b++){
    cnt=0;
    for(unsigned int r = 0; r < n; r++){
      if (Z(r,b-1)==0 ) { // need to increase counter by the number of upper-triangle elements we skip
        tmpX.submat(cnt,b,cnt+(n-r)-1,b).fill(0); cnt += n-r;
        continue; }
      for(unsigned int c = r; c < n; c++){
        tmpX(cnt,b) = exp(parameters[2+b] - tmpX(cnt,b) ) * Z(r,b-1) * Z(c,b-1); //lambda[b]
        cnt++;
      }
    }
    // being done with the elements in tmpX we use it to construct the kernel matrices in an efficient way
    Ks.slice(b) = uppertri2symmat( tmpX.col(b), n );
    for(unsigned int j = 0; j < n*(n+1)/2; j++){ tmpX(j,0) += tmpX(j,b); }
  }
  Kfull = uppertri2symmat(tmpX.col(0), n); //sum sparse vectorization instead of matrices

  //Rcpp::Rcout << "Kele(0,0,1)" << Ks.slice(1)(0,0) << std::endl;
  return Rcpp::List::create(_("full") =  Kfull,
                            _("elements") = Ks);
}


///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// gradients

double evid_Matern_grad(const arma::mat& Kaa, const arma::mat& dK) {
  return - 0.5 * arma::trace( Kaa * dK);
}

arma::mat evid_scale_Matern52_gradients(const arma::mat& X,const arma::mat& Z, const arma::mat& Kaa, const arma::cube& K, arma::vec L,const arma::vec& lambda, unsigned int B){
  // produce a (p x B) matrix "grad" with the gradients of L
  unsigned int n = X.n_rows;
  unsigned int p = X.n_cols;

  arma::rowvec tmprow(n);
  arma::mat tmpX2(n,n);
  arma::cube tmpX(n,n,B); tmpX.zeros(); //reuse for all additive elements

  for(unsigned int i = 0; i < p; i++){ // for every element of x
    for(unsigned int r = 0; r < n; r++){ //for every output row
      tmprow = arma::pow(X(r,i) - conv_to<rowvec>::from(X.col(i)), 2);
      for(unsigned int b = 0; b < B; b++){ //3rd dimension of array is called slice
        tmpX.slice(b).row(r) += tmprow * exp(- L[b + B*i]  );//L(i,b)
      }
    }
  }

  arma::mat tmpZ(n,n); tmpZ.ones();
  for(unsigned int b = 0; b < B; b++){
    if(b>0){ //for b=0 take intitial values
      tmpZ = Z.col(b-1) * Z.col(b-1).t();
    }
    tmpX2 = sqrt(5*tmpX.slice(b)) + 5*tmpX.slice(b)/3;
    tmpX.slice(b) =  0.5 * K.slice(b) % tmpX2 / ( (1+tmpX2) % sqrt(tmpX.slice(b)) );
    //as tmpX is at least 0 on the diagonal and he division gives nan
    tmpX.slice(b).elem(find(tmpX2==0)).zeros();
    //Rcpp::Rcout << tmpZ(0,0) <<" | "<< tmpZ(0,1) <<" | "<< tmpZ(0,2) << std::endl;
    tmpX.slice(b) -= lambda(b) *sqrt(5)/6 * tmpZ;
  }

  for(unsigned int i=0; i < p; i++){
    for(unsigned int r=0; r < n; r++){
      tmpX2.col(r) = pow(X(r, i) - X.col(i),2);
    }
    for(unsigned int b = 0; b < B; b++){
      L[b + B*i] = - 0.5 * sqrt(5) * arma::trace( Kaa * ( tmpX.slice(b) % tmpX2 ) ) * exp( - L[b + B*i] );
    }
  }
  return L;
}

arma::mat evid_scale_Matern32_gradients(const arma::mat& X, const arma::mat& Kaa, const arma::cube& K, arma::vec L, unsigned int B){
  // produce a (p x B) matrix "grad" with the gradients of L
  unsigned int n = X.n_rows;
  unsigned int p = X.n_cols;

  //Rcpp::Rcout << "Matern 32 gradients" << std::endl;

  arma::rowvec tmprow(n);
  arma::mat tmpX2(n,n);
  arma::cube tmpX(n,n,B); tmpX.zeros(); //reuse for all additive elements

  for(unsigned int i = 0; i < p; i++){ // for every element of x
    for(unsigned int r = 0; r < n; r++){ //for every output row
      tmprow = arma::pow(X(r,i) - conv_to<rowvec>::from(X.col(i)), 2);
      for(unsigned int b = 0; b < B; b++){ //3rd dimension of array is called slice
        tmpX.slice(b).row(r) += tmprow * exp(- L[b + B*i]  );//L(i,b)
      }
    }
  }

  //precompute matrices used in higher order loops
  for(unsigned int b = 0; b < B; b++){
    tmpX.slice(b) =  K.slice(b) / ( 1+sqrt(3*tmpX.slice(b)) );
  }

  //tmpX.slice(b) = 0;
  for(unsigned int i=0; i < p; i++){
    for(unsigned int r=0; r < n; r++){
      tmpX2.col(r) = pow(X(r, i) - X.col(i),2);
    }
    for(unsigned int b = 0; b < B; b++){
      L[b + B*i] = - 0.25 * 9 * arma::trace( Kaa * (tmpX.slice(b) % tmpX2) ) * exp( - L[b + B*i] );
    }
  }
  return L;
}


arma::mat evid_scale_Matern12_gradients(const arma::mat& X, const arma::mat& Kaa, const arma::cube& K,arma::vec L, unsigned int B){
  // produce a (p x B) matrix "grad" with the gradients of L
  unsigned int n = X.n_rows;
  unsigned int p = X.n_cols;

  //Rcpp::Rcout << "Matern 12 gradients" << std::endl;

  arma::rowvec tmprow(n);
  arma::mat tmpX2(n,n);
  arma::cube tmpX(n,n,B); tmpX.zeros(); //reuse for all additive elements

  for(unsigned int i = 0; i < p; i++){ // for every element of x
    for(unsigned int r = 0; r < n; r++){ //for every output row
      tmprow = arma::pow(X(r,i) - conv_to<rowvec>::from(X.col(i)), 2);
      for(unsigned int b = 0; b < B; b++){ //3rd dimension of array is called slice
        tmpX.slice(b).row(r) += tmprow * exp(- L[b + B*i]  );//L(i,b)
      }
    }
  }
  //tmpX =

  for(unsigned int b=0; b < B; b++){
    tmpX.slice(b) =  K.slice(b) / arma::sqrt(tmpX.slice(b));
  }

  for(unsigned int i=0; i < p; i++){
    for(unsigned int r=0; r < n; r++){
      tmpX2.col(r) = pow(X(r, i) - X.col(i),2);
    }
    for(unsigned int b=0; b < B; b++){
      //tmpX.slice(b) = K.slice(b) % tmpX.slice(b);
    L[b + B*i] = - 0.25 * arma::trace( Kaa * (tmpX.slice(b) % tmpX2) ) * exp( - L[b + B*i] );
    }
  }
  return L;
}

// reduced to only an output list due to specification of the optimizers
// [[Rcpp::export]]
arma::vec grad_Matern_cpp(const arma::vec& y, const arma::mat& X, const arma::mat& Z,
                            arma::mat& Kfull, arma::cube& K, arma::mat& invKmatn, arma::vec& eigenval,
                            const arma::vec& parameters,
                            arma::vec& stats, const unsigned int& B,unsigned int rho) {

  unsigned int n = X.n_rows;
  unsigned int px = X.n_cols;
  //unsigned int B = Z.n_cols+1;

  arma::vec gradients(parameters.size());

  arma::colvec ybar(n);
  arma::colvec alpha(n);
  arma::mat tmpK(n,n);

  ybar = y - parameters[1];
  alpha = invKmatn * ybar;
  tmpK = invKmatn - alpha * alpha.t();

  //sigma
  gradients[0] = - 0.5 * arma::trace( tmpK ) * exp( parameters[0] ); // thus avoid allocating dK

  //lambda
  for(unsigned int b=0;b<B;b++){
    gradients[2+b] = evid_Matern_grad(tmpK, K.slice(b));
  }

  //L
  switch(rho){
    case 0: gradients.rows(2+B,1+B*(px+1)) = evid_scale_Matern12_gradients(X, tmpK, K, parameters.rows(2+B,1+B*(px+1)), B);
      break;
    case 1: gradients.rows(2+B,1+B*(px+1)) = evid_scale_Matern32_gradients(X, tmpK, K, parameters.rows(2+B,1+B*(px+1)), B);
      break;
    case 2: gradients.rows(2+B,1+B*(px+1)) = evid_scale_Matern52_gradients(X, Z, tmpK, K, parameters.rows(2+B,1+B*(px+1)), parameters.rows(2,1+B), B);
  }

  //mu - gradient approach
  gradients[1]=0; //arma::sum( invKmatn * ybar )

  //RMSE
  stats(0) = pow(arma::norm(ybar - (Kfull * alpha)),2);
  //Evidence
  stats(1) = - 0.5 * (n * log( 2.0 * arma::datum::pi ) + sum(log(eigenval)) + arma::dot( ybar, alpha ) ) ;
  //stats is returned as a reference

  //clip gradients
  //umat idx=(gradients>20);
  /*for(unsigned int g=0; g<gradients.size();g++){
    if(gradients(g)>20){
      gradients.row(g) = 20;
    } else if(gradients(g)<-20) {
      gradients.row(g) = -20;
    }
  }*/
  return gradients;
}
