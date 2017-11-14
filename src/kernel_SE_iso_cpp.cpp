// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <cmath>

using namespace arma;
using namespace Rcpp;

arma::mat uppertri2symmat_iso(arma::vec matvec,unsigned int dim){
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
Rcpp::List kernmat_SE_iso_cpp(const arma::mat& X1,const arma::mat& X2,const arma::mat& Z1,const arma::mat& Z2,
                              const arma::vec& parameters) {
  //input X1,X2,Z1,Z2,parameters
  unsigned int n1 = X1.n_rows;
  unsigned int n2 = X2.n_rows;
  unsigned int p = X2.n_cols;
  unsigned int B = 2;//Z1.n_cols + 1; //including nuisance term

  //(help) storage variables
  arma::mat Kfull(n1, n2);
  arma::rowvec tmprow(n2);
  arma::cube tmpX(n1,n2,B); tmpX.zeros(); //reuse for all additive elements
  // storage cost for training set: O(n^2*B)

  for(unsigned int i = 0; i < p; i++){ // for every element of x
    for(unsigned int r = 0; r < n1; r++){ //for every output row
      tmprow = arma::pow(X1(r,i) - conv_to<rowvec>::from(X2.col(i)), 2);

      tmpX.slice(0).row(r) += tmprow * exp(- parameters[1+0+B*(i+1)]  );
      tmpX.slice(1).row(r) += tmprow * exp(- parameters[1+1+B*(i+1)]  );
    }
  }

  //outest-most iteration is basis dimension
  tmpX.slice(0) = exp(parameters[2] - tmpX.slice(0)); //lambda[0]
  Kfull = tmpX.slice(0);
  unsigned int b = 1; //instead of the loop
  for(unsigned int r = 0; r < n1; r++){
    if (Z1(r,b-1)==0 ) {tmpX.slice(b).row(r).fill(0); continue;}
    for(unsigned int c = 0; c < n2; c++){
      if (Z2(c,b-1)==0 ) {tmpX(r,c,b)=0; continue;}
      tmpX(r,c,b) = exp(parameters[2+b] - tmpX(r,c,b)) * arma::dot(Z1.row(r),Z2.row(c)); //lambda[b]
    }
  }
  Kfull += tmpX.slice(b);

  return Rcpp::List::create(_("full") =  Kfull,
                            _("elements") = tmpX);
}


// [[Rcpp::export]]
Rcpp::List kernmat_SE_iso_symmetric_cpp(const arma::mat& X, const arma::mat& Z, const arma::vec& parameters) {
  //input X1,X2,Z1,Z2,parameters
  unsigned int n = X.n_rows;
  unsigned int p = X.n_cols;
  unsigned int B = 2;//Z.n_cols + 1; //including nuisance term
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
        tmpX(cnt,0) += tmp * exp(- parameters[1+0+B*(i+1)] ); //L(i,b)
        tmpX(cnt,1) += tmp * exp(- parameters[1+1+B*(i+1)] ); //L(i,b)
        //for(unsigned int b = 0; b < B; b++){ //3rd dimension of array is called slice
        //  tmpX(cnt,b) += tmp * exp(- parameters[1+b+B*(i+1)] ); //L(i,b)
        //}
        cnt++;
      }
    }
  }

  tmpX.col(0) = exp(parameters[2] - tmpX.col(0) ); Ks.slice(0) = uppertri2symmat_iso( tmpX.col(0), n);
  unsigned int b = 1;
  cnt=0;
  for(unsigned int r = 0; r < n; r++){
    if (Z(r,b-1)==0 ) { // need to increase counter by the number of upper-triangle elements we skip
      tmpX.submat(cnt,b,cnt+(n-r)-1,b).fill(0); cnt += n-r;
      continue; }
    for(unsigned int c = r; c < n; c++){
      tmpX(cnt,b) = exp(parameters[2+b] - tmpX(cnt,b) ) * Z(r,b-1) * Z(c,b-1); // replace each element of tmpX
      //tmpX(cnt,b) = exp(parameters[2+b] - tmpX(cnt,b) - pow(Z(r,b-1) - Z(c,b-1),2) ); // replace each element of tmpX
      cnt++;
    }
  }
  // being done with the elements in tmpX we use it to construct the kernel matrices in an efficient way
  Ks.slice(b) = uppertri2symmat_iso( tmpX.col(b), n );
  for(unsigned int j = 0; j < n*(n+1)/2; j++){ tmpX(j,0) += tmpX(j,1); }

  Kfull = uppertri2symmat_iso(tmpX.col(0), n); //sum sparse vectorization instead of matrices

  return Rcpp::List::create(_("full") =  Kfull,
                            _("elements") = Ks);
}


/*
double factorial(double x, double result = 1) {
if (x == 1) return result; else return factorial(x - 1, x * result);
}

double binomial_coef(double n, double k, double result=1){
if(k > 1){
result = (n/k) * binomial_coef(n-1,k-1,result);
} else {
return result;
}
}

double RELU(double x){
if(x>0){
return x;
} else {
return 0;
}
}

double KBn_spline_kernel(double d,unsigned int n){
n = 2 * n + 1;
// See Vapnik et al. (1996) K = int B_n B_n dt = B_{2n+1}
// d = x_i - x_j

double result = 0;
for(unsigned int r = 1; r <= (n+1); r++){
result += (pow(-1,r)/factorial(n)) * binomial_coef(n+1,r) * pow(RELU(d + ((n+1)/2) - r),n);
}
return result;
}
*/
