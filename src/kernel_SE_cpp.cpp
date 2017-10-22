// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <cmath>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List kernmat_SE_cpp(const arma::mat& X1,const arma::mat& X2,const arma::mat& Z1,const arma::mat& Z2,
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
  tmpX.slice(0) = exp(parameters[2] - tmpX.slice(0)); //lambda[0]
  Kfull = tmpX.slice(0);
  for(unsigned int b = 1; b < B; b++){
    for(unsigned int r = 0; r < n1; r++){
      if (Z1(r,b-1)==0 ) {tmpX.slice(b).row(r).fill(0); continue;}
      for(unsigned int c = 0; c < n2; c++){
        if (Z2(c,b-1)==0 ) {tmpX(r,c,b)=0; continue;}
        tmpX(r,c,b) = exp(parameters[2+b] - tmpX(r,c,b)) * Z1(r,b-1) * Z2(c,b-1); //lambda[b]
        //tmpX(r,c,b) = exp(parameters[2+b] - tmpX(r,c,b) - pow(Z1(r,b-1) - Z2(c,b-1),2) ); // replace each element of tmpX
      }
    }
    Kfull += tmpX.slice(b);
  }
  return Rcpp::List::create(_("full") =  Kfull,
                            _("elements") = tmpX);
}

//more efficient with pointers and a class but well
arma::mat uppertri2symmat(arma::vec matvec,unsigned int dim){
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
Rcpp::List kernmat_SE_symmetric_cpp(const arma::mat& X, const arma::mat& Z, const arma::vec& parameters) {
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

  tmpX.col(0) = exp(parameters[2] - tmpX.col(0) ); Ks.slice(0) = uppertri2symmat( tmpX.col(0), n);
  for(unsigned int b = 1; b < B; b++){
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
    Ks.slice(b) = uppertri2symmat( tmpX.col(b), n );
    for(unsigned int j = 0; j < n*(n+1)/2; j++){ tmpX(j,0) += tmpX(j,b); }
  }
  Kfull = uppertri2symmat(tmpX.col(0), n); //sum sparse vectorization instead of matrices

  return Rcpp::List::create(_("full") =  Kfull,
                            _("elements") = Ks);
}


// [[Rcpp::export]]
Rcpp::List invkernel_cpp(arma::mat pdmat, const double sigma){
  unsigned int n = pdmat.n_cols;
  arma::vec eigval(n); eigval.ones();
  arma::mat eigvec(n,n); eigvec.eye();
  pdmat.diag() += exp(sigma);

  if(!arma::eig_sym( eigval, eigvec, pdmat)){
    Rcout << "Eigenvalue decomp. not completed." << std::endl;
  }

  //get inverse and eigenvalues
  arma::mat invKmat(n,n);
  invKmat = eigvec;
  for(unsigned int i = 0; i < n; i++){
    invKmat.col(i) = invKmat.col(i) / eigval[i];
  }
  invKmat = invKmat * eigvec.t();

  if( eigval[0]<=0 ){
    Rcout << "Smallest eigenvalue is negative or zero! - This should not happen as sigma puts a lower bound on the eigenvalues." << std::endl;
    Rcout << "Sigma^2: " << exp(sigma) << std::endl;
  }

  return Rcpp::List::create(_("eigenval") =  eigval,
                            _("inv") = invKmat); //out["eigenvec"] = eigvec;
}

///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// gradients

// [[Rcpp::export]]
double mu_solution_cpp(arma::colvec& y, arma::mat& invKmat) {
  //calculates the exact solutions to the maximization problem
  return  0.5 * arma::sum(invKmat * y ) / arma::sum(arma::sum(invKmat));
}

double evid_grad(arma::mat& Kaa, arma::mat& dK) {
  return - 0.5 * arma::trace( Kaa * dK);
}

arma::mat evid_scale_gradients(const arma::mat& X, const arma::mat& Kaa, const arma::cube& K, arma::vec L, unsigned int B){
  // produce a (p x B) matrix "grad" with the gradients of L
  unsigned int n = X.n_rows;
  unsigned int p = X.n_cols;
  arma::mat tmpX(n,n);

  for(unsigned int i=0; i < p; i++){
    for(unsigned int r=0; r < n; r++){
      tmpX.col(r) = pow(X(r, i) - X.col(i),2);
    }
    for(unsigned int b=0; b < B; b++){
      // replace parameter with its gradient
      L[b + B*i] = - 0.5 * arma::trace( Kaa * (K.slice(b) % tmpX) * exp( - L[b + B*i] ) );
    }
  }
  return L;
}

// reduced to only an output list due to specification of the optimizers
// [[Rcpp::export]]
arma::vec grad_SE_cpp(arma::vec& y, arma::mat& X, arma::mat& Z,
                          arma::mat& Kfull, arma::cube& K, arma::mat& invKmatn, arma::vec& eigenval,
                          const arma::vec& parameters,
                          arma::vec& stats) {

  unsigned int n = X.n_rows;
  unsigned int px = X.n_cols;
  unsigned int B = Z.n_cols+1;

  if(B != ((parameters.size()-2)/(px+1)) ) {
    throw std::range_error("Inconsistent parameter vector dimensions.");
  }

  arma::vec gradients(parameters.size());

  arma::colvec ybar(n);
  arma::colvec alpha(n);
  arma::mat tmpK(n,n);

  ybar = y - parameters[1];
  alpha = invKmatn * ybar;
  tmpK = invKmatn - alpha * alpha.t();

  //sigma
  gradients[0] = - 0.5 * arma::sum( tmpK.diag() ) * exp( parameters[0] ); // thus avoid allocating dK

  //lambda
  for(unsigned int b=0;b<B;b++){
    gradients[2+b] = evid_grad(tmpK, K.slice(b));
  }

  //L
  gradients.rows(2+B,1+B*(px+1)) = evid_scale_gradients(X, tmpK, K, parameters.rows(2+B,1+B*(px+1)), B);

  //mu - gradient approach
  gradients[1]=0; //arma::sum( invKmatn * ybar )

  //RMSE
  stats(0) = pow(arma::norm(ybar - (Kfull * alpha)),2);
  //Evidence
  stats(1) = - 0.5 * (n * log( 2.0 * arma::datum::pi ) + sum(log(eigenval)) + arma::dot( ybar, alpha ) ) ;
  //stats is returned as a reference

  return gradients;
}

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

///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// credible intervals
// [[Rcpp::export]]
Rcpp::List pred_cpp(const arma::vec& y_X, const double sigma, const double mu,
                    const arma::mat& invK_XX, arma::mat& K_xX,arma::mat K_xx,
                    double mean_y, double std_y){
  unsigned int nx = K_xx.n_rows;
  unsigned int nX = invK_XX.n_rows;

  arma::vec y_x(nx);
  arma::mat ci(nx,2);
  arma::mat tmp(nx,nX);

  tmp = K_xX * invK_XX;
  y_x = mean_y + std_y *tmp * (y_X - mu) + mu;
  //y_x = mean_y + std_y * y_x;

  K_xx = K_xx - tmp * K_xX.t();
  K_xx.diag() += exp(sigma);

  ci.col(1) = std_y*1.96 * sqrt(K_xx.diag());
  ci.col(0) = y_x - ci.col(1);
  ci.col(1) = y_x + ci.col(1);

  return Rcpp::List::create(_("map") = y_x,
                            _("ci") = ci);
}

// [[Rcpp::export]]
Rcpp::List pred_marginal_cpp(const arma::vec& y_X, const double sigma, const double mu,
                             const arma::mat& invK_XX,
                             const arma::cube& K_xX, const arma::cube& K_xx,
                             const double& mean_y, const double& std_y,
                             bool calculate_ate){ // cube argins for x-kernel matrices
  unsigned int nx = K_xx.slice(0).n_rows;
  unsigned int nX = invK_XX.n_rows;
  unsigned int B = K_xx.n_slices;

  arma::vec y_x(nx);
  arma::mat ci(nx,2);
  arma::mat tmp(nx,nX);

  arma::mat Kmarg_xX(nx,nX); Kmarg_xX = K_xX.slice(1);
  arma::mat Kmarg_xx(nx,nx); Kmarg_xx = K_xx.slice(1);
  for(unsigned int b = 2; b < B; b++){ // not the "constant" nuisance term and b=1
    Kmarg_xX += K_xX.slice(b);
    Kmarg_xx += K_xx.slice(b);
  }
  //invK_XX.diag() = invK_XX.diag() + 0.00001;
  tmp = Kmarg_xX * invK_XX;
  y_x = std_y * tmp * (y_X - mu);
  //y_x = std_y * y_x;

  Kmarg_xx = Kmarg_xx - tmp * Kmarg_xX.t(); //saving storage

  //taking absolute value as this can be numerically unstable:
  ci.col(1) = std_y * 1.96 * sqrt(abs(Kmarg_xx.diag()));
  ci.col(0) = y_x - ci.col(1);
  ci.col(1) = y_x + ci.col(1);

  if(!calculate_ate){
    return Rcpp::List::create(_("map") = y_x,
                              _("ci") = ci);
  } else {
    double ate = arma::mean(y_x);
    arma::vec ate_ci(2);

    //eficient ci calculation
    ate_ci(1) = std_y * sqrt(arma::sum(arma::sum(Kmarg_xx))/pow(nx,2));
    ate_ci(0) = ate - 1.96*ate_ci(1);
    ate_ci(1) = ate + 1.96*ate_ci(1);

    return Rcpp::List::create(_("map") = y_x,
                              _("ci") = ci,
                              _("ate_map") = ate,
                              _("ate_ci") = ate_ci);
  }
}
