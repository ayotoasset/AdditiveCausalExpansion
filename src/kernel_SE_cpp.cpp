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
Rcpp::List invkernel_no_eigen_cpp(arma::mat pdmat, const double sigma){
  unsigned int n = pdmat.n_cols;
  arma::vec eigval(n); eigval.ones();
  arma::mat eigvec(n,n); eigvec.eye();
  pdmat.diag() += exp(sigma);

  //arma::mat invKmat(n,n);
  if(inv_sympd( pdmat, pdmat )){
    Rcout << "Inverse of Kernel not completed." << std::endl;
  }
  return Rcpp::List::create(_("eigenval") = 0,
                            _("inv") = pdmat);
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

  return Rcpp::List::create(_("eigenval") =  eigval,
                            _("inv") = invKmat); //out["eigenvec"] = eigvec;
}

///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// gradients

// [[Rcpp::export]]
double mu_solution_cpp(arma::colvec& y, arma::mat& invKmat) {
  //calculates the exact solutions to the maximization problem
  return  0.5 * arma::sum(invKmat * y ) / arma::sum(arma::sum(invKmat));
}

double evid_grad(const arma::mat& Kaa, const arma::mat& dK) {
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
      L[b + B*i] = - 0.5 * arma::trace( Kaa * (K.slice(b) % tmpX)) * exp( - L[b + B*i] );
    }
  }
  return L;
}

// reduced to only an output list due to specification of the optimizers
// [[Rcpp::export]]
arma::vec grad_SE_cpp(const arma::vec& y, const arma::mat& X, const arma::mat& Z,
                      arma::mat& Kfull, arma::cube& K, arma::mat& invKmatn, arma::vec& eigenval,
                      const arma::vec& parameters,
                      arma::vec& stats, const unsigned int& B) {

  unsigned int n = X.n_rows;
  unsigned int px = X.n_cols;
  //unsigned int B = Z.n_cols+1;
  /*
  if(B != ((parameters.size()-2)/(px+1)) ) {
  throw std::range_error("Inconsistent parameter vector dimensions.");
  }
  */

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

  //clip gradients
  //umat idx=(gradients>20);
  for(unsigned int g=0; g<gradients.size();g++){
    if(gradients(g)>20){
      gradients.row(g) = 20;
    } else if(gradients(g)<-20) {
      gradients.row(g) = -20;
    }
  }


  return gradients;
}

arma::vec evid_hess_lambda(arma::mat mat, const arma::cube& dK, unsigned int b) {
  unsigned int B = dK.n_slices;
  arma::vec hessians(B-b);

  for(unsigned int q=0; q < hessians.size(); q++){
    hessians(q) = - 0.5 * arma::trace( mat * dK.slice(q+b));
  }
  //Rcpp::Rcout << hessians << endl;
  return hessians;
}

arma::mat evid_scale_gradients_hessian(const arma::mat& X,const arma::mat& invK, const arma::mat& Kaa, const arma::cube& K,
                                       arma::vec L, unsigned int B,const arma::mat& Kaa_sigma,const arma::mat& Kaa_lambda, arma::mat& Hessian){
  // produce a (p x B) matrix "grad" with the gradients of L
  unsigned int n = X.n_rows;
  unsigned int p = X.n_cols;
  arma::mat tmpX(n,n);
  arma::mat tmpdK(n,n);
  double tmpL;

  arma::mat tmp(n,n);
  arma::mat tmpX2(n,n);
  //arma::mat tmpdK2(n,n);

  for(unsigned int i=0; i < p; i++){
    for(unsigned int r=0; r < n; r++){
      tmpX.col(r) = pow(X(r, i) - X.col(i),2);
    }
    for(unsigned int b=0; b < B; b++){
      // replace parameter with its gradient
      tmpdK = (K.slice(b) % tmpX);

      tmpL = L[b + B*i];
      L[b + B*i] = - 0.5 * arma::trace( Kaa * tmpdK) * exp( - tmpL );

      tmp = Kaa_lambda * tmpdK * invK;

      //Hessian L-L
      /*
      for(unsigned int j = (p-1); j > i; j--){
      for(unsigned int r2 = 0; r2 < n; r2++){
      tmpX2.col(r2) = pow(X(r2, j) - X.col(j),2);
      }
      for(unsigned int q = b; q < B; q++){
      //Hessian on self (L) for differnt i across all b
      Hessian(2 + q + B*(j+1),2 + b + B*(i+1)) = Hessian(2 + b + B*(i+1),2 + q + B*(j+1))  = - 0.5 * arma::trace( invK *(K.slice(q) % tmpX2) * tmp) * exp( - tmpL ) * exp( - L(q+B*j) );
      }
      }
      */

      //Hessian wrt to sigma (too large)
      //Hessian(0, 2 + b + B*(i+1)) = Hessian(2 + b + B*(i+1),0) = - 0.5 * arma::trace( Kaa_sigma * tmpdK) * exp( - tmpL );

      for(unsigned int q=0; q<B; q++){
      //Hessian wrt to lambda
      Hessian(2+q,2 + b + B*(i+1)) = Hessian(2 + b + B*(i+1),2+q) = - 0.5 * arma::trace( K.slice(q) * tmp) * exp( - tmpL );
      }


      /*
      for(unsigned int q=(b+1); q<B; q++){
        //Hessian on self (L) for differnt b and same i !
        Hessian(2 + q + B*(i+1), 2 + b + B*(i+1)) = Hessian(2 + b + B*(i+1),2 + q + B*(i+1))  = - 0.5 * arma::trace( (K.slice(q) % tmpX) * tmp) * exp( - tmpL ) * exp( - L(q+B*i) );
      }
      */

      //Hessian on self (L)
      Hessian(2 + b + B*(i+1),2 + b + B*(i+1)) = L(b + B*i) - 0.5 * arma::trace( tmpdK * tmp) * exp( - 2 * tmpL );

    }
    }
  return L;
  }

// [[Rcpp::export]]
Rcpp::List grad_SE_Hessian_cpp(const arma::vec& y, const arma::mat& X, const arma::mat& Z,
                               const arma::mat& Kfull, const arma::cube& K, const arma::mat& invKmatn,
                               arma::vec& eigenval,const arma::vec& parameters,
                               arma::vec& stats,const unsigned int& B) {

  unsigned int n = X.n_rows;
  unsigned int px = X.n_cols;
  //unsigned int B = K.n_slices;

  arma::vec gradients(parameters.size());
  arma::mat Hessian(parameters.size(),parameters.size()); Hessian.eye();

  arma::colvec ybar(n);
  arma::colvec alpha(n);

  //Hessian tmp matrices
  arma::mat tmpKH(n,n);
  arma::mat tmpKlambda(n,n);
  arma::mat tmpKL(n,n);

  //Jacobian tmp matrix
  arma::mat tmpK(n,n);

  ybar = y - parameters[1];
  alpha = invKmatn * ybar;
  tmpK = invKmatn - alpha * alpha.t();
  //tmpK = invKmatn - invKmatn * ybar * ybar.t() * invKmat.();

  tmpKlambda = (invKmatn - 2 * alpha * alpha.t());
  tmpKH = invKmatn * tmpKlambda;

  //Jacobian: sigma
  gradients[0] = - 0.5 * arma::sum( tmpK.diag() ) * exp( parameters[0] ); // thus avoid allocating dK
  //Hessian sigma-sigma (too large)
  //Hessian(0,0) = gradients[0] - 0.5 * exp( 2 * parameters[0] ) * arma::sum(tmpKH.diag());

  //Jacobian: lambda
  for(unsigned int b=0;b<B;b++){
    gradients[2+b] = evid_grad(tmpK, K.slice(b));
    //Hessian lambda-lambda
    Hessian.submat(2+b,2+b,2+(B-1),2+b) = evid_hess_lambda(invKmatn * K.slice(b) * tmpKlambda, K , b);
    Hessian.submat(2+b,2+b,2+b,2+(B-1)) = Hessian.submat(2+b,2+b,2+(B-1),2+b).t();
    Hessian(2+b,2+b) += gradients[2+b]; //add gradient to diagonal

    //Hessian: sigma - lambda (too large)
    //Hessian(0,2+b) = Hessian(2+b,0) = exp( parameters[0] ) * evid_grad(tmpKH, K.slice(b));
  }

  //Rcpp::Rcout << 4 << std::endl;

  //Jacobian: Kernel ARD lengths and Hessian: sigma-L and lambda - L
  gradients.rows(2+B,1+B*(px+1)) = evid_scale_gradients_hessian(X, invKmatn ,tmpK, K, parameters.rows(2+B,1+B*(px+1)), B,tmpKH, tmpKlambda, Hessian);

  //Rcpp::Rcout << 5 << std::endl;
  //mu - gradient approach
  gradients[1]=0; //arma::sum( invKmatn * ybar )

  //RMSE
  stats(0) = pow(arma::norm(ybar - (Kfull * alpha)),2);
  //Evidence
  //stats(1) = datum::nan;
  stats(1) = - 0.5 * (n * log( 2.0 * arma::datum::pi ) + sum(log(eigenval)) + arma::dot( ybar, alpha ) ) ;

  //stats is returned as a reference

  //Rcpp::Rcout << Hessian << endl;

  return Rcpp::List::create(_("gradients") =  gradients,
                            _("Hessian") = Hessian);

}




