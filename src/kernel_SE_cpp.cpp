// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <cmath>

using namespace arma;
using namespace Rcpp;

double evid_grad(arma::mat& mymat, arma::mat& dK) {
  return - 0.5 * trace( mymat * dK);
}

// [[Rcpp::export]]
Rcpp::List kernmat_SE_cpp(const arma::mat& X1,const arma::mat& X2,const arma::mat& Z1,const arma::mat& Z2,
                          arma::vec lambda, arma::mat& L) {
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
      tmprow = pow(X1(r,i) - conv_to<rowvec>::from(X2.col(i)), 2);
      for(unsigned int b = 0; b < B; b++){ //3rd dimension of array is called slice
        tmpX.slice(b).row(r) += tmprow * exp(- L(i,b) );
      }
    }
  }

  //outest-most iteration is basis dimension
  tmpX.slice(0) = exp(lambda[0] - tmpX.slice(0));
  Kfull = tmpX.slice(0);
  for(unsigned int b = 1; b < B; b++){
    for(unsigned int r = 0; r < n1; r++){
      if (Z1(r,b-1)==0 ) {tmpX.slice(b).row(r).fill(0); continue;}
      for(unsigned int c = 0; c < n2; c++){
        if (Z2(c,b-1)==0 ) {tmpX(r,c,b)=0; continue;}
        tmpX(r,c,b) = exp(lambda[b] - tmpX(r,c,b)) * Z1(r,b-1) * Z2(c,b-1);
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
      out(r,c) = matvec(cnt);
      out(c,r) = matvec(cnt);
      cnt++;
    }
  }
  return out;
}


// [[Rcpp::export]]
Rcpp::List kernmat_SE_symmetric_cpp(const arma::mat& X, const arma::mat& Z, arma::vec lambda, arma::mat& L) {
  //input X1,X2,Z1,Z2,parameters
  unsigned int n = X.n_rows;
  unsigned int p = X.n_cols;
  unsigned int B = Z.n_cols + 1; //including nuisance term
  unsigned int cnt;

  //parameters
  //arma::mat L = para["L"]; // p x B
  //arma::vec lambda = para["lambda"];

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
          tmpX(cnt,b) += tmp * exp(- L(i,b) );
        } cnt++;
      }
    }
  }

  tmpX.col(0) = exp(lambda[0] - tmpX.col(0) ); Ks.slice(0) = uppertri2symmat( tmpX.col(0), n);
  for(unsigned int b = 1; b < B; b++){
    cnt=0;
    for(unsigned int r = 0; r < n; r++){
      if (Z(r,b-1)==0 ) { // need to increase counter by the number of upper-triangle elements we skip
        tmpX.submat(cnt,b,cnt+(n-r)-1,b).fill(0); cnt += n-r;
        continue; }
      for(unsigned int c = r; c < n; c++){
        tmpX(cnt,b) = exp(lambda[b] - tmpX(cnt,b) ) * Z(r,b-1) * Z(c,b-1); // replace each element of tmpX
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
arma::rowvec stats_SE(arma::colvec y, arma::mat& Kmat, arma::mat& invKmatn, arma::vec& eigenval, double mu) {
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
  double val = sum(log(eigenval));
  //double sign; arma::log_det(val,sign,invKmatn);  //logdet of INVERSE -> -val
  stats(1) = - 0.5 * (n * log( 2.0 * arma::datum::pi ) + val + arma::dot( y, alpha ) ) ;

  return stats;
}

// [[Rcpp::export]]
Rcpp::List invkernel_cpp(arma::mat pdmat, const double& sigma){
  unsigned int n = pdmat.n_cols;
  arma::vec eigval(n);
  arma::mat eigvec(n,n);
  pdmat.diag() += exp(sigma);
  arma::eig_sym( eigval, eigvec, pdmat);

  arma::mat invKmat(n,n); invKmat = eigvec;
  for(unsigned int i = 0; i < n; i++){
    invKmat.col(i) = invKmat.col(i) / eigval[i];
  }
  invKmat = invKmat * eigvec.t();
  return Rcpp::List::create(_("eigenval") =  eigval,
                            _("inv") = invKmat); //out["eigenvec"] = eigvec;
}

// [[Rcpp::export]]
double mu_solution_cpp(arma::colvec& y, arma::mat& invKmat) {
  //calculates the exact solutions to the maximization problem
  //double mu = 0.5 * arma::sum(invKmat * y ) / arma::sum(arma::sum(invKmat));
  //return mu;
  return  0.5 * arma::sum(invKmat * y ) / arma::sum(arma::sum(invKmat));
}

///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///// /////


/*

 Rcpp::List evid_scale_gradients(arma::mat X, arma::mat tmpK,
 Rcpp::List K,
Rcpp::List parameters){
// produce a (p x 2) matrix "grad" with the gradients of Lm and La in the first and second column, respectively
unsigned int n = X.n_rows;
unsigned int p = X.n_cols;

//predefine matrices and temporary storage files
arma::vec grad_m(p);
arma::vec grad_a(p);
arma::mat tmpX(n,n);
arma::mat tmpM(n,n);
arma::mat tmpA(n,n);

arma::mat L = para["L"]; // p x B
arma::vec lambda = para["lambda"];

for(unsigned int i=0; i<p; i++){
//tmpX = abs(matrix(rep(X[,i],n),n,n) - t(matrix(rep(X[,i],n),n,n)))
for(unsigned int r=0; r< n; r++){
tmpX.col(r) = pow(X(r, i) - X.col(i),2);
}

//Klist$Km * exp( 2*log(tmpX) - parameters$Lm[i] )
// NEED TO REWRITE FOR MATRIX L
tmpM = Km % tmpX * exp(- as<arma::vec>(L[b])[i]);
//tmpA = Ka % tmpX * exp(- as<arma::vec>(parameters["La"])[i]);

grad_m(i) = evid_grad(tmpK,tmpM);
//Rcpp::Rcout << grad_m(i) << std::endl;
grad_a(i) = evid_grad(tmpK,tmpA);
}

return Rcpp::List::create(Named("m") = grad_m,Named("a") = grad_a );
}


*/

/*
// [[Rcpp::export]]
Rcpp::List grad_GP_SE_cpp(arma::colvec y, arma::mat X, arma::colvec z,arma::colvec w, arma::mat Kmat, arma::mat Km, arma::mat Ka, arma::mat invKmatn, Rcpp::List parameters) {
  Rcpp::List gradients = clone(parameters); //for the same list structure

  unsigned int n = X.n_rows;
  //unsigned int p = X.n_cols;

  //preallocate memory
  arma::colvec stats(2);
  arma::colvec ybar(n);
  arma::colvec alpha(n);
  arma::mat tmpK(n,n);
  arma::mat dK(n,n);

  ybar = y - as<double>(parameters["mu"]); // - parameters["mu_z"]
  alpha = invKmatn * ybar;
  tmpK = invKmatn - alpha * alpha.t();

  //sigma
  dK = arma::diagmat( exp( as<double>(parameters["sigma"]) ) / w ); //+ as<double>(parameters["sigma_z"]) * z
  gradients["sigma"] = evid_grad(tmpK, dK);

  //sigma_z
  //dK.diag() = dK.diag() % z;
  //gradients["sigma_z"] = evid_grad(tmpK, dK);

  //lambdam
  gradients["lambdam"] = evid_grad(tmpK, Km);

  //lambdaa
  gradients["lambdaa"] = evid_grad(tmpK, Ka);

  //Lm and La
  Rcpp::List tmp = evid_scale_gradients(X, tmpK, Km, Ka, parameters);
  //arma::vec tmp(2,fill::zeros);

  gradients["Lm"] = as<std::vector<double> >(tmp["m"]);
  gradients["La"] = as<std::vector<double> >(tmp["a"]);

  //mu - gradient approach

  gradients["mu"] = arma::sum( invKmatn * ybar );

  //mu_z?
  //gradients["mu"] = z.t() * invKmatn * ybar ;

  //Might remove stats but it is more efficient to calculate it here

  //RMSE
  stats(0) = pow(arma::norm(y - (Kmat * alpha + as<double>(parameters["mu"]))),2);

  //Evidence
  double val;
  double sign;
  arma::log_det(val,sign,invKmatn);  //logdet of INVERSE -> -val
  stats(1) = - 0.5 * (n * log( 2.0 * arma::datum::pi ) - val + arma::dot( ybar, alpha ) ) ;

  //output
  return Rcpp::List::create(Named("gradients") = gradients,Named("stats") = stats );
}






*/
