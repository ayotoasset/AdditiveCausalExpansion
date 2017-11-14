// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;
//#include <Rcpp.h>
//using namespace Rcpp;

// [[Rcpp::export]]
bool Nesterov_cpp(double learn_rate, double momentum, arma::vec& nu, arma::vec grad, arma::vec& para){

  bool output_flag = arma::is_finite(grad); // if all gradients finite

  if(output_flag==false){
    Rcout << "Element that is not finite: " << find_nonfinite(grad) << std::endl;
    grad.elem( find_nonfinite(grad) ).zeros();}

  nu = momentum * nu + learn_rate * grad;
  para = para + nu;

  return output_flag;
  //return para;
}

// [[Rcpp::export]]
bool Newton_cpp(double iter, double learn_rate, double momentum, arma::vec& nu, arma::vec grad, arma::mat Hessian, arma::vec& para){

  bool output_flag = arma::is_finite(grad); // if all gradients finite
  // if one element of the Hessian is not finite, use the identity matrix
  if(arma::is_finite(Hessian)){
    if(~arma::inv_sympd(Hessian,Hessian)){
      Rcout << "Hessian not invertible: using identity matrix" << find_nonfinite(grad) << std::endl;
      Hessian.eye();
    }
  }

  if(output_flag==false){
    Rcout << "Element that is not finite: " << find_nonfinite(grad) << std::endl;
    grad.elem( find_nonfinite(grad) ).zeros();}

  //Hessian.diag() += 0.001;
  nu = momentum * nu + learn_rate * Hessian * grad;
  para = para + nu;

  return output_flag;
  //return para;
}

// [[Rcpp::export]]
bool Nadam_cpp(double iter,double learn_rate,double beta1, double beta2, double eps, arma::vec& m, arma::vec& v, const arma::vec& grad, arma::vec& para){

  bool output_flag = arma::is_finite(grad);

  //use references to update
  m = (beta1 * m) + (1-beta1) * grad;
  v = beta2 * v + (1-beta2) * arma::pow(grad,2);
  para = para + learn_rate * ( ( beta1 * m  + (1-beta1)*grad ) / (1-pow(beta1,iter)) ) / ( sqrt(v / (1-pow(beta2,iter)) ) + eps);

  return output_flag;
}

// [[Rcpp::export]]
bool Adam_cpp(double iter,double learn_rate,double beta1, double beta2, double eps, arma::vec& m, arma::vec& v, const arma::vec& grad, arma::vec& para){

  bool output_flag = arma::is_finite(grad);

  //use references to update
  m = (beta1 * m) + (1-beta1) * grad;
  v = beta2 * v + (1-beta2) * arma::pow(grad,2);
  para = para + learn_rate * ( m  / (1-pow(beta1,iter)) ) / ( sqrt(v / (1-pow(beta2,iter)) ) + eps);

  return output_flag;
}
