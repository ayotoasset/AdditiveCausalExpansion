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
bool Newton_cpp(double iter, double learn_rate, double momentum, arma::vec& nu, arma::vec grad, const arma::mat& Hessian, arma::vec& para, const unsigned int& B){

  bool output_flag = arma::is_finite(grad); // if all gradients finite


  if(output_flag==false){
    Rcout << "Element of the gradient is not finite: " << find_nonfinite(grad) << std::endl;
    grad.elem( find_nonfinite(grad) ).zeros();}

  grad.rows(2,1+B) = solve(Hessian.submat(2,2,1+B,1+B),grad.rows(2,1+B));
  for(unsigned int j= 2+B; j < grad.size(); j++){
    grad(j) = grad(j) / Hessian(j,j);
  }
  //grad = solve(Hessian,grad);


  nu = momentum * nu + learn_rate * grad;
  para = para + nu;

  return output_flag;
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
