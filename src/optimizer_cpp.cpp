// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
bool Nesterov_cpp(double learn_rate,
                  double momentum,
                  arma::vec& nu,
                  const arma::vec& grad,
                  arma::vec& para) {

  bool output_flag = arma::is_finite(grad); // if all gradients finite

  nu = momentum * nu + learn_rate * grad;
  para = para + nu;

  return output_flag;
}

// [[Rcpp::export]]
bool Nadam_cpp(double iter,
               double learn_rate,
               double beta1,
               double beta2,
               double eps,
               arma::vec& m,
               arma::vec& v,
               const arma::vec& grad,
               arma::vec& para) {

  bool output_flag = arma::is_finite(grad);

  //use references to update
  m = beta1 * m + (1 - beta1) * grad;
  v = beta2 * v + (1 - beta2) * arma::pow(grad, 2);
  para = para + learn_rate * ((beta1 * m  + (1 - beta1) * grad ) / (1 - pow(beta1, iter)) ) /
                                              (sqrt(v / (1 - pow(beta2, iter))) + eps);

  return output_flag;
}

// [[Rcpp::export]]
bool Adam_cpp(double iter,
              double learn_rate,
              double beta1,
              double beta2,
              double eps,
              arma::vec& m,
              arma::vec& v,
              const arma::vec& grad,
              arma::vec& para) {

  bool output_flag = arma::is_finite(grad);

  //use references to update
  m = (beta1 * m) + (1-beta1) * grad;
  v = beta2 * v + (1-beta2) * arma::pow(grad,2);
  para = para + learn_rate * (m  / (1-pow(beta1,iter)) ) / ( sqrt(v / (1-pow(beta2,iter)) ) + eps);

  return output_flag;
}
