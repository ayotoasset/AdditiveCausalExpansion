// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

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
  y_x = mean_y + std_y * (tmp * (y_X - mu) + mu);

  K_xx = K_xx - tmp * K_xX.t();
  K_xx.diag() += exp(sigma);
  //absolute value for nuemric stability
  ci.col(1) = std_y * 1.96 * sqrt(abs(K_xx.diag()));
  ci.col(0) = y_x - ci.col(1);
  ci.col(1) = y_x + ci.col(1);

  return Rcpp::List::create(_("map") = y_x,
                            _("ci") = ci);
}

// [[Rcpp::export]]
Rcpp::List pred_marginal_cpp(const arma::vec& y_X, const double sigma, const double mu,
                             const arma::mat& invK_XX,
                             const arma::cube& K_xX, const arma::cube& K_xx,
                             const double& mean_y, const double& std_y, const double& std_Z,
                             bool calculate_ate){ // cube argins for x-kernel matrices
  unsigned int nx = K_xx.slice(0).n_rows;
  unsigned int nX = invK_XX.n_rows;
  unsigned int B = K_xx.n_slices;

  arma::vec y_x(nx);
  arma::mat ci(nx,2);
  arma::mat tmp(nx,nX);

  arma::mat Kmarg_xX(nx,nX);
  arma::mat Kmarg_xx(nx,nx);

  if(B>1){ //non-additive kernel (regular GP)
    Kmarg_xX = K_xX.slice(1);
    Kmarg_xx = K_xx.slice(1);
    if(B>2){
      for(unsigned int b = 2; b < B; b++){ // not the "constant" nuisance term and b=1
        Kmarg_xX += K_xX.slice(b);
        Kmarg_xx += K_xx.slice(b);
      }
    }
  } else if (B==1){
    Kmarg_xX = K_xX.slice(0);
    Kmarg_xx = K_xx.slice(0);
  }

  tmp = Kmarg_xX * invK_XX;
  y_x = std_y * tmp * (y_X - mu) / std_Z;

  Kmarg_xx = Kmarg_xx - tmp * Kmarg_xX.t(); //saving storage

  //taking absolute value as this can be numerically unstable:
  ci.col(1) =  std_y * 1.96 * sqrt(abs(Kmarg_xx.diag())) / std_Z;
  ci.col(0) = y_x - ci.col(1);
  ci.col(1) = y_x + ci.col(1);

  if(!calculate_ate){
    return Rcpp::List::create(_("map") = y_x,
                              _("ci") = ci);
  } else {
    double ate = arma::mean(y_x);
    arma::vec ate_ci(2);

    //eficient ci calculation
    ate_ci(1) = std_y * sqrt(arma::accu(Kmarg_xx)) / nx;
    ate_ci(0) = ate - 1.96*ate_ci(1);
    ate_ci(1) = ate + 1.96*ate_ci(1);

    return Rcpp::List::create(_("map") = y_x,
                              _("ci") = ci,
                              _("ate_map") = ate,
                              _("ate_ci") = ate_ci);
  }
}
