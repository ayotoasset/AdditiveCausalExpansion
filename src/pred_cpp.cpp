// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List pred_cpp(const arma::vec& y_X, const double sigma, const double mu,
                    const arma::mat& invK_XX, arma::mat& K_xX, arma::mat K_xx,
                    double mean_y, double std_y){
  unsigned int nx = K_xx.n_rows;
  unsigned int nX = invK_XX.n_rows;

  arma::vec y_x(nx);
  arma::vec var(nx);
  arma::mat ci(nx,2);
  arma::mat tmp(nx,nX);

  tmp = K_xX * invK_XX;
  y_x = mean_y + std_y * (tmp * (y_X - mu) + mu);

  K_xx = K_xx - tmp * K_xX.t();
  K_xx.diag() += exp(sigma);
  //absolute value for nuemric stability

  // This is due to numeric instability!
  // TODO: Change calculation of inverse to more numeric stable, e.g. cholesky + solving system
  arma::vec diag_var = K_xx.diag();
  if (arma::any(diag_var < 0)) {
    Rcpp::Rcout << "New lower bound on variance" << std::endl;
    diag_var -= arma::min(diag_var.elem(find(diag_var < 0)));
  }

  var = std_y * sqrt(K_xx.diag());
  ci.col(0) = y_x - 1.96 * var;
  ci.col(1) = y_x + 1.96 * var;
  var = pow(var, 2);

  return Rcpp::List::create(_("map") = y_x,
                            _("ci")  = ci,
                            _("var") = var);
}

// [[Rcpp::export]]
Rcpp::List pred_marginal_cpp(const arma::vec& y_X, const arma::colvec& Z_x, const double sigma, const double mu,
                             const arma::mat& invK_XX,
                             const arma::cube& K_xX, const arma::cube& K_xx,
                             const double& mean_y, const double& std_y, const double& std_Z,
                             bool calculate_ate) {
  // cube argins for x-kernel matrices
  unsigned int nx = K_xx.slice(0).n_rows;
  unsigned int nX = invK_XX.n_rows;
  unsigned int B = K_xx.n_slices;

  arma::vec y_x(nx);
  arma::vec var(nx);
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



  arma::vec diag_var = Kmarg_xx.diag();
  if (arma::any(diag_var < 0)) {
    Rcpp::Rcout << "New lower bound on variance" << std::endl;
    diag_var -= arma::min(diag_var.elem(find(diag_var < 0)));
  }

  //taking absolute value as this can be numerically unstable:
  var = std_y * sqrt(diag_var) / std_Z;
  ci.col(0) = y_x - 1.96 * var;
  ci.col(1) = y_x + 1.96 * var;
  var = pow(var, 2);

  if(!calculate_ate){
    return Rcpp::List::create(_("map") = y_x,
                              _("ci")  = ci,
                              _("var") = var);
  } else {
    // average treatment effect
    double ate = arma::mean(y_x);
    arma::vec ate_ci(2);
    //eficient ci calculation -- assign temporary standard deviation:
    double ate_var = std_y * sqrt(arma::accu(Kmarg_xx)) / nx;
    ate_ci(0) = ate - 1.96 * ate_var;
    ate_ci(1) = ate + 1.96 * ate_var;
    ate_var = pow(ate_var, 2);

    // average treatment on the treated
    unsigned int ntx = arma::sum(Z_x);
    double att = arma::dot(y_x, Z_x) / ntx;
    double att_var = std_y * sqrt(arma::dot(Kmarg_xx * Z_x, Z_x)) / ntx;
    arma::vec att_ci(2);
    att_ci(0) = att - 1.96 * att_var;
    att_ci(1) = att + 1.96 * att_var;
    att_var = pow(att_var, 2);

    // average treatment on the untreated
    unsigned int nux = nx - ntx;
    double atu = (ate * nx - att * ntx) / nux;
    double atu_var = std_y * sqrt(arma::dot(Kmarg_xx * (Z_x == 0), Z_x == 0)) / nux;
    arma::vec atu_ci(2);
    atu_ci(0) = atu - 1.96 * atu_var;
    atu_ci(1) = atu + 1.96 * atu_var;
    atu_var = pow(atu_var, 2);

    return Rcpp::List::create(_("map") = y_x,
                              _("ci")  = ci,
                              _("var") = var,
                              _("ate") = Rcpp::List::create(_("map") = ate,
                                                            _("ci")  = ate_ci,
                                                            _("var") = ate_var),
                              _("att") = Rcpp::List::create(_("map") = att,
                                                            _("ci")  = att_ci,
                                                            _("var") = att_var),
                              _("atu") = Rcpp::List::create(_("map") = atu,
                                                            _("ci")  = atu_ci,
                                                            _("var") = atu_var)
                                );
  }
}
