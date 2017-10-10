// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
using namespace Rcpp;

arma::mat d(arma::colvec x, arma::vec knots) {
  //returns a vector of d_1 ... d_k

  unsigned int n = x.size();
  unsigned int K = knots.size(); //unique elements
  arma::mat d(n,K);

  //temporarily store constant terms in K
  d.col(K) = (x > knots[K]) * arma::pow(x - knots[K],3);
  for(unsigned int i=1; i < (K - 1); i++){
    d.col(i) = (x > knots[i]) * arma::pow(x - knots[i],3) - d.col(K) ;
    d.col(i) = d.col(i) / (knots[K] - knots[i]);
  }
  //overwrite K-th term with actual value
  d.col(K) = d.col(K-1);

  return d;
}

arma::mat dd(arma::colvec x, arma::vec knots) {
  //returns a vector of the derivatives of d_1 ... d_k

  unsigned int n = x.size();
  knots = unique(knots);
  if(!knots.is_sorted("ascend")){
    knots = sort(knots,"ascend");
  }
  unsigned int K = knots.size(); //unique elements
  arma::mat dd(n,K);

  //temporarily store constant terms in K
  dd.col(K) = 3 * (x > knots[K]) * arma::pow(x - knots[K],2);
  for(unsigned int i=0; i < (K - 1); i++){
    dd.col(i) = 3 * (x > knots[i]) * arma::pow(x - knots[i],2) - dd.col(K) ;
    dd.col(i) = dd.col(i) / (knots[K] - knots[i]);
  }
  //overwrite K-th term with actual value
  dd.col(K) = dd.col(K-1);

  return dd;
}

// [[Rcpp::export]]
arma::mat ncs_basis(arma::colvec x,arma::vec knots) {
  // returns natural cubic spline basis based on the truncated power basis (see p. 145 in Hastie et al., 2001)

  knots = unique(knots);
  if(!knots.is_sorted("ascend")){
    knots = sort(knots,"ascend");
  }

  // get sample size
  unsigned int n = x.size();
  // get basis dimension
  unsigned int B = knots.size() + 1; //(no intercept)
  //allocate memory for the design matrix
  arma::mat design(n,B);
  design.col(0) = x;
  design.submat( 0, 1, n, B ) = d(x,knots);

  return design;
}

// [[Rcpp::export]]
arma::mat ncs_basis_deriv(arma::colvec x,arma::vec knots) {

  knots = unique(knots);
  if(!knots.is_sorted("ascend")){
    knots = sort(knots,"ascend");
  }

  // get sample size
  unsigned int n = x.size();
  // get basis dimension
  unsigned int B = knots.size() + 1; //(no intercept)
  //allocate memory for the design matrix
  arma::mat design(n,B);

  design.col(0).ones(); //was x
  design.submat( 0, 1, n, B ) = dd(x,knots);

  return design;
}
