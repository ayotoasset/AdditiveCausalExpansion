#ifndef ACE_KERNELUTILITIES_H
#define ACE_KERNELUTILITIES_H

#include <RcppArmadillo.h>

// maybe change to template instead of inline:
inline arma::mat uppertri2symmat(const arma::vec& matvec,
                                 const unsigned int& dim){
  //register output matrix
  arma::mat out(dim, dim);

  unsigned int cnt = 0;
  for (unsigned int r = 0; r < dim; r++) {
    for (unsigned int c = r; c < dim; c++) {
      out(r, c) = out(c, r) = matvec(cnt);
      cnt++;
    }
  }
  return out;
}

// general gradient equation, dK is the matrix of elementwise derivatives wrt to a parameter, kernel agnostic
inline double evid_grad(const arma::mat& Kaa,
                        const arma::mat& dK) {
  return - 0.5 * arma::trace(Kaa * dK);
}

// gradient wrt to sigma, kernel agnostic
inline double sigma_gradient(const arma::mat& Kaa, double sigma) {
  return - 0.5 * arma::trace(Kaa) * exp(sigma);
}

inline double logevidence(const arma::vec& y, const arma::colvec& alpha, const arma::vec& eigenval, const unsigned int& n) {
 return - 0.5 * (n * log(2.0 * arma::datum::pi) + sum(log(eigenval)) + arma::dot(y, alpha));

}

inline double sign(double x) {
  return (0 < x) - (x < 0);
}

#endif //ACE_KERNELUTILITIES_H
