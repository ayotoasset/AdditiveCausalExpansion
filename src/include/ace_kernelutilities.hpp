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

#endif //ACE_KERNELUTILITIES_H
