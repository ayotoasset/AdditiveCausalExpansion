#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _ace_Adam_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ace_grad_Matern_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ace_grad_SE_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ace_invkernel_cpp(SEXP, SEXP);
extern SEXP _ace_kernmat_Matern32_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ace_kernmat_Matern32_symmetric_cpp(SEXP, SEXP, SEXP);
extern SEXP _ace_kernmat_SE_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ace_kernmat_SE_symmetric_cpp(SEXP, SEXP, SEXP);
extern SEXP _ace_mu_solution_cpp(SEXP, SEXP);
extern SEXP _ace_Nadam_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ace_ncs_basis(SEXP, SEXP);
extern SEXP _ace_ncs_basis_deriv(SEXP, SEXP);
extern SEXP _ace_Nesterov_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ace_norm_clip_cpp(SEXP, SEXP, SEXP);
extern SEXP _ace_normalize_test(SEXP, SEXP, SEXP);
extern SEXP _ace_normalize_train(SEXP, SEXP, SEXP);
extern SEXP _ace_pred_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ace_pred_marginal_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ace_stats_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_ace_Adam_cpp",                       (DL_FUNC) &_ace_Adam_cpp,                        9},
  {"_ace_grad_Matern_cpp",                (DL_FUNC) &_ace_grad_Matern_cpp,                11},
  {"_ace_grad_SE_cpp",                    (DL_FUNC) &_ace_grad_SE_cpp,                    11},
  {"_ace_invkernel_cpp",                  (DL_FUNC) &_ace_invkernel_cpp,                   2},
  {"_ace_kernmat_Matern32_cpp",           (DL_FUNC) &_ace_kernmat_Matern32_cpp,            5},
  {"_ace_kernmat_Matern32_symmetric_cpp", (DL_FUNC) &_ace_kernmat_Matern32_symmetric_cpp,  3},
  {"_ace_kernmat_SE_cpp",                 (DL_FUNC) &_ace_kernmat_SE_cpp,                  5},
  {"_ace_kernmat_SE_symmetric_cpp",       (DL_FUNC) &_ace_kernmat_SE_symmetric_cpp,        3},
  {"_ace_mu_solution_cpp",                (DL_FUNC) &_ace_mu_solution_cpp,                 2},
  {"_ace_Nadam_cpp",                      (DL_FUNC) &_ace_Nadam_cpp,                       9},
  {"_ace_ncs_basis",                      (DL_FUNC) &_ace_ncs_basis,                       2},
  {"_ace_ncs_basis_deriv",                (DL_FUNC) &_ace_ncs_basis_deriv,                 2},
  {"_ace_Nesterov_cpp",                   (DL_FUNC) &_ace_Nesterov_cpp,                    5},
  {"_ace_norm_clip_cpp",                  (DL_FUNC) &_ace_norm_clip_cpp,                   3},
  {"_ace_normalize_test",                 (DL_FUNC) &_ace_normalize_test,                  3},
  {"_ace_normalize_train",                (DL_FUNC) &_ace_normalize_train,                 3},
  {"_ace_pred_cpp",                       (DL_FUNC) &_ace_pred_cpp,                        8},
  {"_ace_pred_marginal_cpp",              (DL_FUNC) &_ace_pred_marginal_cpp,              11},
  {"_ace_stats_cpp",                      (DL_FUNC) &_ace_stats_cpp,                       6},
  {NULL, NULL, 0}
};

void R_init_ace(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
