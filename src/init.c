#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _roll_roll_cov(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _roll_roll_eigen(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _roll_roll_lm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _roll_roll_mean(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _roll_roll_pcr(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _roll_roll_prod(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _roll_roll_scale(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _roll_roll_sd(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _roll_roll_sum(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _roll_roll_var(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _roll_roll_vif(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_roll_roll_cov",   (DL_FUNC) &_roll_roll_cov,    9},
  {"_roll_roll_eigen", (DL_FUNC) &_roll_roll_eigen,  9},
  {"_roll_roll_lm",    (DL_FUNC) &_roll_roll_lm,    13},
  {"_roll_roll_mean",  (DL_FUNC) &_roll_roll_mean,   7},
  {"_roll_roll_pcr",   (DL_FUNC) &_roll_roll_pcr,   14},
  {"_roll_roll_prod",  (DL_FUNC) &_roll_roll_prod,   7},
  {"_roll_roll_scale", (DL_FUNC) &_roll_roll_scale,  9},
  {"_roll_roll_sd",    (DL_FUNC) &_roll_roll_sd,     8},
  {"_roll_roll_sum",   (DL_FUNC) &_roll_roll_sum,    7},
  {"_roll_roll_var",   (DL_FUNC) &_roll_roll_var,    8},
  {"_roll_roll_vif",   (DL_FUNC) &_roll_roll_vif,    9},
  {NULL, NULL, 0}
};

void R_init_roll(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}