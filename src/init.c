#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _roll_roll_all(void *, void *, void *, void *, void *, void *);
extern SEXP _roll_roll_any(void *, void *, void *, void *, void *, void *);
extern SEXP _roll_roll_cov(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _roll_roll_crossprod(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _roll_roll_idxquantile(void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _roll_roll_lm(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _roll_roll_mean(void *, void *, void *, void *, void *, void *, void *);
extern SEXP _roll_roll_prod(void *, void *, void *, void *, void *, void *, void *);
extern SEXP _roll_roll_quantile(void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _roll_roll_scale(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _roll_roll_sd(void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _roll_roll_sum(void *, void *, void *, void *, void *, void *, void *);
extern SEXP _roll_roll_var(void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
  {"_roll_roll_all",         (DL_FUNC) &_roll_roll_all,          6},
  {"_roll_roll_any",         (DL_FUNC) &_roll_roll_any,          6},
  {"_roll_roll_cov",         (DL_FUNC) &_roll_roll_cov,         10},
  {"_roll_roll_crossprod",   (DL_FUNC) &_roll_roll_crossprod,   10},
  {"_roll_roll_idxquantile", (DL_FUNC) &_roll_roll_idxquantile,  8},
  {"_roll_roll_lm",          (DL_FUNC) &_roll_roll_lm,           9},
  {"_roll_roll_mean",        (DL_FUNC) &_roll_roll_mean,         7},
  {"_roll_roll_prod",        (DL_FUNC) &_roll_roll_prod,         7},
  {"_roll_roll_quantile",    (DL_FUNC) &_roll_roll_quantile,     8},
  {"_roll_roll_scale",       (DL_FUNC) &_roll_roll_scale,        9},
  {"_roll_roll_sd",          (DL_FUNC) &_roll_roll_sd,           8},
  {"_roll_roll_sum",         (DL_FUNC) &_roll_roll_sum,          7},
  {"_roll_roll_var",         (DL_FUNC) &_roll_roll_var,          8},
  {NULL, NULL, 0}
};

void R_init_roll(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}