#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Function declarations
SEXP integrate_rigorous_c(SEXP f, SEXP a, SEXP b, SEXP precision,
                          SEXP abs_tol, SEXP rel_tol, SEXP max_eval,
                          SEXP max_depth, SEXP verbose);

SEXP test_flint_basic(void);

// Method registration
static const R_CallMethodDef CallEntries[] = {
  {"integrate_rigorous_c", (DL_FUNC) &integrate_rigorous_c, 9},
  {"test_flint_basic", (DL_FUNC) &test_flint_basic, 0},
  {NULL, NULL, 0}
};

void R_init_integralARB(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
