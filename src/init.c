#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Forward declarations for functions implemented in r_interface.c
SEXP integrate_expression_c(SEXP expression, SEXP a_str, SEXP b_str, SEXP precision,
                            SEXP abs_tol, SEXP rel_tol, SEXP max_eval,
                            SEXP max_depth, SEXP verbose);
SEXP list_builtin_functions_c(void);
SEXP test_flint_basic(void);

// Registration table for C functions callable from R
static const R_CallMethodDef CallEntries[] = {
  {"integrate_expression_c",     (DL_FUNC) &integrate_expression_c,     9},
  {"list_builtin_functions_c",   (DL_FUNC) &list_builtin_functions_c,   0},
  {"test_flint_basic",           (DL_FUNC) &test_flint_basic,           0},
  {NULL, NULL, 0}
};

void R_init_integralARB(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
