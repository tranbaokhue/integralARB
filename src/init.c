#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Function declarations
SEXP integrate_rigorous_c(SEXP f, SEXP a, SEXP b, SEXP precision,
                          SEXP abs_tol, SEXP rel_tol, SEXP max_eval,
                          SEXP max_depth, SEXP verbose);

SEXP integrate_rigorous_builtin_c(SEXP func_name, SEXP a, SEXP b, SEXP precision,
                                  SEXP abs_tol, SEXP rel_tol, SEXP max_eval,
                                  SEXP max_depth, SEXP verbose, SEXP params);

SEXP test_flint_basic(void);
SEXP list_builtin_functions_c(void);
SEXP benchmark_integration_c(SEXP func_name, SEXP precision_levels);

// Method registration
static const R_CallMethodDef CallEntries[] = {
  {"integrate_rigorous_c", (DL_FUNC) &integrate_rigorous_c, 9},
  {"integrate_rigorous_builtin_c", (DL_FUNC) &integrate_rigorous_builtin_c, 10},
  {"test_flint_basic", (DL_FUNC) &test_flint_basic, 0},
  {"list_builtin_functions_c", (DL_FUNC) &list_builtin_functions_c, 0},
  {"benchmark_integration_c", (DL_FUNC) &benchmark_integration_c, 2},
  {NULL, NULL, 0}
};

void R_init_integralARB(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

// Updated src/Makevars.in to include builtin_functions.c
// PKG_CPPFLAGS = @PKG_CFLAGS@ -I./include -DHAVE_FLINT
// PKG_LIBS = @PKG_LIBS@
//
// # C source files to compile - now includes builtin_functions.c
// OBJECTS = init.o r_interface.o integration_engine.o builtin_functions.o
//
// # Use C11 standard to avoid warnings about C11 extensions
// PKG_CFLAGS = -std=c11 -Wno-unused-parameter
