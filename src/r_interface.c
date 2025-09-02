#include <R.h>
#include <Rinternals.h>
#include "include/integralARB.h"

// Main integration function called from R
SEXP integrate_rigorous_c(SEXP f, SEXP a, SEXP b, SEXP precision,
                          SEXP abs_tol, SEXP rel_tol, SEXP max_eval,
                          SEXP max_depth, SEXP verbose) {

  // Input validation and conversion
  double a_val = asReal(a);
  double b_val = asReal(b);
  int prec_bits = asInteger(precision);
  double abs_tol_val = asReal(abs_tol);
  double rel_tol_val = asReal(rel_tol);
  int max_eval_val = asInteger(max_eval);
  int max_depth_val = asInteger(max_depth);
  int verbose_val = asLogical(verbose);

  // Set up integration parameters
  integration_params_t params = {
    .precision = prec_bits,
    .abs_tolerance = abs_tol_val,
    .rel_tolerance = rel_tol_val,
    .max_evaluations = max_eval_val,
    .max_depth = max_depth_val,
    .verbose = verbose_val
  };

  // Set up integrand based on input type
  integrand_t integrand = {0};

  if (TYPEOF(f) == CLOSXP || TYPEOF(f) == BUILTINSXP || TYPEOF(f) == SPECIALSXP) {
    // R function
    integrand.type = INTEGRAND_R_FUNCTION;
    integrand.r_function = f;
  } else if (TYPEOF(f) == STRSXP) {
    // String - for now treat as builtin function name
    const char* fname = CHAR(STRING_ELT(f, 0));
    if (strcmp(fname, "arctangent") == 0 || strcmp(fname, "builtin") == 0) {
      integrand.type = INTEGRAND_BUILTIN;
    } else {
      // For other strings, we'd need expression parsing (future enhancement)
      integrand.type = INTEGRAND_R_FUNCTION;
      // Create a function that parses the expression (simplified for now)
      integrand.r_function = f; // This will need proper expression parsing
    }
  } else {
    // Default to builtin for testing
    integrand.type = INTEGRAND_BUILTIN;
  }

  // Perform integration
  integration_result_t result = perform_integration(&integrand, a_val, b_val, &params);

  // Convert result to R list
  SEXP r_result = PROTECT(allocVector(VECSXP, 8));
  SEXP names = PROTECT(allocVector(STRSXP, 8));

  SET_STRING_ELT(names, 0, mkChar("value"));
  SET_STRING_ELT(names, 1, mkChar("error_bound"));
  SET_STRING_ELT(names, 2, mkChar("abs_error"));
  SET_STRING_ELT(names, 3, mkChar("rel_error"));
  SET_STRING_ELT(names, 4, mkChar("evaluations"));
  SET_STRING_ELT(names, 5, mkChar("subdivisions"));
  SET_STRING_ELT(names, 6, mkChar("status"));
  SET_STRING_ELT(names, 7, mkChar("message"));

  SET_VECTOR_ELT(r_result, 0, ScalarReal(result.value));
  SET_VECTOR_ELT(r_result, 1, ScalarReal(result.error_bound));
  SET_VECTOR_ELT(r_result, 2, ScalarReal(result.abs_error));
  SET_VECTOR_ELT(r_result, 3, ScalarReal(result.rel_error));
  SET_VECTOR_ELT(r_result, 4, ScalarInteger(result.evaluations));
  SET_VECTOR_ELT(r_result, 5, ScalarInteger(result.subdivisions));
  SET_VECTOR_ELT(r_result, 6, ScalarInteger(result.status));
  SET_VECTOR_ELT(r_result, 7, mkString(result.message));

  setAttrib(r_result, R_NamesSymbol, names);

  // Cleanup
  cleanup_integrand(&integrand);
  cleanup_integration_result(&result);

  UNPROTECT(2);
  return r_result;
}

// Test function to verify FLINT is working
SEXP test_flint_basic(void) {
  // Simple test: compute π using ARB
  arb_t pi;
  arb_init(pi);

  // Compute π with 64-bit precision
  arb_const_pi(pi, 64);

  // Convert to double for R
  double pi_value = arf_get_d(arb_midref(pi), ARF_RND_NEAR);

  arb_clear(pi);

  return ScalarReal(pi_value);
}
