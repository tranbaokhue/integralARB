#include <R.h>
#include <Rinternals.h>
#include <string.h>
#include "include/integralARB.h"

// External function from builtin_functions.c
extern const builtin_function_entry_t* find_builtin_function(const char* name);
extern void list_builtin_functions(char* buffer, size_t buffer_size);

// Built-in function integration interface
SEXP integrate_rigorous_builtin_c(SEXP func_name, SEXP a, SEXP b, SEXP precision,
                                  SEXP abs_tol, SEXP rel_tol, SEXP max_eval,
                                  SEXP max_depth, SEXP verbose, SEXP params) {

  // Input validation and conversion
  const char* fname = CHAR(STRING_ELT(func_name, 0));
  double a_val = asReal(a);
  double b_val = asReal(b);
  int prec_bits = asInteger(precision);
  double abs_tol_val = asReal(abs_tol);
  double rel_tol_val = asReal(rel_tol);
  int max_eval_val = asInteger(max_eval);
  int max_depth_val = asInteger(max_depth);
  int verbose_val = asLogical(verbose);

  // Find the built-in function
  const builtin_function_entry_t* func_entry = find_builtin_function(fname);
  if (!func_entry) {
    error("Unknown built-in function: %s", fname);
  }

  // Set up integration parameters
  integration_params_t integration_params = {
    .precision = prec_bits,
    .abs_tolerance = abs_tol_val,
    .rel_tolerance = rel_tol_val,
    .max_evaluations = max_eval_val,
    .max_depth = max_depth_val,
    .verbose = verbose_val
  };

  // Set up integrand for built-in function
  integrand_t integrand = {
    .type = INTEGRAND_BUILTIN,
    .data = NULL,  // Most functions don't need parameters
    .r_function = R_NilValue,
    .builtin_func = func_entry->func
  };

  // Handle parameterized functions (like gamma_n)
  int param_value = 0;
  if (strcmp(fname, "gamma_0") == 0) {
    param_value = 0;
    integrand.data = &param_value;
  } else if (strcmp(fname, "gamma_1") == 0) {
    param_value = 1;
    integrand.data = &param_value;
  } else if (strcmp(fname, "gamma_2") == 0) {
    param_value = 2;
    integrand.data = &param_value;
  }

  // Perform integration
  integration_result_t result = perform_integration(&integrand, a_val, b_val, &integration_params);

  // Convert result to R list
  SEXP r_result = PROTECT(allocVector(VECSXP, 12));
  SEXP names = PROTECT(allocVector(STRSXP, 12));

  SET_STRING_ELT(names, 0, mkChar("value"));
  SET_STRING_ELT(names, 1, mkChar("error_bound"));
  SET_STRING_ELT(names, 2, mkChar("abs_error"));
  SET_STRING_ELT(names, 3, mkChar("rel_error"));
  SET_STRING_ELT(names, 4, mkChar("evaluations"));
  SET_STRING_ELT(names, 5, mkChar("subdivisions"));
  SET_STRING_ELT(names, 6, mkChar("status"));
  SET_STRING_ELT(names, 7, mkChar("message"));
  SET_STRING_ELT(names, 8, mkChar("value_str"));
  SET_STRING_ELT(names, 9, mkChar("error_str"));
  SET_STRING_ELT(names, 10, mkChar("function_name"));
  SET_STRING_ELT(names, 11, mkChar("function_description"));

  SET_VECTOR_ELT(r_result, 0, ScalarReal(result.value));
  SET_VECTOR_ELT(r_result, 1, ScalarReal(result.error_bound));
  SET_VECTOR_ELT(r_result, 2, ScalarReal(result.abs_error));
  SET_VECTOR_ELT(r_result, 3, ScalarReal(result.rel_error));
  SET_VECTOR_ELT(r_result, 4, ScalarInteger(result.evaluations));
  SET_VECTOR_ELT(r_result, 5, ScalarInteger(result.subdivisions));
  SET_VECTOR_ELT(r_result, 6, ScalarInteger(result.status));
  SET_VECTOR_ELT(r_result, 7, mkString(result.message));
  SET_VECTOR_ELT(r_result, 8, mkString(result.value_str ? result.value_str : ""));
  SET_VECTOR_ELT(r_result, 9, mkString(result.error_str ? result.error_str : ""));
  SET_VECTOR_ELT(r_result, 10, mkString(func_entry->name));
  SET_VECTOR_ELT(r_result, 11, mkString(func_entry->description));

  setAttrib(r_result, R_NamesSymbol, names);

  // Cleanup
  cleanup_integrand(&integrand);
  cleanup_integration_result(&result);

  UNPROTECT(2);
  return r_result;
}

// Legacy R function interface (now deprecated)
SEXP integrate_rigorous_c(SEXP f, SEXP a, SEXP b, SEXP precision,
                          SEXP abs_tol, SEXP rel_tol, SEXP max_eval,
                          SEXP max_depth, SEXP verbose) {

  // For backward compatibility - redirect to built-in functions if possible
  if (TYPEOF(f) == STRSXP) {
    const char* fname = CHAR(STRING_ELT(f, 0));
    if (find_builtin_function(fname)) {
      return integrate_rigorous_builtin_c(f, a, b, precision, abs_tol, rel_tol,
                                          max_eval, max_depth, verbose, R_NilValue);
    }
  }

  // If not a built-in function, return error with helpful message
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

  SET_VECTOR_ELT(r_result, 0, ScalarReal(0.0));
  SET_VECTOR_ELT(r_result, 1, ScalarReal(R_PosInf));
  SET_VECTOR_ELT(r_result, 2, ScalarReal(R_PosInf));
  SET_VECTOR_ELT(r_result, 3, ScalarReal(R_PosInf));
  SET_VECTOR_ELT(r_result, 4, ScalarInteger(0));
  SET_VECTOR_ELT(r_result, 5, ScalarInteger(0));
  SET_VECTOR_ELT(r_result, 6, ScalarInteger(-1));
  SET_VECTOR_ELT(r_result, 7, mkString("Only built-in analytical functions supported. Use integrate_rigorous() with function name."));

  setAttrib(r_result, R_NamesSymbol, names);

  UNPROTECT(2);
  return r_result;
}

// Function to list all available built-in functions
SEXP list_builtin_functions_c() {
  char buffer[4096];
  list_builtin_functions(buffer, sizeof(buffer));
  return mkString(buffer);
}

// Test function to verify FLINT is working
SEXP test_flint_basic(void) {
  arb_t pi;
  arb_init(pi);

  // Compute π with high precision
  arb_const_pi(pi, 128);

  // Convert to string for high precision
  char* pi_str = arb_get_str(pi, 50, ARB_STR_MORE);
  SEXP result = mkString(pi_str);

  flint_free(pi_str);
  arb_clear(pi);

  return result;
}

// Benchmark function to test integration accuracy
SEXP benchmark_integration_c(SEXP func_name, SEXP precision_levels) {
  const char* fname = CHAR(STRING_ELT(func_name, 0));
  int* precisions = INTEGER(precision_levels);
  int n_prec = LENGTH(precision_levels);

  // Known exact result for arctangent(0,1) = π/4
  if (strcmp(fname, "arctangent") != 0) {
    error("Benchmark only supports 'arctangent' function currently");
  }

  SEXP results = PROTECT(allocVector(VECSXP, n_prec));
  SEXP result_names = PROTECT(allocVector(STRSXP, n_prec));

  for (int i = 0; i < n_prec; i++) {
    int prec = precisions[i];

    // Set up parameters for high accuracy
    integration_params_t params = {
      .precision = prec,
      .abs_tolerance = pow(2, -prec - 10),
      .rel_tolerance = pow(2, -prec - 10),
      .max_evaluations = 1000000,
      .max_depth = 100,
      .verbose = 0
    };

    // Get arctangent function
    const builtin_function_entry_t* func_entry = find_builtin_function("arctangent");

    integrand_t integrand = {
      .type = INTEGRAND_BUILTIN,
      .data = NULL,
      .r_function = R_NilValue,
      .builtin_func = func_entry->func
    };

    // Integrate from 0 to 1
    integration_result_t result = perform_integration(&integrand, 0.0, 1.0, &params);

    // Create result list for this precision level
    SEXP prec_result = PROTECT(allocVector(VECSXP, 4));
    SEXP prec_names = PROTECT(allocVector(STRSXP, 4));

    SET_STRING_ELT(prec_names, 0, mkChar("precision"));
    SET_STRING_ELT(prec_names, 1, mkChar("value_str"));
    SET_STRING_ELT(prec_names, 2, mkChar("error_str"));
    SET_STRING_ELT(prec_names, 3, mkChar("status"));

    SET_VECTOR_ELT(prec_result, 0, ScalarInteger(prec));
    SET_VECTOR_ELT(prec_result, 1, mkString(result.value_str ? result.value_str : ""));
    SET_VECTOR_ELT(prec_result, 2, mkString(result.error_str ? result.error_str : ""));
    SET_VECTOR_ELT(prec_result, 3, ScalarInteger(result.status));

    setAttrib(prec_result, R_NamesSymbol, prec_names);

    SET_VECTOR_ELT(results, i, prec_result);

    char name_buf[32];
    snprintf(name_buf, sizeof(name_buf), "precision_%d", prec);
    SET_STRING_ELT(result_names, i, mkChar(name_buf));

    cleanup_integration_result(&result);
    UNPROTECT(2);
  }

  setAttrib(results, R_NamesSymbol, result_names);
  UNPROTECT(2);
  return results;
}
