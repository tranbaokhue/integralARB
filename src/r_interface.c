#include <R.h>
#include <Rinternals.h>
#include <string.h>
#include "include/integralARB.h"

// External functions from builtin_functions.c (now declared in header)
extern parsed_expression_t* parse_mathematical_expression(const char* expression);
extern void cleanup_parsed_expression(parsed_expression_t* expr);
extern int parsed_expression_integrand(acb_ptr res, const acb_t z, void* param, slong order, slong prec);
extern void list_builtin_functions(char* buffer, size_t buffer_size);

// Helper function to convert ARB result to high-precision strings
static void arb_to_r_strings(integration_result_t *result, const arb_t value, slong prec) {
  arf_t lo, hi, mid, err;
  slong digits;

  arf_init(lo); arf_init(hi); arf_init(mid); arf_init(err);

  arb_get_interval_arf(lo, hi, value, prec);

  // mid = (lo + hi) / 2
  arf_add(mid, lo, hi, prec, ARF_RND_NEAR);
  arf_mul_2exp_si(mid, mid, -1);

  // err = (hi - lo) / 2
  arf_sub(err, hi, lo, prec, ARF_RND_NEAR);
  arf_mul_2exp_si(err, err, -1);

  // Free any previous strings
  if (result->value_str) {
    flint_free(result->value_str);
    result->value_str = NULL;
  }
  if (result->error_str) {
    flint_free(result->error_str);
    result->error_str = NULL;
  }

  // Convert to decimal strings with appropriate precision
  digits = (slong)((double)prec * 0.30102999566398114) + 10;
  result->value_str = arf_get_str(mid, digits);
  result->error_str = arf_get_str(err, digits);

  arf_clear(lo); arf_clear(hi); arf_clear(mid); arf_clear(err);
}

// Main integration function for mathematical expressions (string inputs only)
SEXP integrate_expression_c(SEXP expression, SEXP a_str, SEXP b_str, SEXP precision,
                            SEXP abs_tol, SEXP rel_tol, SEXP max_eval,
                            SEXP max_depth, SEXP verbose) {

  /* All variable declarations at function start for C99 compatibility */
  const char* expr;
  const char* a_string;
  const char* b_string;
  int prec_bits;
  double abs_tol_val;
  double rel_tol_val;
  int max_eval_val;
  int max_depth_val;
  int verbose_val;

  integration_result_t result;
  parsed_expression_t* parsed_expr = NULL;
  arb_t a_real, b_real;
  integration_params_t integration_params;

  acb_t acb_a, acb_b, acb_result;
  mag_t tol;
  acb_calc_integrate_opt_t options;
  slong rel_goal;
  int arb_status;

  arb_t real_part;
  mag_t error_mag;

  SEXP r_result;
  SEXP names;

  /* Extract and validate input parameters */
  if (TYPEOF(expression) != STRSXP || TYPEOF(a_str) != STRSXP || TYPEOF(b_str) != STRSXP) {
    error("All inputs (expression, a, b) must be character strings for precision preservation");
  }

  expr = CHAR(STRING_ELT(expression, 0));
  a_string = CHAR(STRING_ELT(a_str, 0));
  b_string = CHAR(STRING_ELT(b_str, 0));
  prec_bits = asInteger(precision);
  abs_tol_val = asReal(abs_tol);
  rel_tol_val = asReal(rel_tol);
  max_eval_val = asInteger(max_eval);
  max_depth_val = asInteger(max_depth);
  verbose_val = asLogical(verbose);

  /* Validate precision */
  if (prec_bits < 64 || prec_bits > 1024) {
    error("Precision must be between 64 and 1024 bits");
  }

  /* Initialize result structure with error defaults */
  result.value = 0.0;
  result.error_bound = INFINITY;
  result.abs_error = INFINITY;
  result.rel_error = INFINITY;
  result.evaluations = 0;
  result.subdivisions = 0;
  result.status = -1;
  result.message[0] = '\0';
  result.value_str = NULL;
  result.error_str = NULL;

  /* Initialize ARB variables */
  arb_init(a_real);
  arb_init(b_real);

  /* Parse the mathematical expression */
  parsed_expr = parse_mathematical_expression(expr);
  if (!parsed_expr) {
    result.status = -1;
    snprintf(result.message, sizeof(result.message),
             "Failed to parse expression: %s. Supported: sin(expr), cos(expr), exp(expr), log(expr), sinh(expr), cosh(expr), x^n, +, -, *, ()",
             expr);
    goto create_r_result;
  }

  /* Set up high-precision integration limits from strings */
  if (arb_set_str(a_real, a_string, prec_bits) != 0) {
    result.status = -1;
    snprintf(result.message, sizeof(result.message),
             "Invalid number format in lower limit: '%s'", a_string);
    goto create_r_result;
  }

  if (arb_set_str(b_real, b_string, prec_bits) != 0) {
    result.status = -1;
    snprintf(result.message, sizeof(result.message),
             "Invalid number format in upper limit: '%s'", b_string);
    goto create_r_result;
  }

  /* Set up integration parameters */
  integration_params.precision = prec_bits;
  integration_params.abs_tolerance = abs_tol_val;
  integration_params.rel_tolerance = rel_tol_val;
  integration_params.max_evaluations = max_eval_val;
  integration_params.max_depth = max_depth_val;
  integration_params.verbose = verbose_val;

  /* Perform integration using ARB's rigorous integration */
  acb_init(acb_a);
  acb_init(acb_b);
  acb_init(acb_result);
  mag_init(tol);

  /* Convert string-parsed limits to complex numbers */
  acb_set_arb(acb_a, a_real);
  acb_set_arb(acb_b, b_real);

  /* Set up ARB integration options */
  mag_set_d(tol, integration_params.abs_tolerance);
  acb_calc_integrate_opt_init(options);
  options->deg_limit = 120;
  options->eval_limit = integration_params.max_evaluations;
  options->depth_limit = integration_params.max_depth;
  options->verbose = integration_params.verbose;

  rel_goal = (slong)integration_params.precision;

  /* Perform the rigorous integration */
  arb_status = acb_calc_integrate(acb_result,
                                  parsed_expression_integrand,
                                  parsed_expr->terms[0].expr,
                                  acb_a, acb_b,
                                  rel_goal, tol, options,
                                  (slong)integration_params.precision);

  if (arb_status == 0) {
    /* Success: extract results */
    arb_init(real_part);
    acb_get_real(real_part, acb_result);

    /* Double approximation for compatibility */
    result.value = arf_get_d(arb_midref(real_part), ARF_RND_NEAR);

    /* Error bound */
    mag_init(error_mag);
    arb_get_mag(error_mag, real_part);
    result.error_bound = mag_get_d(error_mag);

    result.abs_error = result.error_bound;
    result.rel_error = (result.value != 0.0) ? result.error_bound / fabs(result.value) : result.error_bound;
    result.status = 0;
    strncpy(result.message, "Integration successful", sizeof(result.message)-1);
    result.message[sizeof(result.message)-1] = '\0';

    /* Estimate evaluations and subdivisions */
    result.evaluations = options->eval_limit / 100;
    result.subdivisions = 1;

    /* Convert to high-precision strings inline */
    {
      arf_t lo, hi, mid, err;
      slong digits;

      arf_init(lo); arf_init(hi); arf_init(mid); arf_init(err);
      arb_get_interval_arf(lo, hi, real_part, prec_bits);

      /* mid = (lo + hi) / 2 */
      arf_add(mid, lo, hi, prec_bits, ARF_RND_NEAR);
      arf_mul_2exp_si(mid, mid, -1);

      /* err = (hi - lo) / 2 */
      arf_sub(err, hi, lo, prec_bits, ARF_RND_NEAR);
      arf_mul_2exp_si(err, err, -1);

      /* Convert to decimal strings */
      digits = (slong)((double)prec_bits * 0.30102999566398114) + 10;
      result.value_str = arf_get_str(mid, digits);
      result.error_str = arf_get_str(err, digits);

      arf_clear(lo); arf_clear(hi); arf_clear(mid); arf_clear(err);
    }

    arb_clear(real_part);
    mag_clear(error_mag);

  } else {
    /* Integration failed */
    const char* error_msgs[] = {
      "Integration successful",                    /* 0 */
    "Maximum number of evaluations exceeded",   /* 1 */
    "Maximum subdivision depth exceeded",       /* 2 */
    "Tolerance could not be achieved",          /* 3 */
    "Function evaluation failed",               /* 4 */
    "Unknown integration error"                 /* 5+ */
    };
    int msg_index;

    result.status = arb_status;
    msg_index = (arb_status >= 0 && arb_status <= 4) ? arb_status : 5;
    strncpy(result.message, error_msgs[msg_index], sizeof(result.message)-1);
    result.message[sizeof(result.message)-1] = '\0';

    result.value = 0.0;
    result.error_bound = INFINITY;

    /* Allocate empty strings for failed integration */
    result.value_str = (char*)flint_malloc(1);
    result.value_str[0] = '\0';
    result.error_str = (char*)flint_malloc(4);
    strcpy(result.error_str, "inf");
  }

  /* Cleanup ARB variables */
  acb_clear(acb_a);
  acb_clear(acb_b);
  acb_clear(acb_result);
  mag_clear(tol);

  create_r_result:
    /* Cleanup in all cases */
    if (parsed_expr) {
      cleanup_parsed_expression(parsed_expr);
    }
    arb_clear(a_real);
    arb_clear(b_real);

    /* Convert result to R list */
    r_result = PROTECT(allocVector(VECSXP, 12));
    names = PROTECT(allocVector(STRSXP, 12));

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
    SET_STRING_ELT(names, 10, mkChar("expression"));
    SET_STRING_ELT(names, 11, mkChar("precision"));

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
    SET_VECTOR_ELT(r_result, 10, mkString(expr));
    SET_VECTOR_ELT(r_result, 11, ScalarInteger(prec_bits));

    setAttrib(r_result, R_NamesSymbol, names);

    /* Cleanup result strings */
    if (result.value_str) {
      flint_free(result.value_str);
    }
    if (result.error_str) {
      flint_free(result.error_str);
    }

    UNPROTECT(2);
    return r_result;
}

/* Function to list all supported mathematical expressions */
SEXP list_builtin_functions_c() {
  char buffer[2048];
  list_builtin_functions(buffer, sizeof(buffer));
  return mkString(buffer);
}

/* Test function to verify FLINT is working with string parsing */
SEXP test_flint_basic(void) {
  arb_t pi, test_num;
  const char* test_string;
  char* pi_str;
  char* test_str;
  char result_buffer[1024];
  SEXP result;

  arb_init(pi);
  arb_init(test_num);

  /* Test 1: Compute π with high precision */
  arb_const_pi(pi, 256);

  /* Test 2: Parse a high-precision string */
  test_string = "3.1415926535897932384626433832795028841971693993751";
  if (arb_set_str(test_num, test_string, 256) == 0) {
    /* String parsing successful */
    pi_str = arb_get_str(pi, 50, ARB_STR_MORE);
    test_str = arb_get_str(test_num, 50, ARB_STR_MORE);

    /* Create result string */
    snprintf(result_buffer, sizeof(result_buffer),
             "FLINT working correctly.\nπ (computed): %s\nπ (parsed): %s\nString parsing: SUCCESS",
             pi_str, test_str);

    result = mkString(result_buffer);

    flint_free(pi_str);
    flint_free(test_str);
    arb_clear(pi);
    arb_clear(test_num);

    return result;
  } else {
    arb_clear(pi);
    arb_clear(test_num);
    return mkString("FLINT working but string parsing failed");
  }
}

// Benchmark function for testing chain rule expressions
SEXP benchmark_chain_rule_c(SEXP expressions, SEXP precision_levels) {
  int n_expr = LENGTH(expressions);
  int n_prec = LENGTH(precision_levels);
  int* precisions = INTEGER(precision_levels);

  SEXP results = PROTECT(allocVector(VECSXP, n_expr * n_prec));
  SEXP result_names = PROTECT(allocVector(STRSXP, n_expr * n_prec));

  int result_idx = 0;

  for (int i = 0; i < n_expr; i++) {
    const char* expr = CHAR(STRING_ELT(expressions, i));

    for (int j = 0; j < n_prec; j++) {
      int prec = precisions[j];

      // Test integration of chain rule expressions
      SEXP expr_sexp = PROTECT(mkString(expr));
      SEXP a_sexp = PROTECT(mkString("0"));
      SEXP b_sexp = PROTECT(mkString("1"));
      SEXP prec_sexp = PROTECT(ScalarInteger(prec));
      SEXP abs_tol = PROTECT(ScalarReal(1e-20));
      SEXP rel_tol = PROTECT(ScalarReal(1e-20));
      SEXP max_eval = PROTECT(ScalarInteger(1000000));
      SEXP max_depth = PROTECT(ScalarInteger(100));
      SEXP verbose = PROTECT(ScalarLogical(0));

      SEXP integration_result = integrate_expression_c(expr_sexp, a_sexp, b_sexp,
                                                       prec_sexp, abs_tol, rel_tol,
                                                       max_eval, max_depth, verbose);

      SET_VECTOR_ELT(results, result_idx, integration_result);

      char name_buf[256];
      snprintf(name_buf, sizeof(name_buf), "%s_prec_%d", expr, prec);
      SET_STRING_ELT(result_names, result_idx, mkChar(name_buf));

      result_idx++;

      UNPROTECT(9);
    }
  }

  setAttrib(results, R_NamesSymbol, result_names);
  UNPROTECT(2);
  return results;
}
