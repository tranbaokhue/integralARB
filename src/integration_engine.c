#include "include/integralARB.h"

#include <flint/flint.h>
#include <flint/arb.h>
#include <flint/arf.h>
#include <flint/acb.h>
#include <flint/acb_calc.h>

#include <R.h>
#include <Rinternals.h>

#include <math.h>
#include <string.h>
#include <stdlib.h>

/* Helper: convert arb_t to midpoint/error decimal strings and store in result.
 Strings are allocated by arf_get_str and must be freed with flint_free. */
static void integration_result_set_strings_from_arb(integration_result_t *result,
                                                    const arb_t x,
                                                    slong prec_bits,
                                                    slong digits)
{
  arf_t lo, hi, mid, err;
  arf_init(lo); arf_init(hi);
  arf_init(mid); arf_init(err);

  /* Get interval endpoints to given precision */
  arb_get_interval_arf(lo, hi, x, prec_bits);

  /* mid = (lo + hi) / 2  -- use explicit rounding mode */
  {
    int _ret = arf_add(mid, lo, hi, prec_bits, ARF_RND_NEAR);
    (void)_ret; /* ignore inexactness flag if you don't need it */
  arf_mul_2exp_si(mid, mid, -1);
  }

  /* err = (hi - lo) / 2  -- use explicit rounding mode */
  {
    int _ret = arf_sub(err, hi, lo, prec_bits, ARF_RND_NEAR);
    (void)_ret;
    arf_mul_2exp_si(err, err, -1);
  }

  /* Free any previous strings stored in result (defensive) */
  if (result->value_str) { flint_free(result->value_str); result->value_str = NULL; }
  if (result->error_str) { flint_free(result->error_str); result->error_str = NULL; }

  /* Convert ARF -> decimal strings (arf_get_str uses FLINT allocator) */
  result->value_str = arf_get_str(mid, digits);
  result->error_str = arf_get_str(err, digits);

  arf_clear(lo); arf_clear(hi); arf_clear(mid); arf_clear(err);
}

/* Free FLINT-allocated result strings (safe to call even if NULL). */
void cleanup_integration_result(integration_result_t* result) {
  if (!result) return;
  if (result->value_str) {
    flint_free(result->value_str);
    result->value_str = NULL;
  }
  if (result->error_str) {
    flint_free(result->error_str);
    result->error_str = NULL;
  }
}

/* Minimal integrand cleanup placeholder */
void cleanup_integrand(integrand_t* integrand) {
  (void)integrand;
}

/* Data passed to the R integrand wrapper */
typedef struct {
  SEXP r_function;
  SEXP r_env;
  integration_params_t* params;
} r_function_data_t;

/* R integrand wrapper:
 - calls the R function with the midpoint of the REAL part of z as a double.
 - expects a numeric scalar back (converted to double).
 - sets res accordingly (real arb with imag = 0).
 */
int r_function_integrand(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
  r_function_data_t *rdata = (r_function_data_t*) param;
  if (!rdata || rdata->r_function == R_NilValue) {
    acb_indeterminate(res);
    return 1;
  }

  if (order > 1) {             /* derivatives not supported here */
acb_indeterminate(res);
    return 1;
  }

  /* Convert midpoint of real(z) -> double (user agreed double inputs for now) */
  double x = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);

  /* Call R: f(x) */
  SEXP call = PROTECT(lang2(rdata->r_function, ScalarReal(x)));
  int r_error = 0;
  SEXP rval = R_tryEval(call, rdata->r_env, &r_error);
  UNPROTECT(1);

  if (r_error || rval == R_NilValue) {
    acb_indeterminate(res);
    return 1;
  }

  /* Coerce to real vector and extract first element (protect the coerced SEXP) */
  SEXP rnum = PROTECT(coerceVector(rval, REALSXP));
  double y = REAL(rnum)[0];
  UNPROTECT(1);

  /* Put the result into an arb and set res = that real arb */
  arb_t tmp;
  arb_init(tmp);
  arb_set_d(tmp, y);
  acb_set_arb(res, tmp);
  arb_clear(tmp);

  return 0;
}

/* Simple builtin integrands (examples). These are suitable argument examples
 for testing or default builtin use. */
int builtin_arctangent(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
  (void)param;
  if (order > 1) { acb_indeterminate(res); return 1; }
  acb_atan(res, z, prec);
  return 0;
}

int builtin_gaussian(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
  (void)param;
  if (order > 1) { acb_indeterminate(res); return 1; }
  acb_t t;
  acb_init(t);
  acb_mul(t, z, z, prec);
  acb_neg(t, t);
  acb_exp(res, t, prec);
  acb_clear(t);
  return 0;
}

/* Main integration function. */
integration_result_t perform_integration(integrand_t* integrand,
                                         double a, double b,
                                         integration_params_t* params)
{
  integration_result_t result;

  /* safe defaults */
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

  if (!params) {
    result.status = -1;
    strncpy(result.message, "Null params passed", sizeof(result.message)-1);
    result.message[sizeof(result.message)-1] = '\0';
    return result;
  }

  /* Arb variables */
  acb_t acb_a, acb_b, acb_result;
  mag_t tol;
  acb_calc_integrate_opt_t options;

  acb_init(acb_a);
  acb_init(acb_b);
  acb_init(acb_result);
  mag_init(tol);

  /* endpoints */
  acb_set_d(acb_a, a);
  acb_set_d(acb_b, b);

  /* absolute tolerance (mag_t) */
  mag_set_d(tol, params->abs_tolerance);

  /* options */
  acb_calc_integrate_opt_init(options);
  options->deg_limit = 60;
  options->eval_limit = params->max_evaluations;
  options->depth_limit = params->max_depth;
  options->verbose = params->verbose;

  /* choose integrand */
  acb_calc_func_t integrand_func = NULL;
  void *integrand_param = NULL;
  r_function_data_t rdata;
  if (integrand && integrand->type == INTEGRAND_R_FUNCTION) {
    integrand_func = r_function_integrand;
    rdata.r_function = integrand->r_function;
    rdata.r_env = R_GlobalEnv;
    rdata.params = params;
    integrand_param = &rdata;
  } else {
    integrand_func = builtin_arctangent; /* default builtin */
  integrand_param = NULL;
  }

  /* relative goal: pass precision as a heuristic (common usage) */
  slong rel_goal = (slong) params->precision;

  int status = acb_calc_integrate(acb_result,
                                  integrand_func,
                                  integrand_param,
                                  acb_a,
                                  acb_b,
                                  rel_goal,
                                  tol,
                                  options,
                                  (slong) params->precision);

  if (status == 0) {
    /* Success: extract real part (arb) and populate result */
    arb_t real_part;
    arb_init(real_part);
    acb_get_real(real_part, acb_result);

    result.value = arf_get_d(arb_midref(real_part), ARF_RND_NEAR);

    mag_t error_mag;
    mag_init(error_mag);
    arb_get_mag(error_mag, real_part);
    result.error_bound = mag_get_d(error_mag);

    result.abs_error = result.error_bound;
    result.rel_error = (result.value != 0.0) ? result.error_bound / fabs(result.value)
      : result.error_bound;
    result.status = 0;
    strncpy(result.message, "Integration successful", sizeof(result.message)-1);
    result.message[sizeof(result.message)-1] = '\0';

    /* diagnostics (basic) */
    result.evaluations = options->eval_limit;
    result.subdivisions = 1;

    /* decimal digits ≈ bits * log10(2), plus margin */
    slong digits = (slong)((double)params->precision * 0.30102999566398114) + 5;

    integration_result_set_strings_from_arb(&result, real_part, (slong)params->precision, digits);

    arb_clear(real_part);
    mag_clear(error_mag);
  } else {
    /* Failure: set defaults and give FLINT-allocated strings so cleanup is safe */
    result.status = status;
    strncpy(result.message, "Integration failed - check function or increase precision",
            sizeof(result.message)-1);
    result.message[sizeof(result.message)-1] = '\0';
    result.value = 0.0;
    result.error_bound = INFINITY;

    result.value_str = (char *) flint_malloc(1);
    result.value_str[0] = '\0';
    result.error_str = (char *) flint_malloc(4);
    strcpy(result.error_str, "inf");
  }

  /* locals cleanup */
  acb_clear(acb_a);
  acb_clear(acb_b);
  acb_clear(acb_result);
  mag_clear(tol);

  return result;
}
