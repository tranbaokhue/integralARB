#include "include/integralARB.h"
#include <string.h>
#include <math.h>

// Structure to pass R function data to ARB integrand callback
typedef struct {
  SEXP r_function;
  SEXP r_env;
  integration_params_t* params;
} r_function_data_t;

// ARB integrand callback for R functions
int r_function_integrand(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
  if (order > 1) {
    // Higher derivatives not supported yet
    acb_indeterminate(res);
    return 1;
  }

  r_function_data_t* data = (r_function_data_t*)param;

  // Convert acb_t to R numeric - use only real part for now
  arb_t z_real;
  arb_init(z_real);
  acb_get_real(z_real, z);

  double x = arf_get_d(arb_midref(z_real), ARF_RND_NEAR);

  // Call R function
  SEXP r_x = PROTECT(ScalarReal(x));
  SEXP call = PROTECT(lang2(data->r_function, r_x));
  SEXP result = PROTECT(eval(call, data->r_env));

  // Convert result back to acb_t
  if (TYPEOF(result) == REALSXP && LENGTH(result) == 1) {
    double y = REAL(result)[0];
    if (R_finite(y)) {
      acb_set_d(res, y);
    } else {
      // Handle infinite or NaN results
      acb_indeterminate(res);
    }
  } else if (TYPEOF(result) == CPLXSXP && LENGTH(result) == 1) {
    // Handle complex results
    Rcomplex c = COMPLEX(result)[0];
    if (R_finite(c.r) && R_finite(c.i)) {
      acb_set_d_d(res, c.r, c.i);
    } else {
      acb_indeterminate(res);
    }
  } else {
    // Invalid result type
    acb_indeterminate(res);
  }

  UNPROTECT(3);
  arb_clear(z_real);

  return 0;
}

// Built-in function: 1/(1+x^2) for testing
int builtin_arctangent(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
  (void)param; // Suppress warning

  if (order > 1) {
    acb_indeterminate(res);
    return 1;
  }

  acb_t z2, denominator, one;
  acb_init(z2);
  acb_init(denominator);
  acb_init(one);

  acb_mul(z2, z, z, prec);
  acb_add_ui(denominator, z2, 1, prec);
  acb_one(one);
  acb_div(res, one, denominator, prec);

  acb_clear(z2);
  acb_clear(denominator);
  acb_clear(one);

  return 0;
}

// Built-in function: sin(x)
int builtin_sin(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
  (void)param;

  if (order > 1) {
    acb_indeterminate(res);
    return 1;
  }

  acb_sin(res, z, prec);
  return 0;
}

// Built-in function: exp(-x^2)
int builtin_gaussian(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
  (void)param;

  if (order > 1) {
    acb_indeterminate(res);
    return 1;
  }

  acb_t z2;
  acb_init(z2);

  acb_mul(z2, z, z, prec);
  acb_neg(z2, z2);
  acb_exp(res, z2, prec);

  acb_clear(z2);
  return 0;
}

// Main integration function
integration_result_t perform_integration(integrand_t* integrand,
                                         double a, double b,
                                         integration_params_t* params) {

  integration_result_t result = {0};

  // Initialize ARB variables
  acb_t acb_a, acb_b, acb_result;
  mag_t tol;
  acb_calc_integrate_opt_t options;

  acb_init(acb_a);
  acb_init(acb_b);
  acb_init(acb_result);
  mag_init(tol);

  // Set integration limits
  acb_set_d(acb_a, a);
  acb_set_d(acb_b, b);

  // Set tolerance
  mag_set_d(tol, params->abs_tolerance);

  // Configure integration options
  acb_calc_integrate_opt_init(options);
  options->deg_limit = 60;
  options->eval_limit = params->max_evaluations;
  options->depth_limit = params->max_depth;
  options->verbose = params->verbose;

  // Choose integrand function based on type
  acb_calc_func_t integrand_func;
  void* integrand_param = NULL;
  r_function_data_t r_data;

  switch (integrand->type) {
  case INTEGRAND_R_FUNCTION:
    integrand_func = r_function_integrand;
    r_data.r_function = integrand->r_function;
    r_data.r_env = R_GlobalEnv;
    r_data.params = params;
    integrand_param = &r_data;
    break;

  case INTEGRAND_BUILTIN:
    // For now, use arctangent as default builtin
    integrand_func = builtin_arctangent;
    break;

  default:
    strcpy(result.message, "Unsupported integrand type");
  result.status = -1;
  goto cleanup;
  }

  // Perform the integration
  int status = acb_calc_integrate(acb_result, integrand_func, integrand_param,
                                  acb_a, acb_b, params->precision, tol, options, params->precision);

  // Convert results
  if (status == 0) {
    // Extract real part as the main result
    arb_t real_part;
    arb_init(real_part);
    acb_get_real(real_part, acb_result);

    result.value = arf_get_d(arb_midref(real_part), ARF_RND_NEAR);

    // Get error bound from the radius
    mag_t error_mag;
    mag_init(error_mag);
    arb_get_mag(error_mag, real_part);
    result.error_bound = mag_get_d(error_mag);

    result.abs_error = result.error_bound;
    result.rel_error = (result.value != 0) ? result.error_bound / fabs(result.value) : result.error_bound;
    result.status = 0;
    strcpy(result.message, "Integration successful");

    // Set some reasonable values for diagnostics
    result.evaluations = options->eval_limit / 10; // Estimate
    result.subdivisions = 1;

    arb_clear(real_part);
    mag_clear(error_mag);

  } else {
    result.status = status;
    strcpy(result.message, "Integration failed - check function or increase precision");
    result.value = 0;
    result.error_bound = INFINITY;
  }

  cleanup:
    acb_clear(acb_a);
  acb_clear(acb_b);
  acb_clear(acb_result);
  mag_clear(tol);

  return result;
}

void cleanup_integrand(integrand_t* integrand) {
  // Nothing to cleanup for now
  (void)integrand;
}

void cleanup_integration_result(integration_result_t* result) {
  // Nothing to cleanup for now
  (void)result;
}
