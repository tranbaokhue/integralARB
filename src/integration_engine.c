#include "include/integralARB.h"
#include <string.h>
#include <math.h>

// Helper function to convert arb result to strings with high precision
static void integration_result_set_strings_from_arb(integration_result_t *result,
                                                    const arb_t x,
                                                    slong prec_bits,
                                                    slong digits)
{
  // Get interval endpoints
  arf_t lo, hi, mid, err;
  arf_init(lo); arf_init(hi);
  arf_init(mid); arf_init(err);

  arb_get_interval_arf(lo, hi, x, prec_bits);

  // mid = (lo + hi) / 2
  arf_add(mid, lo, hi, prec_bits, ARF_RND_NEAR);
  arf_mul_2exp_si(mid, mid, -1);

  // err = (hi - lo) / 2
  arf_sub(err, hi, lo, prec_bits, ARF_RND_NEAR);
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

  // Convert to decimal strings
  result->value_str = arf_get_str(mid, digits);
  result->error_str = arf_get_str(err, digits);

  arf_clear(lo); arf_clear(hi); arf_clear(mid); arf_clear(err);
}

// Cleanup functions
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

void cleanup_integrand(integrand_t* integrand) {
  // Nothing to cleanup for built-in functions
  (void)integrand;
}

// Main integration function
integration_result_t perform_integration(integrand_t* integrand,
                                         double a, double b,
                                         integration_params_t* params) {

  integration_result_t result;

  // Initialize result with safe defaults
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

  if (!integrand) {
    result.status = -1;
    strncpy(result.message, "Null integrand passed", sizeof(result.message)-1);
    result.message[sizeof(result.message)-1] = '\0';
    return result;
  }

  // ARB variables
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

  // Configure integration options for high accuracy
  acb_calc_integrate_opt_init(options);
  options->deg_limit = 120;  // Higher degree limit for better accuracy
  options->eval_limit = params->max_evaluations;
  options->depth_limit = params->max_depth;
  options->verbose = params->verbose;

  // Choose integrand function
  acb_calc_func_t integrand_func = NULL;
  void* integrand_param = NULL;

  if (integrand->type == INTEGRAND_BUILTIN) {
    integrand_func = integrand->builtin_func;
    integrand_param = integrand->data;
  } else {
    result.status = -1;
    strncpy(result.message, "Only built-in functions supported", sizeof(result.message)-1);
    result.message[sizeof(result.message)-1] = '\0';
    goto cleanup;
  }

  if (!integrand_func) {
    result.status = -1;
    strncpy(result.message, "Invalid integrand function", sizeof(result.message)-1);
    result.message[sizeof(result.message)-1] = '\0';
    goto cleanup;
  }

  // Perform the integration with high precision
  slong rel_goal = (slong)params->precision;

  int status = acb_calc_integrate(acb_result,
                                  integrand_func,
                                  integrand_param,
                                  acb_a,
                                  acb_b,
                                  rel_goal,
                                  tol,
                                  options,
                                  (slong)params->precision);

  if (status == 0) {
    // Success: extract real part and populate result
    arb_t real_part;
    arb_init(real_part);
    acb_get_real(real_part, acb_result);

    // Double approximation for compatibility
    result.value = arf_get_d(arb_midref(real_part), ARF_RND_NEAR);

    // Error bound
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

    // Basic diagnostics - these are estimates
    result.evaluations = options->eval_limit / 100;  // Rough estimate
    result.subdivisions = 1;

    // High-precision string representations
    slong digits = (slong)((double)params->precision * 0.30102999566398114) + 10;
    integration_result_set_strings_from_arb(&result, real_part,
                                            (slong)params->precision, digits);

    arb_clear(real_part);
    mag_clear(error_mag);

  } else {
    // Integration failed
    result.status = status;

    const char* error_msgs[] = {
      "Integration successful",                    // 0
      "Maximum number of evaluations exceeded",   // 1
      "Maximum subdivision depth exceeded",       // 2
      "Tolerance could not be achieved",          // 3
      "Function evaluation failed",               // 4
      "Unknown integration error"                 // 5+
    };

    int msg_index = (status >= 0 && status <= 4) ? status : 5;
    strncpy(result.message, error_msgs[msg_index], sizeof(result.message)-1);
    result.message[sizeof(result.message)-1] = '\0';

    result.value = 0.0;
    result.error_bound = INFINITY;

    // Allocate empty strings to make cleanup safe
    result.value_str = (char*)flint_malloc(1);
    result.value_str[0] = '\0';
    result.error_str = (char*)flint_malloc(4);
    strcpy(result.error_str, "inf");
  }

  cleanup:
    acb_clear(acb_a);
  acb_clear(acb_b);
  acb_clear(acb_result);
  mag_clear(tol);

  return result;
}
