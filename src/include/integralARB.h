#ifndef INTEGRALARB_H
#define INTEGRALARB_H

#include <R.h>
#include <Rinternals.h>
#include "flint/flint.h"
#include "flint/arb.h"
#include "flint/acb.h"
#include "flint/acb_calc.h"

// Integration result structure
typedef struct {
  double value;           // Main result value
  double error_bound;     // Rigorous error bound
  double abs_error;       // Absolute error estimate
  double rel_error;       // Relative error estimate
  int evaluations;        // Number of function evaluations
  int subdivisions;       // Number of subdivisions used
  int status;             // Success/failure status (0 = success)
  char message[256];      // Status message
  // High-precision strings (decimal) to preserve precision when returning to R
  char *value_str;   // Midpoint value as string (precision ~= params->precision)
  char *error_str;   // Error bound (radius/ubound) as string
} integration_result_t;

// Integration parameters
typedef struct {
  int precision;          // Working precision in bits
  double abs_tolerance;   // Absolute tolerance
  double rel_tolerance;   // Relative tolerance
  int max_evaluations;    // Maximum function evaluations
  int max_depth;         // Maximum subdivision depth
  int verbose;           // Verbose output flag
} integration_params_t;

// Function type definitions
typedef enum {
  INTEGRAND_EXPRESSION,   // String expression
  INTEGRAND_R_FUNCTION,   // R function object
  INTEGRAND_BUILTIN      // Built-in optimized function
} integrand_type_t;

// Integrand structure
typedef struct {
  integrand_type_t type;
  void* data;            // Function-specific data
  SEXP r_function;       // R function (if type is R_FUNCTION)
} integrand_t;

// Function declarations
integration_result_t perform_integration(integrand_t* integrand,
                                         double a, double b,
                                         integration_params_t* params);

void cleanup_integrand(integrand_t* integrand);
void cleanup_integration_result(integration_result_t* result);

#endif // INTEGRALARB_H
