// src/include/integralARB.h - Updated header with built-in function support

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
  double value;           // Double approximation of result
  double error_bound;     // Double approximation of error bound
  double abs_error;       // Absolute error estimate
  double rel_error;       // Relative error estimate
  int evaluations;        // Number of function evaluations
  int subdivisions;       // Number of subdivisions used
  int status;            // Success/failure status (0 = success)
  char message[256];     // Status message
  // High-precision string representations
  char *value_str;       // High-precision value as string
  char *error_str;       // High-precision error bound as string
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
  INTEGRAND_EXPRESSION,   // String expression (future)
  INTEGRAND_R_FUNCTION,   // R function object (deprecated)
  INTEGRAND_BUILTIN      // Built-in optimized function
} integrand_type_t;

// Built-in function entry structure
typedef struct {
  const char* name;
  const char* description;
  const char* latex_form;
  acb_calc_func_t func;
  const char* exact_antiderivative;
} builtin_function_entry_t;

// Integrand structure
typedef struct {
  integrand_type_t type;
  void* data;              // Function-specific parameters
  SEXP r_function;         // R function (deprecated)
  acb_calc_func_t builtin_func; // Built-in function pointer
} integrand_t;

// Function declarations
integration_result_t perform_integration(integrand_t* integrand,
                                         double a, double b,
                                         integration_params_t* params);

void cleanup_integrand(integrand_t* integrand);
void cleanup_integration_result(integration_result_t* result);

// Built-in function declarations
extern const builtin_function_entry_t* find_builtin_function(const char* name);
extern void list_builtin_functions(char* buffer, size_t buffer_size);

#endif // INTEGRALARB_H
