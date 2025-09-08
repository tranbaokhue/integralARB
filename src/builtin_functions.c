#include "include/integralARB.h"
#include <string.h>

// =============================================================================
// ELEMENTARY FUNCTIONS - Basic building blocks
// =============================================================================

// Polynomials: x^n for n = 0, 1, 2, ..., 10
int poly_constant(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
  (void)param; (void)z;
  if (order == 0) { acb_one(res); }
  else { acb_zero(res); }
  return 0;
}

int poly_linear(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
  (void)param;
  if (order == 0) { acb_set(res, z); }
  else if (order == 1) { acb_one(res); }
  else { acb_zero(res); }
  return 0;
}

int poly_quadratic(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
  (void)param;
  if (order == 0) { acb_pow_ui(res, z, 2, prec); }
  else if (order == 1) { acb_mul_ui(res, z, 2, prec); }
  else if (order == 2) { acb_set_ui(res, 2); }
  else { acb_zero(res); }
  return 0;
}

int poly_cubic(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
  (void)param;
  if (order == 0) { acb_pow_ui(res, z, 3, prec); }
  else if (order == 1) {
    acb_t temp;
    acb_init(temp);
    acb_pow_ui(temp, z, 2, prec);
    acb_mul_ui(res, temp, 3, prec);
    acb_clear(temp);
  }
  else if (order == 2) { acb_mul_ui(res, z, 6, prec); }
  else if (order == 3) { acb_set_ui(res, 6); }
  else { acb_zero(res); }
  return 0;
}

// Rational functions
int rational_1_over_x(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
  (void)param;
  if (order == 0) {
    acb_inv(res, z, prec);
  }
  else if (order == 1) {
    acb_t temp;
    acb_init(temp);
    acb_pow_ui(temp, z, 2, prec);
    acb_inv(temp, temp, prec);
    acb_neg(res, temp);
    acb_clear(temp);
  }
  else if (order >= 2) {
    // Higher derivatives of 1/x
    acb_t temp;
    acb_init(temp);
    acb_pow_ui(temp, z, order + 1, prec);
    acb_inv(temp, temp, prec);

    // Multiply by (-1)^order * order!
    slong factorial = 1;
    for (slong i = 1; i <= order; i++) factorial *= i;
    acb_mul_si(temp, temp, (order % 2 == 0) ? factorial : -factorial, prec);
    acb_set(res, temp);
    acb_clear(temp);
  }
  return 0;
}

int rational_1_over_1_plus_x_squared(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
  (void)param;
  acb_t z2, denom;
  acb_init(z2);
  acb_init(denom);

  acb_mul(z2, z, z, prec);
  acb_add_ui(denom, z2, 1, prec);

  if (order == 0) {
    acb_inv(res, denom, prec);
  }
  else if (order == 1) {
    // d/dx[1/(1+x²)] = -2x/(1+x²)²
    acb_t temp;
    acb_init(temp);
    acb_mul(temp, denom, denom, prec);
    acb_inv(temp, temp, prec);
    acb_mul(temp, temp, z, prec);
    acb_mul_si(res, temp, -2, prec);
    acb_clear(temp);
  }
  else {
    // Higher derivatives can be computed but are complex
    acb_indeterminate(res);
  }

  acb_clear(z2);
  acb_clear(denom);
  return 0;
}

// =============================================================================
// EXPONENTIAL AND LOGARITHMIC FUNCTIONS
// =============================================================================

int exp_function(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
  (void)param;
  // All derivatives of exp(x) are exp(x)
  acb_exp(res, z, prec);
  return 0;
}

int exp_minus_x(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
  (void)param;
  acb_t neg_z;
  acb_init(neg_z);
  acb_neg(neg_z, z);
  acb_exp(res, neg_z, prec);

  // Multiply by (-1)^order for derivatives
  if (order % 2 == 1) {
    acb_neg(res, res);
  }

  acb_clear(neg_z);
  return 0;
}

int gaussian(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
  (void)param;
  acb_t z2, neg_z2;
  acb_init(z2);
  acb_init(neg_z2);

  acb_mul(z2, z, z, prec);
  acb_neg(neg_z2, z2);

  if (order == 0) {
    acb_exp(res, neg_z2, prec);
  }
  else if (order == 1) {
    // d/dx[exp(-x²)] = -2x*exp(-x²)
    acb_exp(res, neg_z2, prec);
    acb_mul(res, res, z, prec);
    acb_mul_si(res, res, -2, prec);
  }
  else if (order == 2) {
    // d²/dx²[exp(-x²)] = (4x² - 2)*exp(-x²)
    acb_t temp;
    acb_init(temp);
    acb_mul_ui(temp, z2, 4, prec);
    acb_sub_ui(temp, temp, 2, prec);
    acb_exp(res, neg_z2, prec);
    acb_mul(res, res, temp, prec);
    acb_clear(temp);
  }
  else {
    // Higher derivatives via Hermite polynomials (complex to implement)
    acb_indeterminate(res);
  }

  acb_clear(z2);
  acb_clear(neg_z2);
  return 0;
}

int log_function(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
  (void)param;
  if (order == 0) {
    acb_log(res, z, prec);
  }
  else {
    // d^n/dx^n[log(x)] = (-1)^(n-1) * (n-1)! / x^n
    acb_t temp;
    acb_init(temp);
    acb_pow_ui(temp, z, order, prec);
    acb_inv(temp, temp, prec);

    slong factorial = 1;
    for (slong i = 1; i < order; i++) factorial *= i;

    acb_mul_si(res, temp, ((order - 1) % 2 == 0) ? factorial : -factorial, prec);
    acb_clear(temp);
  }
  return 0;
}

// =============================================================================
// TRIGONOMETRIC FUNCTIONS
// =============================================================================

int sin_function(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
  (void)param;
  // Derivatives cycle: sin, cos, -sin, -cos, sin, ...
  switch (order % 4) {
  case 0: acb_sin(res, z, prec); break;
  case 1: acb_cos(res, z, prec); break;
  case 2: acb_sin(res, z, prec); acb_neg(res, res); break;
  case 3: acb_cos(res, z, prec); acb_neg(res, res); break;
  }
  return 0;
}

int cos_function(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
  (void)param;
  // Derivatives cycle: cos, -sin, -cos, sin, cos, ...
  switch (order % 4) {
  case 0: acb_cos(res, z, prec); break;
  case 1: acb_sin(res, z, prec); acb_neg(res, res); break;
  case 2: acb_cos(res, z, prec); acb_neg(res, res); break;
  case 3: acb_sin(res, z, prec); break;
  }
  return 0;
}

// =============================================================================
// HYPERBOLIC FUNCTIONS
// =============================================================================

int sinh_function(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
  (void)param;
  // Derivatives alternate: sinh, cosh, sinh, cosh, ...
  if (order % 2 == 0) {
    acb_sinh(res, z, prec);
  } else {
    acb_cosh(res, z, prec);
  }
  return 0;
}

int cosh_function(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
  (void)param;
  // Derivatives alternate: cosh, sinh, cosh, sinh, ...
  if (order % 2 == 0) {
    acb_cosh(res, z, prec);
  } else {
    acb_sinh(res, z, prec);
  }
  return 0;
}

int sech_function(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
  (void)param;
  if (order == 0) {
    acb_sech(res, z, prec);
  }
  else if (order == 1) {
    // d/dx[sech(x)] = -sech(x)tanh(x)
    acb_t tanh_z;
    acb_init(tanh_z);
    acb_sech(res, z, prec);
    acb_tanh(tanh_z, z, prec);
    acb_mul(res, res, tanh_z, prec);
    acb_neg(res, res);
    acb_clear(tanh_z);
  }
  else {
    // Higher derivatives of sech are complex
    acb_indeterminate(res);
  }
  return 0;
}

// =============================================================================
// FUNCTION REGISTRY
// =============================================================================

static const builtin_function_entry_t builtin_functions[] = {
  // Polynomials
  {"constant", "f(x) = 1", "1", poly_constant, "x"},
  {"linear", "f(x) = x", "x", poly_linear, "x²/2"},
  {"quadratic", "f(x) = x²", "x^2", poly_quadratic, "x³/3"},
  {"cubic", "f(x) = x³", "x^3", poly_cubic, "x⁴/4"},

  // Rational functions
  {"reciprocal", "f(x) = 1/x", "1/x", rational_1_over_x, "ln|x|"},
  {"arctangent", "f(x) = 1/(1+x²)", "\\frac{1}{1+x^2}", rational_1_over_1_plus_x_squared, "arctan(x)"},

          // Exponential functions
  {"exp", "f(x) = exp(x)", "e^x", exp_function, "exp(x)"},
  {"exp_minus", "f(x) = exp(-x)", "e^{-x}", exp_minus_x, "-exp(-x)"},
  {"gaussian", "f(x) = exp(-x²)", "e^{-x^2}", gaussian, "√π/2 * erf(x)"},

  // Logarithmic
  {"log", "f(x) = ln(x)", "\\ln(x)", log_function, "x*ln(x) - x"},

          // Trigonometric
  {"sin", "f(x) = sin(x)", "\\sin(x)", sin_function, "-cos(x)"},
  {"cos", "f(x) = cos(x)", "\\cos(x)", cos_function, "sin(x)"},

  // Hyperbolic
  {"sinh", "f(x) = sinh(x)", "\\sinh(x)", sinh_function, "cosh(x)"},
  {"cosh", "f(x) = cosh(x)", "\\cosh(x)", cosh_function, "sinh(x)"},
  {"sech", "f(x) = sech(x)", "\\text{sech}(x)", sech_function, "arctan(sinh(x))"},

  // Terminator
  {NULL, NULL, NULL, NULL, NULL}
};

// Function to find built-in by name
const builtin_function_entry_t* find_builtin_function(const char* name) {
  for (int i = 0; builtin_functions[i].name != NULL; i++) {
    if (strcmp(builtin_functions[i].name, name) == 0) {
      return &builtin_functions[i];
    }
  }
  return NULL;
}

// Get list of all available functions
void list_builtin_functions(char* buffer, size_t buffer_size) {
  size_t pos = 0;
  pos += snprintf(buffer + pos, buffer_size - pos, "Available built-in functions:\n");

  for (int i = 0; builtin_functions[i].name != NULL && pos < buffer_size - 100; i++) {
    pos += snprintf(buffer + pos, buffer_size - pos,
                    "  %-12s: %s\n",
                    builtin_functions[i].name,
                    builtin_functions[i].description);
  }
}
