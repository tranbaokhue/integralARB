#ifndef INTEGRALARB_H
#define INTEGRALARB_H

#include <R.h>
#include <Rinternals.h>
#include "flint/flint.h"
#include "flint/arb.h"
#include "flint/acb.h"
#include "flint/acb_calc.h"

// Forward declaration for ACB calc function type
typedef int (*acb_calc_func_t)(acb_ptr, const acb_t, void*, slong, slong);

// Integration result structure
typedef struct {
  int evaluations;        // Number of function evaluations (estimate)
  int subdivisions;       // Number of subdivisions used (estimate)
  int status;            // Success/failure status (0 = success)
  char message[256];     // Status message
  // High-precision string representations (rigorous)
  char *value_str;       // High-precision value as string
  char *error_str;       // Rigorous error bound as string
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

// Built-in function entry structure
typedef struct {
  const char* name;
  const char* description;
  const char* latex_form;
  acb_calc_func_t func;
  const char* exact_antiderivative;
} builtin_function_entry_t;

// Function type definitions
typedef enum {
  INTEGRAND_EXPRESSION,   // String expression (current)
  INTEGRAND_R_FUNCTION,   // R function object (deprecated)
  INTEGRAND_BUILTIN      // Built-in optimized function
} integrand_type_t;

// Integrand structure
typedef struct {
  integrand_type_t type;
  void* data;              // Function-specific parameters
  SEXP r_function;         // R function (deprecated)
  acb_calc_func_t builtin_func; // Built-in function pointer
} integrand_t;

// Tokenizer types (needed for debugging)
typedef enum {
  TOKEN_NUMBER,
  TOKEN_IDENTIFIER,
  TOKEN_PLUS,
  TOKEN_MINUS,
  TOKEN_MULTIPLY,
  TOKEN_DIVIDE,
  TOKEN_POWER,
  TOKEN_LPAREN,
  TOKEN_RPAREN,
  TOKEN_EOF,
  TOKEN_ERROR
} token_type_t;

typedef struct {
  token_type_t type;
  char* value;
  double pos;
} token_t;

typedef struct {
  const char* input;
  int position;
  int length;
  token_t current_token;
} tokenizer_t;

// Expression parser types - exposed from builtin_functions.c
typedef enum {
  EXPR_NUMBER,
  EXPR_VARIABLE,      // x
  EXPR_FUNCTION,      // sin, cos, exp, etc.
  EXPR_BINARY_OP,     // +, -, *, /, ^
  EXPR_UNARY_MINUS,
  EXPR_CONSTANT
} expr_type_t;

typedef struct expression_node {
  expr_type_t type;
  char* value;                          // Number string or function name
  struct expression_node* left;
  struct expression_node* right;
  struct expression_node* argument;     // For functions
} expression_node_t;

typedef struct {
  arb_t coefficient;
  expression_node_t* expr;
  int sign;           // +1 or -1
} term_t;

typedef struct {
  term_t* terms;
  int num_terms;
  char* original_expression;
} parsed_expression_t;

// Tokenizer function declarations (for debugging)
tokenizer_t* create_tokenizer(const char* input);
void free_tokenizer(tokenizer_t* tok);
void next_token(tokenizer_t* tok);

// Function declarations
integration_result_t perform_integration(integrand_t* integrand,
                                         double a, double b,
                                         integration_params_t* params);

void cleanup_integrand(integrand_t* integrand);
void cleanup_integration_result(integration_result_t* result);

// Expression parser function declarations
parsed_expression_t* parse_mathematical_expression(const char* expression);
void cleanup_parsed_expression(parsed_expression_t* expr);
int parsed_expression_integrand(acb_ptr res, const acb_t z, void* param, slong order, slong prec);
int evaluate_expression_tree(acb_ptr result, expression_node_t* node,
                             const acb_t x, slong order, slong prec);

// Built-in function declarations
const builtin_function_entry_t* find_builtin_function(const char* name);
void list_builtin_functions(char* buffer, size_t buffer_size);

#endif // INTEGRALARB_H
