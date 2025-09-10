#include "include/integralARB.h"
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>

// Add these specific FLINT headers
#include "flint/flint.h"
#include "flint/arb.h"
#include "flint/acb.h"
#include "flint/arf.h"

// =============================================================================
// TOKENIZER AND PARSER STRUCTURES
// =============================================================================

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
  char* value;        // String representation for all tokens
  double pos;         // Position in input string
} token_t;

typedef struct {
  const char* input;
  int position;
  int length;
  token_t current_token;
} tokenizer_t;

// =============================================================================
// TOKENIZER IMPLEMENTATION
// =============================================================================

static void skip_whitespace(tokenizer_t* tok) {
  while (tok->position < tok->length && isspace(tok->input[tok->position])) {
    tok->position++;
  }
}

static char* extract_number(tokenizer_t* tok) {
  int start = tok->position;

  // Handle decimal numbers and scientific notation
  while (tok->position < tok->length &&
         (isdigit(tok->input[tok->position]) ||
         tok->input[tok->position] == '.' ||
         tok->input[tok->position] == 'e' ||
         tok->input[tok->position] == 'E' ||
         (tok->input[tok->position] == '-' && tok->position > start) ||
         (tok->input[tok->position] == '+' && tok->position > start))) {
    tok->position++;
  }

  int length = tok->position - start;
  char* result = malloc(length + 1);
  strncpy(result, &tok->input[start], length);
  result[length] = '\0';
  return result;
}

static char* extract_identifier(tokenizer_t* tok) {
  int start = tok->position;

  while (tok->position < tok->length &&
         (isalnum(tok->input[tok->position]) || tok->input[tok->position] == '_')) {
    tok->position++;
  }

  int length = tok->position - start;
  char* result = malloc(length + 1);
  strncpy(result, &tok->input[start], length);
  result[length] = '\0';
  return result;
}

static void next_token(tokenizer_t* tok) {
  if (tok->current_token.value) {
    free(tok->current_token.value);
    tok->current_token.value = NULL;
  }

  skip_whitespace(tok);

  if (tok->position >= tok->length) {
    tok->current_token.type = TOKEN_EOF;
    return;
  }

  char c = tok->input[tok->position];

  if (isdigit(c) || c == '.') {
    tok->current_token.type = TOKEN_NUMBER;
    tok->current_token.value = extract_number(tok);
  }
  else if (isalpha(c)) {
    tok->current_token.type = TOKEN_IDENTIFIER;
    tok->current_token.value = extract_identifier(tok);
  }
  else {
    tok->current_token.value = malloc(2);
    tok->current_token.value[0] = c;
    tok->current_token.value[1] = '\0';
    tok->position++;

    switch (c) {
    case '+': tok->current_token.type = TOKEN_PLUS; break;
    case '-': tok->current_token.type = TOKEN_MINUS; break;
    case '*': tok->current_token.type = TOKEN_MULTIPLY; break;
    case '/': tok->current_token.type = TOKEN_DIVIDE; break;
    case '^': tok->current_token.type = TOKEN_POWER; break;
    case '(': tok->current_token.type = TOKEN_LPAREN; break;
    case ')': tok->current_token.type = TOKEN_RPAREN; break;
    default: tok->current_token.type = TOKEN_ERROR; break;
    }
  }
}

static tokenizer_t* create_tokenizer(const char* input) {
  tokenizer_t* tok = malloc(sizeof(tokenizer_t));
  tok->input = input;
  tok->position = 0;
  tok->length = strlen(input);
  tok->current_token.value = NULL;
  next_token(tok);
  return tok;
}

static void free_tokenizer(tokenizer_t* tok) {
  if (tok->current_token.value) {
    free(tok->current_token.value);
  }
  free(tok);
}

// =============================================================================
// EXPRESSION TREE FUNCTIONS
// =============================================================================

static expression_node_t* create_node(expr_type_t type, const char* value) {
  expression_node_t* node = malloc(sizeof(expression_node_t));
  node->type = type;
  node->value = value ? strdup(value) : NULL;
  node->left = NULL;
  node->right = NULL;
  node->argument = NULL;
  return node;
}

static void free_expression_tree(expression_node_t* node) {
  if (!node) return;

  free(node->value);
  free_expression_tree(node->left);
  free_expression_tree(node->right);
  free_expression_tree(node->argument);
  free(node);
}

// =============================================================================
// RECURSIVE DESCENT PARSER
// =============================================================================

// Forward declarations
static expression_node_t* parse_expression(tokenizer_t* tok);
static expression_node_t* parse_term(tokenizer_t* tok);
static expression_node_t* parse_factor(tokenizer_t* tok);
static expression_node_t* parse_primary(tokenizer_t* tok);

static expression_node_t* parse_primary(tokenizer_t* tok) {
  if (tok->current_token.type == TOKEN_NUMBER) {
    expression_node_t* node = create_node(EXPR_NUMBER, tok->current_token.value);
    next_token(tok);
    return node;
  }

  if (tok->current_token.type == TOKEN_IDENTIFIER) {
    char* identifier = strdup(tok->current_token.value);
    next_token(tok);

    // Check if it's a function call
    if (tok->current_token.type == TOKEN_LPAREN) {
      next_token(tok); // consume '('
      expression_node_t* arg = parse_expression(tok);

      if (tok->current_token.type != TOKEN_RPAREN) {
        free(identifier);
        free_expression_tree(arg);
        return NULL; // Error: expected ')'
      }
      next_token(tok); // consume ')'

      expression_node_t* func_node = create_node(EXPR_FUNCTION, identifier);
      func_node->argument = arg;
      free(identifier);
      return func_node;
    }
    else if (strcmp(identifier, "x") == 0) {
      expression_node_t* var_node = create_node(EXPR_VARIABLE, identifier);
      free(identifier);
      return var_node;
    }
    else if (strcmp(identifier, "pi") == 0) {
      expression_node_t* const_node = create_node(EXPR_CONSTANT, "pi");
      free(identifier);
      return const_node;
    } else if (strcmp(identifier, "e") == 0) {
      expression_node_t* const_node = create_node(EXPR_CONSTANT, "e");
      free(identifier);
      return const_node;
    } else if (strcmp(identifier, "ln") == 0) {
      // Handle ln(number) as a constant evaluation
      if (tok->current_token.type == TOKEN_LPAREN) {
        next_token(tok); // consume '('
        if (tok->current_token.type == TOKEN_NUMBER) {
          char* number = strdup(tok->current_token.value);
          next_token(tok);
          if (tok->current_token.type == TOKEN_RPAREN) {
            next_token(tok); // consume ')'
            char* ln_expr = malloc(strlen(number) + 10);
            sprintf(ln_expr, "ln_%s", number);
            expression_node_t* ln_node = create_node(EXPR_CONSTANT, ln_expr);
            free(number);
            free(ln_expr);
            free(identifier);
            return ln_node;
          }
        }
      }
    }
    else {
      free(identifier);
      return NULL; // Error: unknown identifier
    }
  }

  if (tok->current_token.type == TOKEN_LPAREN) {
    next_token(tok); // consume '('
    expression_node_t* node = parse_expression(tok);
    if (tok->current_token.type != TOKEN_RPAREN) {
      free_expression_tree(node);
      return NULL; // Error: expected ')'
    }
    next_token(tok); // consume ')'
    return node;
  }

  if (tok->current_token.type == TOKEN_MINUS) {
    next_token(tok);
    expression_node_t* operand = parse_primary(tok);
    if (!operand) return NULL;

    expression_node_t* unary_node = create_node(EXPR_UNARY_MINUS, NULL);
    unary_node->right = operand;
    return unary_node;
  }

  return NULL; // Error
}

static expression_node_t* parse_factor(tokenizer_t* tok) {
  expression_node_t* left = parse_primary(tok);
  if (!left) return NULL;

  while (tok->current_token.type == TOKEN_POWER) {
    char* op = strdup(tok->current_token.value);
    next_token(tok);
    expression_node_t* right = parse_primary(tok);
    if (!right) {
      free(op);
      free_expression_tree(left);
      return NULL;
    }

    expression_node_t* op_node = create_node(EXPR_BINARY_OP, op);
    op_node->left = left;
    op_node->right = right;
    left = op_node;
    free(op);
  }

  return left;
}

static expression_node_t* parse_term(tokenizer_t* tok) {
  expression_node_t* left = parse_factor(tok);
  if (!left) return NULL;

  while (tok->current_token.type == TOKEN_MULTIPLY || tok->current_token.type == TOKEN_DIVIDE) {
    char* op = strdup(tok->current_token.value);
    next_token(tok);
    expression_node_t* right = parse_factor(tok);
    if (!right) {
      free(op);
      free_expression_tree(left);
      return NULL;
    }

    expression_node_t* op_node = create_node(EXPR_BINARY_OP, op);
    op_node->left = left;
    op_node->right = right;
    left = op_node;
    free(op);
  }

  return left;
}

static expression_node_t* parse_expression(tokenizer_t* tok) {
  expression_node_t* left = parse_term(tok);
  if (!left) return NULL;

  while (tok->current_token.type == TOKEN_PLUS || tok->current_token.type == TOKEN_MINUS) {
    char* op = strdup(tok->current_token.value);
    next_token(tok);
    expression_node_t* right = parse_term(tok);
    if (!right) {
      free(op);
      free_expression_tree(left);
      return NULL;
    }

    expression_node_t* op_node = create_node(EXPR_BINARY_OP, op);
    op_node->left = left;
    op_node->right = right;
    left = op_node;
    free(op);
  }

  return left;
}

// =============================================================================
// SUPPORTED FUNCTION VALIDATION
// =============================================================================

static int is_supported_function(const char* name) {
  const char* supported[] = {
    "sin", "cos", "exp", "log", "ln", "sinh", "cosh", "atan", "arctan", NULL
  };

  for (int i = 0; supported[i] != NULL; i++) {
    if (strcmp(name, supported[i]) == 0) {
      return 1;
    }
  }
  return 0;
}

static int validate_expression_tree(expression_node_t* node) {
  if (!node) return 0;

  switch (node->type) {
  case EXPR_NUMBER:
  case EXPR_VARIABLE:
    return 1;

  case EXPR_FUNCTION:
    if (!is_supported_function(node->value)) {
      return 0;
    }
    return validate_expression_tree(node->argument);

  case EXPR_BINARY_OP:
    return validate_expression_tree(node->left) &&
      validate_expression_tree(node->right);

  case EXPR_UNARY_MINUS:
    return validate_expression_tree(node->right);

  default:
    return 0;
  }
}

// =============================================================================
// ARB EVALUATION ENGINE
// =============================================================================

int evaluate_expression_tree(acb_ptr result, expression_node_t* node,
                                    const acb_t x, slong order, slong prec);

static int evaluate_function_derivative(acb_ptr result, const char* func_name, const acb_t arg, slong order, slong prec) {
  if (strcmp(func_name, "sin") == 0) {
    switch (order % 4) {
    case 0: acb_sin(result, arg, prec); break;
    case 1: acb_cos(result, arg, prec); break;
    case 2: acb_sin(result, arg, prec); acb_neg(result, result); break;
    case 3: acb_cos(result, arg, prec); acb_neg(result, result); break;
    }
    return 1;
  }
  else if (strcmp(func_name, "cos") == 0) {
    switch (order % 4) {
    case 0: acb_cos(result, arg, prec); break;
    case 1: acb_sin(result, arg, prec); acb_neg(result, result); break;
    case 2: acb_cos(result, arg, prec); acb_neg(result, result); break;
    case 3: acb_sin(result, arg, prec); break;
    }
    return 1;
  }
  else if (strcmp(func_name, "atan") == 0 || strcmp(func_name, "arctan") == 0) {
    if (order == 0) {
      acb_atan(result, arg, prec);
    } else if (order == 1) {
      // d/dx[atan(x)] = 1/(1+x²)
      acb_t temp;
      acb_init(temp);
      acb_mul(temp, arg, arg, prec);
      acb_add_ui(temp, temp, 1, prec);
      acb_inv(result, temp, prec);
      acb_clear(temp);
    } else {
      // Higher derivatives: d^n/dx^n[atan(x)] = derivatives of 1/(1+x²)
      acb_indeterminate(result);
    }
    return 1;
  }
  else if (strcmp(func_name, "exp") == 0) {
    // All derivatives of exp(x) are exp(x)
    acb_exp(result, arg, prec);
    return 1;
  }
  else if (strcmp(func_name, "log") == 0 || strcmp(func_name, "ln") == 0) {
    if (order == 0) {
      acb_log(result, arg, prec);
    } else if (order == 1) {
      // d/dx[log(x)] = 1/x
      acb_inv(result, arg, prec);
    } else {
      // d^n/dx^n[log(x)] = (-1)^(n-1) * (n-1)! / x^n
      acb_t temp;
      acb_init(temp);
      acb_pow_ui(temp, arg, order, prec);
      acb_inv(temp, temp, prec);

      slong factorial = 1;
      for (slong i = 1; i < order; i++) factorial *= i;

      acb_mul_si(result, temp, ((order - 1) % 2 == 0) ? factorial : -factorial, prec);
      acb_clear(temp);
    }
    return 1;
  }
  else if (strcmp(func_name, "sinh") == 0) {
    // Derivatives alternate: sinh, cosh, sinh, cosh, ...
    if (order % 2 == 0) {
      acb_sinh(result, arg, prec);
    } else {
      acb_cosh(result, arg, prec);
    }
    return 1;
  }
  else if (strcmp(func_name, "cosh") == 0) {
    // Derivatives alternate: cosh, sinh, cosh, sinh, ...
    if (order % 2 == 0) {
      acb_cosh(result, arg, prec);
    } else {
      acb_sinh(result, arg, prec);
    }
    return 1;
  }

  return 0; // Unsupported function
}

int evaluate_expression_tree(acb_ptr result, expression_node_t* node,
                                    const acb_t x, slong order, slong prec) {
  if (!node) return 0;

  switch (node->type) {
  case EXPR_NUMBER: {
    if (order == 0) {
    arb_t real_part;
    arb_init(real_part);

    if (arb_set_str(real_part, node->value, prec) != 0) {
      arb_clear(real_part);
      acb_indeterminate(result);
      return 0;
    }

    acb_set_arb(result, real_part);
    arb_clear(real_part);
  } else {
    acb_zero(result);
  }
  return 1;
  }

  case EXPR_VARIABLE: {
    if (order == 0) {
    acb_set(result, x);
  } else if (order == 1) {
    acb_one(result);
  } else {
    acb_zero(result);
  }
  return 1;
  }

  case EXPR_FUNCTION: {
    if (order == 0) {
    // Function evaluation
    acb_t arg_val;
    acb_init(arg_val);
    if (!evaluate_expression_tree(arg_val, node->argument, x, 0, prec)) {
      acb_clear(arg_val);
      return 0;
    }
    int success = evaluate_function_derivative(result, node->value, arg_val, 0, prec);
    acb_clear(arg_val);
    return success;
  }
    else if (order == 1) {
      // Case 1: f(c*x) where c is constant
      if (node->argument->type == EXPR_BINARY_OP &&
          strcmp(node->argument->value, "*") == 0) {

        // Check both orders: c*x or x*c
        expression_node_t* coeff_node = NULL;
        expression_node_t* var_node = NULL;

        if (node->argument->left->type == EXPR_NUMBER &&
            node->argument->right->type == EXPR_VARIABLE) {
          coeff_node = node->argument->left;
          var_node = node->argument->right;
        } else if (node->argument->left->type == EXPR_VARIABLE &&
          node->argument->right->type == EXPR_NUMBER) {
          coeff_node = node->argument->right;
          var_node = node->argument->left;
        }

        if (coeff_node && var_node) {
          arb_t coeff;
          acb_t scaled_x;
          arb_init(coeff);
          acb_init(scaled_x);

          arb_set_str(coeff, coeff_node->value, prec);
          acb_mul_arb(scaled_x, x, coeff, prec);

          // f'(c*x) * c
          evaluate_function_derivative(result, node->value, scaled_x, 1, prec);
          acb_mul_arb(result, result, coeff, prec);

          arb_clear(coeff);
          acb_clear(scaled_x);
          return 1;
        }
      }

      // Case 2: f(-x)
      else if (node->argument->type == EXPR_UNARY_MINUS &&
               node->argument->right->type == EXPR_VARIABLE) {

        acb_t neg_x;
        acb_init(neg_x);
        acb_neg(neg_x, x);

        evaluate_function_derivative(result, node->value, neg_x, 1, prec);
        acb_neg(result, result);  // f'(-x) * (-1)

        acb_clear(neg_x);
        return 1;
      }

      // Case 3: f(x) - simple variable
      else if (node->argument->type == EXPR_VARIABLE) {
        return evaluate_function_derivative(result, node->value, x, 1, prec);
      }

      // All other cases: NOT SUPPORTED
      else {
        acb_indeterminate(result);
        return 0;
      }
    }
    else {
      acb_indeterminate(result);
      return 0;
    }
  }

  case EXPR_BINARY_OP: {
    if (strcmp(node->value, "+") == 0) {
    acb_t left_val, right_val;
    acb_init(left_val); acb_init(right_val);

    if (!evaluate_expression_tree(left_val, node->left, x, order, prec) ||
        !evaluate_expression_tree(right_val, node->right, x, order, prec)) {
        acb_clear(left_val); acb_clear(right_val);
        return 0;
    }

    acb_add(result, left_val, right_val, prec);
    acb_clear(left_val); acb_clear(right_val);
    return 1;
  }
    else if (strcmp(node->value, "-") == 0) {
      acb_t left_val, right_val;
      acb_init(left_val); acb_init(right_val);

      if (!evaluate_expression_tree(left_val, node->left, x, order, prec) ||
          !evaluate_expression_tree(right_val, node->right, x, order, prec)) {
          acb_clear(left_val); acb_clear(right_val);
          return 0;
      }

      acb_sub(result, left_val, right_val, prec);
      acb_clear(left_val); acb_clear(right_val);
      return 1;
    }
    else if (strcmp(node->value, "*") == 0) {
      if (order == 0) {
        // Multiplication: (f * g)
        acb_t left_val, right_val;
        acb_init(left_val); acb_init(right_val);

        if (!evaluate_expression_tree(left_val, node->left, x, 0, prec) ||
            !evaluate_expression_tree(right_val, node->right, x, 0, prec)) {
            acb_clear(left_val); acb_clear(right_val);
            return 0;
        }

        acb_mul(result, left_val, right_val, prec);
        acb_clear(left_val); acb_clear(right_val);
        return 1;
      }
      else if (strcmp(node->value, "^") == 0) {
        if (node->left->type == EXPR_VARIABLE && node->right->type == EXPR_NUMBER) {
          // Handle x^n case
          if (order == 0) {
            acb_t exp_val;
            acb_init(exp_val);

            arb_t exp_arb;
            arb_init(exp_arb);
            arb_set_str(exp_arb, node->right->value, prec);
            acb_set_arb(exp_val, exp_arb);

            acb_pow(result, x, exp_val, prec);  // x^n

            arb_clear(exp_arb);
            acb_clear(exp_val);
            return 1;
          }
          else if (order == 1) {
            // Power rule: d/dx[x^n] = n*x^(n-1)
            arb_t n, n_minus_1;
            acb_t x_power;
            arb_init(n); arb_init(n_minus_1); acb_init(x_power);

            arb_set_str(n, node->right->value, prec);
            arb_sub_ui(n_minus_1, n, 1, prec);

            acb_pow_arb(x_power, x, n_minus_1, prec);  // x^(n-1)
            acb_mul_arb(result, x_power, n, prec);     // n*x^(n-1)

            arb_clear(n); arb_clear(n_minus_1); acb_clear(x_power);
            return 1;
          }
        }
      }
      else if (strcmp(node->value, "/") == 0) {
        if (order == 0) {
          // General division: a/b
          acb_t left_val, right_val;
          acb_init(left_val); acb_init(right_val);

          if (!evaluate_expression_tree(left_val, node->left, x, 0, prec) ||
              !evaluate_expression_tree(right_val, node->right, x, 0, prec)) {
              acb_clear(left_val); acb_clear(right_val);
              return 0;
          }

          acb_div(result, left_val, right_val, prec);
          acb_clear(left_val); acb_clear(right_val);
          return 1;
        }
        else if (order == 1) {
          // Only 1/x derivative supported
          if (node->left->type == EXPR_NUMBER &&
              strcmp(node->left->value, "1") == 0 &&
              node->right->type == EXPR_VARIABLE) {

            acb_t x_squared;
            acb_init(x_squared);
            acb_mul(x_squared, x, x, prec);
            acb_inv(result, x_squared, prec);
            acb_neg(result, result);
            acb_clear(x_squared);
            return 1;
          } else {
            acb_indeterminate(result);
            return 0;
          }
        }
      }

      else if (order == 1) {
        // SPECIAL CASE: c*f(x) where c is constant
        if (node->left->type == EXPR_NUMBER &&
            node->right->type == EXPR_FUNCTION) {

          arb_t coeff;
          acb_t func_derivative;
          arb_init(coeff);
          acb_init(func_derivative);

          arb_set_str(coeff, node->left->value, prec);

          if (!evaluate_expression_tree(func_derivative, node->right, x, 1, prec)) {
            arb_clear(coeff);
            acb_clear(func_derivative);
            return 0;
          }

          acb_mul_arb(result, func_derivative, coeff, prec);

          arb_clear(coeff);
          acb_clear(func_derivative);
          return 1;
        }

        // SPECIAL CASE: f(x)*c where c is constant
        else if (node->left->type == EXPR_FUNCTION &&
                 node->right->type == EXPR_NUMBER) {

          arb_t coeff;
          acb_t func_derivative;
          arb_init(coeff);
          acb_init(func_derivative);

          arb_set_str(coeff, node->right->value, prec);

          if (!evaluate_expression_tree(func_derivative, node->left, x, 1, prec)) {
            arb_clear(coeff);
            acb_clear(func_derivative);
            return 0;
          }

          acb_mul_arb(result, func_derivative, coeff, prec);

          arb_clear(coeff);
          acb_clear(func_derivative);
          return 1;
        }

        // GENERAL CASE: Product rule (f*g)' = f'*g + f*g'
        else {
          acb_t f, g, f_prime, g_prime, term1, term2;
          acb_init(f); acb_init(g); acb_init(f_prime); acb_init(g_prime);
          acb_init(term1); acb_init(term2);

          if (!evaluate_expression_tree(f, node->left, x, 0, prec) ||
              !evaluate_expression_tree(g, node->right, x, 0, prec) ||
              !evaluate_expression_tree(f_prime, node->left, x, 1, prec) ||
              !evaluate_expression_tree(g_prime, node->right, x, 1, prec)) {
              acb_clear(f); acb_clear(g); acb_clear(f_prime); acb_clear(g_prime);
              acb_clear(term1); acb_clear(term2);
              return 0;
          }

          acb_mul(term1, f_prime, g, prec);    // f' * g
          acb_mul(term2, f, g_prime, prec);    // f * g'
          acb_add(result, term1, term2, prec); // f'*g + f*g'

          acb_clear(f); acb_clear(g); acb_clear(f_prime); acb_clear(g_prime);
          acb_clear(term1); acb_clear(term2);
          return 1;
        }
      }
      else {
        // Higher order multiplication derivatives not implemented
        acb_indeterminate(result);
        return 0;
      }
    }
    else if (strcmp(node->value, "^") == 0) {
      // Handle simple cases like x^n where n is a constant
      if (order == 0) {
        acb_t base_val, exp_val;
        acb_init(base_val); acb_init(exp_val);

        if (!evaluate_expression_tree(base_val, node->left, x, 0, prec) ||
            !evaluate_expression_tree(exp_val, node->right, x, 0, prec)) {
            acb_clear(base_val); acb_clear(exp_val);
            return 0;
        }

        acb_pow(result, base_val, exp_val, prec);
        acb_clear(base_val); acb_clear(exp_val);
        return 1;
      } else {
        // Power rule derivatives are complex - not implemented
        acb_indeterminate(result);
        return 0;
      }
    }
    else {
      // Other operators not supported
      acb_indeterminate(result);
      return 0;
    }
  }

  case EXPR_UNARY_MINUS: {
    if (!evaluate_expression_tree(result, node->right, x, order, prec)) {
    return 0;
  }
    acb_neg(result, result);
    return 1;
  }

  case EXPR_CONSTANT: {
    if (order == 0) {
    if (strcmp(node->value, "pi") == 0) {
      arb_t pi_val;
      arb_init(pi_val);
      arb_const_pi(pi_val, prec);
      acb_set_arb(result, pi_val);
      arb_clear(pi_val);
    } else if (strcmp(node->value, "e") == 0) {
      arb_t e_val;
      arb_init(e_val);
      arb_const_e(e_val, prec);
      acb_set_arb(result, e_val);
      arb_clear(e_val);
    } else if (strncmp(node->value, "ln_", 3) == 0) {
      // Handle ln(number)
      arb_t num_val, ln_val;
      arb_init(num_val); arb_init(ln_val);
      arb_set_str(num_val, &node->value[3], prec); // Skip "ln_"
      arb_log(ln_val, num_val, prec);
      acb_set_arb(result, ln_val);
      arb_clear(num_val); arb_clear(ln_val);
    }
  } else {
    acb_zero(result); // Constants have zero derivative
  }
  return 1;
  }

  default:
    return 0;
  }
}

// =============================================================================
// MAIN INTEGRAND FUNCTION FOR PARSED EXPRESSIONS
// =============================================================================

int parsed_expression_integrand(acb_ptr res, const acb_t z, void* param, slong order, slong prec) {
  expression_node_t* expr_tree = (expression_node_t*)param;

  if (!evaluate_expression_tree(res, expr_tree, z, order, prec)) {
    acb_indeterminate(res);
    return -1;
  }

  return 0;
}

// =============================================================================
// PUBLIC INTERFACE FUNCTIONS
// =============================================================================

parsed_expression_t* parse_mathematical_expression(const char* expression) {
  tokenizer_t* tok = create_tokenizer(expression);
  expression_node_t* tree = parse_expression(tok);

  if (!tree || !validate_expression_tree(tree)) {
    free_tokenizer(tok);
    if (tree) free_expression_tree(tree);
    return NULL;
  }

  parsed_expression_t* result = malloc(sizeof(parsed_expression_t));
  result->terms = malloc(sizeof(term_t));
  result->num_terms = 1;
  result->original_expression = strdup(expression);

  // For now, treat entire expression as single term
  arb_init(result->terms[0].coefficient);
  arb_one(result->terms[0].coefficient);
  result->terms[0].expr = tree;
  result->terms[0].sign = 1;

  free_tokenizer(tok);
  return result;
}

void cleanup_parsed_expression(parsed_expression_t* expr) {
  if (!expr) return;

  for (int i = 0; i < expr->num_terms; i++) {
    arb_clear(expr->terms[i].coefficient);
    free_expression_tree(expr->terms[i].expr);
  }

  free(expr->terms);
  free(expr->original_expression);
  free(expr);
}

const builtin_function_entry_t* find_builtin_function(const char* name) {
  // This is now mainly for validation - actual evaluation uses the parser
  const char* supported[] = {
    "sin", "cos", "exp", "log", "ln", "sinh", "cosh", "atan", "arctan", NULL
  };

  for (int i = 0; supported[i] != NULL; i++) {
    if (strcmp(name, supported[i]) == 0) {
      // Return a dummy entry - not used for parsed expressions
      static builtin_function_entry_t dummy = {"parsed", "Parsed expression", "expr", NULL, "N/A"};
      return &dummy;
    }
  }
  return NULL;
}

void list_builtin_functions(char* buffer, size_t buffer_size) {
  snprintf(buffer, buffer_size,
           "Supported mathematical expressions:\n"
           "Variables: x\n"
           "Functions: sin(expr), cos(expr), exp(expr), log(expr), ln(expr), sinh(expr), cosh(expr), atan(expr)\n"
           "  where expr can be: x, 2*x, -x, etc.\n"
           "Operations: +, -, *, ^, ()\n"
           "Numbers: Any decimal number (use strings for high precision)\n"
           "Examples: 'sin(2*x)', 'cos(x) + exp(-x)', 'atan(x)', '2*sinh(x) - x'\n"
           "Chain rule: First derivatives supported for f(g(x)) where g(x) is simple\n");
}
