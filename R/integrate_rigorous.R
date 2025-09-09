#' Rigorous Numerical Integration with Mathematical Expression Parsing
#'
#' Computes definite integrals using ball arithmetic with mathematically
#' rigorous error bounds. Supports mathematical expressions with automatic
#' parsing and derivative computation for guaranteed accuracy.
#'
#' @param expression Mathematical expression as a character string (required)
#' @param a Lower limit of integration as character string for precision (required)
#' @param b Upper limit of integration as character string for precision (required)
#' @param precision Working precision in bits (default: 128, range: 64-1024)
#' @param abs_tol Absolute error tolerance (default: 2^(-precision-10))
#' @param rel_tol Relative error tolerance (default: 2^(-precision-10))
#' @param max_eval Maximum number of function evaluations (default: 1000000)
#' @param max_depth Maximum subdivision depth (default: 100)
#' @param verbose Print integration progress (default: FALSE)
#'
#' @return An object of class `rigorous_result` containing:
#' \itemize{
#'   \item `value` - Double approximation of integral value
#'   \item `value_str` - High-precision string representation
#'   \item `error_bound` - Rigorous error bound (double)
#'   \item `error_str` - High-precision error bound (string)
#'   \item `abs_error` - Absolute error estimate
#'   \item `rel_error` - Relative error estimate
#'   \item `evaluations` - Number of function evaluations performed
#'   \item `status` - Integration status (0 = success)
#'   \item `message` - Status message
#'   \item `expression` - Original expression string
#'   \item `precision` - Working precision used
#' }
#'
#' @details
#' **Supported Mathematical Expressions:**
#'
#' *Variables:*
#' - `x` - Integration variable
#'
#' *Functions (with chain rule support):*
#' - `sin(expr)` - Sine function
#' - `cos(expr)` - Cosine function
#' - `exp(expr)` - Exponential function
#' - `log(expr)`, `ln(expr)` - Natural logarithm
#' - `sinh(expr)` - Hyperbolic sine
#' - `cosh(expr)` - Hyperbolic cosine
#'
#' *Operations:*
#' - `+`, `-` - Addition, subtraction
#' - `*` - Multiplication (with first derivative product rule)
#' - `^` - Exponentiation (limited support)
#' - `()` - Parentheses for grouping
#'
#' *Numbers:*
#' - Any decimal number (use strings for high precision)
#' - Scientific notation supported
#'
#' **Chain Rule Examples:**
#' - `sin(2*x)` - Sine of scaled argument
#' - `cos(3*x)` - Cosine with coefficient
#' - `exp(-x)` - Negative exponential
#' - `log(2*x)` - Logarithm of scaled argument
#'
#' @section Accuracy Guarantee:
#' This function provides mathematically rigorous error bounds with accuracy
#' typically exceeding 1e-20 at 128-bit precision. All supported functions
#' include complete derivative implementations required by ARB's integration
#' algorithm.
#'
#' @section String Input Requirement:
#' For precision preservation, ALL inputs (expression, a, b) must be character
#' strings. Numeric inputs are automatically converted but may lose precision
#' for high-precision requests (>64 bits).
#'
#' @section Chain Rule Support:
#' First derivatives are supported for compositions f(g(x)) where g(x) is
#' a simple expression. Higher-order derivatives and complex compositions
#' are not yet implemented.
#'
#' @examples
#' \dontrun{
#' # Example 1: Basic trigonometric integral ∫₀^π sin(x) dx = 2
#' result1 <- integrate_rigorous("sin(x)", "0", "3.14159265358979323846")
#' print(result1)
#' cat("High precision result:", result1$value_str, "\n")
#'
#' # Example 2: Chain rule - sine with scaling ∫₀^π sin(2*x) dx = 1
#' result2 <- integrate_rigorous("sin(2*x)", "0", "3.14159265358979323846")
#' print(result2)
#'
#' # Example 3: Exponential decay ∫₀^1 exp(-x) dx = 1 - 1/e
#' result3 <- integrate_rigorous("exp(-x)", "0", "1")
#' print(result3)
#'
#' # Example 4: Polynomial ∫₀^1 x^2 dx = 1/3
#' result4 <- integrate_rigorous("x^2", "0", "1")
#' print(result4)
#'
#' # Example 5: Complex expression with addition
#' result5 <- integrate_rigorous("sin(x) + cos(x)", "0", "1.5707963267948966")
#' print(result5)
#'
#' # Example 6: High precision with string endpoints
#' result6 <- integrate_rigorous("cos(3*x)",
#'                              "0",
#'                              "1.047197551196597746154214461093167628066",
#'                              precision = 256)
#' cat("256-bit result:", result6$value_str, "\n")
#'
#' # Example 7: Hyperbolic functions ∫₀^1 sinh(x) dx = cosh(1) - 1
#' result7 <- integrate_rigorous("sinh(x)", "0", "1")
#' print(result7)
#'
#' # Example 8: Logarithmic integral ∫₁^e log(x) dx = 1
#' result8 <- integrate_rigorous("log(x)", "1", "2.718281828459045235360287471352662498")
#' print(result8)
#'
#' # Example 9: Product expression ∫₀^1 2*sin(x) dx
#' result9 <- integrate_rigorous("2*sin(x)", "0", "3.14159265358979323846")
#' print(result9)
#'
#' # Example 10: Verify against known result
#' arctangent_result <- integrate_rigorous("1/(1+x^2)", "0", "1")
#' pi_quarter <- pi/4
#' cat("∫₀¹ 1/(1+x²) dx =", arctangent_result$value, "\n")
#' cat("π/4 =", pi_quarter, "\n")
#' cat("Difference:", abs(arctangent_result$value - pi_quarter), "\n")
#' }
#'
#' @references
#' Johansson, F. (2017). Arb: Efficient Arbitrary-Precision Midpoint-Radius
#' Interval Arithmetic. IEEE Transactions on Computers.
#'
#' Johansson, F. (2018). Numerical integration in arbitrary-precision ball
#' arithmetic. arXiv:1802.07942.
#'
#' @seealso [list_supported_expressions()], [print.rigorous_result()], [compare_with_exact()]
#' @useDynLib integralARB, .registration=TRUE
#' @export
integrate_rigorous <- function(expression, a, b, precision = 128L,
                               abs_tol = NULL, rel_tol = NULL,
                               max_eval = 1000000L, max_depth = 100L,
                               verbose = FALSE) {

  # Input validation
  if (!is.character(expression) || length(expression) != 1) {
    stop("'expression' must be a single character string")
  }

  if (!is.character(a) || length(a) != 1) {
    if (is.numeric(a) && precision > 64) {
      warning("High precision (", precision, " bits) with numeric lower limit. ",
              "Consider using string: \"", a, "\" for full precision")
    }
    a <- as.character(a)
  }

  if (!is.character(b) || length(b) != 1) {
    if (is.numeric(b) && precision > 64) {
      warning("High precision (", precision, " bits) with numeric upper limit. ",
              "Consider using string: \"", b, "\" for full precision")
    }
    b <- as.character(b)
  }

  # Validate precision
  if (!is.numeric(precision) || precision < 64 || precision > 1024) {
    stop("'precision' must be between 64 and 1024 bits")
  }
  precision <- as.integer(precision)

  # Set default tolerances for high accuracy
  if (is.null(abs_tol)) abs_tol <- 2^(-precision-10)
  if (is.null(rel_tol)) rel_tol <- 2^(-precision-10)

  # Validate other parameters
  if (!is.logical(verbose)) {
    stop("'verbose' must be logical")
  }

  max_eval <- as.integer(max_eval)
  max_depth <- as.integer(max_depth)

  if (max_eval < 1000) {
    warning("max_eval < 1000 may lead to insufficient accuracy")
  }

  # Call C interface with expression parser
  result <- .Call("integrate_expression_c", expression, a, b, precision,
                  abs_tol, rel_tol, max_eval, max_depth, verbose,
                  PACKAGE = "integralARB")

  # Create S3 object with enhanced information
  class(result) <- "rigorous_result"
  result$call <- match.call()
  result$limits <- c(a, b)
  result$input_expression <- expression

  return(result)
}

#' List supported mathematical expressions
#'
#' Returns information about all supported mathematical expressions
#' and functions available for rigorous integration.
#'
#' @return Character string describing supported expressions
#' @export
#' @examples
#' \dontrun{
#' cat(list_supported_expressions())
#' }
list_supported_expressions <- function() {
  result <- .Call("list_builtin_functions_c", PACKAGE = "integralARB")
  return(result)
}

#' Print method for rigorous integration results
#'
#' @param x A rigorous_result object
#' @param digits Number of digits to display for double values
#' @param show_string Show high-precision string representation
#' @param show_call Show the original function call
#' @param ... Additional arguments
#' @export
print.rigorous_result <- function(x, digits = 15, show_string = TRUE,
                                  show_call = FALSE, ...) {
  cat("Rigorous Integration Result\n")
  cat("===========================\n")

  if (show_call && !is.null(x$call)) {
    cat("Call: ")
    print(x$call)
    cat("\n")
  }

  cat("Expression: ", x$expression %||% x$input_expression %||% "unknown", "\n")
  cat("Limits:     [", x$limits[1], ", ", x$limits[2], "]\n")
  cat("Precision:  ", x$precision, " bits\n\n")

  if (x$status == 0) {
    if (show_string && !is.null(x$value_str) && nchar(x$value_str) > 0) {
      cat("Value (high-precision):\n")
      cat("  ", x$value_str, "\n")
      if (!is.null(x$error_str) && nchar(x$error_str) > 0) {
        cat("Error bound:\n")
        cat("  ±", x$error_str, "\n")
      }
      cat("Double approximation: ", format(x$value, digits = digits), "\n")
    } else {
      cat("Value:      ", format(x$value, digits = digits), "\n")
      cat("Error:      ±", format(x$error_bound, digits = 5), "\n")
    }

    if (!is.null(x$rel_error) && is.finite(x$rel_error) && x$rel_error > 0) {
      cat("Rel. Error: ", format(x$rel_error, digits = 3), "\n")
    }

    # Show accuracy achieved
    if (!is.null(x$error_bound) && is.finite(x$error_bound) && x$error_bound > 0) {
      accuracy_digits <- -log10(x$error_bound)
      if (accuracy_digits > 0) {
        cat("Accuracy:   ~", floor(accuracy_digits), " decimal digits\n")
      }
    }

  } else {
    cat("Status:     FAILED (", x$status, ")\n")
    cat("Message:    ", x$message %||% "Integration failed", "\n")
  }

  cat("Evaluations:", x$evaluations %||% 0, "\n")

  invisible(x)
}

#' Get exact analytical result for known integrals
#'
#' Returns the exact analytical value for definite integrals when known.
#' Useful for validating numerical integration results.
#'
#' @param expression Mathematical expression string
#' @param a Lower limit (character string)
#' @param b Upper limit (character string)
#' @return Exact value if known, NA otherwise
#' @export
#' @examples
#' \dontrun{
#' # Known exact results
#' get_exact_result("sin(x)", "0", "3.14159265358979323846")  # Should be 2
#' get_exact_result("cos(x)", "0", "1.5707963267948966")      # Should be 1
#' get_exact_result("exp(x)", "0", "1")                       # Should be e-1
#' }
get_exact_result <- function(expression, a, b) {
  # Convert string limits to numeric for calculation
  a_num <- tryCatch(as.numeric(a), error = function(e) NA)
  b_num <- tryCatch(as.numeric(b), error = function(e) NA)

  if (is.na(a_num) || is.na(b_num)) {
    return(NA)
  }

  # Known exact results for common expressions
  if (expression == "sin(x)") {
    return(-cos(b_num) + cos(a_num))
  } else if (expression == "cos(x)") {
    return(sin(b_num) - sin(a_num))
  } else if (expression == "exp(x)") {
    return(exp(b_num) - exp(a_num))
  } else if (expression == "exp(-x)") {
    return(exp(-a_num) - exp(-b_num))
  } else if (expression == "x") {
    return((b_num^2 - a_num^2) / 2)
  } else if (expression == "x^2") {
    return((b_num^3 - a_num^3) / 3)
  } else if (expression == "sinh(x)") {
    return(cosh(b_num) - cosh(a_num))
  } else if (expression == "cosh(x)") {
    return(sinh(b_num) - sinh(a_num))
  } else if (expression == "sin(2*x)") {
    return((-cos(2*b_num) + cos(2*a_num)) / 2)
  } else if (expression == "cos(2*x)") {
    return((sin(2*b_num) - sin(2*a_num)) / 2)
  }

  # Special case: arctangent integral
  if (expression == "1/(1+x^2)" || expression == "1/(1+x*x)") {
    return(atan(b_num) - atan(a_num))
  }

  return(NA)
}

#' Compare numerical result with exact analytical result
#'
#' Compares the numerical integration result with the known exact value
#' when available, providing error analysis.
#'
#' @param result A rigorous_result object
#' @return List with comparison information
#' @export
#' @examples
#' \dontrun{
#' result <- integrate_rigorous("sin(x)", "0", "3.14159265358979323846")
#' comparison <- compare_with_exact(result)
#' print(comparison)
#' }
compare_with_exact <- function(result) {
  if (!inherits(result, "rigorous_result")) {
    stop("Input must be a rigorous_result object")
  }

  expression <- result$expression %||% result$input_expression
  if (is.null(expression)) {
    return(list(
      exact_available = FALSE,
      message = "No expression found in result object"
    ))
  }

  exact <- get_exact_result(expression, result$limits[1], result$limits[2])

  if (is.na(exact)) {
    return(list(
      exact_available = FALSE,
      message = paste("No exact result available for expression:", expression)
    ))
  }

  numerical_error <- abs(result$value - exact)
  relative_error <- numerical_error / abs(exact)

  list(
    exact_available = TRUE,
    expression = expression,
    exact_value = exact,
    numerical_value = result$value,
    absolute_error = numerical_error,
    relative_error = relative_error,
    error_bound = result$error_bound,
    bound_is_correct = numerical_error <= result$error_bound,
    accuracy_digits = -log10(relative_error),
    limits = result$limits
  )
}
