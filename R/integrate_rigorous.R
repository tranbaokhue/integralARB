#' Rigorous Numerical Integration with Built-in Analytical Functions
#'
#' Computes definite integrals using ball arithmetic with mathematically
#' rigorous error bounds. Only supports built-in analytical functions that
#' can provide complex holomorphic extensions for guaranteed accuracy.
#'
#' @param func_name Name of built-in function (see details for available functions)
#' @param a Lower limit of integration (numeric)
#' @param b Upper limit of integration (numeric)
#' @param precision Working precision in bits (default: 128, range: 64-1024)
#' @param abs_tol Absolute error tolerance (default: 2^(-precision-10))
#' @param rel_tol Relative error tolerance (default: 2^(-precision-10))
#' @param max_eval Maximum number of function evaluations (default: 1000000)
#' @param max_depth Maximum subdivision depth (default: 100)
#' @param verbose Print integration progress (default: FALSE)
#' @param params Additional parameters for parameterized functions (list)
#'
#' @return An object of class `rigorous_result` containing:
#' \itemize{
#'   \item `value` - Double approximation of integral value
#'   \item `value_str` - High-precision string representation
#'   \item `error_bound` - Rigorous error bound (double)
#'   \item `error_str` - High-precision error bound (string)
#'   \item `evaluations` - Number of function evaluations performed
#'   \item `status` - Integration status (0 = success)
#'   \item `precision` - Working precision used
#' }
#'
#' @details
#' **Available Built-in Functions:**
#'
#' *Polynomials:*
#' - `"constant"`: f(x) = 1
#' - `"linear"`: f(x) = x
#' - `"quadratic"`: f(x) = x²
#' - `"cubic"`: f(x) = x³
#'
#' *Rational Functions:*
#' - `"reciprocal"`: f(x) = 1/x
#' - `"arctangent"`: f(x) = 1/(1+x²)
#'
#' *Exponential Functions:*
#' - `"exp"`: f(x) = exp(x)
#' - `"exp_minus"`: f(x) = exp(-x)
#' - `"gaussian"`: f(x) = exp(-x²)
#'
#' *Logarithmic Functions:*
#' - `"log"`: f(x) = ln(x)
#'
#' *Trigonometric Functions:*
#' - `"sin"`: f(x) = sin(x)
#' - `"cos"`: f(x) = cos(x)
#' - `"sinc"`: f(x) = sin(x)/x
#'
#' *Hyperbolic Functions:*
#' - `"sinh"`: f(x) = sinh(x)
#' - `"cosh"`: f(x) = cosh(x)
#' - `"sech"`: f(x) = sech(x)
#'
#' *Special Integrands:*
#' - `"gamma_0"`: f(x) = exp(-x) [for Γ(1)]
#' - `"gamma_1"`: f(x) = x*exp(-x) [for Γ(2)]
#' - `"gamma_2"`: f(x) = x²*exp(-x) [for Γ(3)]
#'
#' @section Accuracy Guarantee:
#' This function provides mathematically rigorous error bounds with accuracy
#' typically exceeding 1e-20 at 128-bit precision. All built-in functions
#' include proper complex holomorphic extensions required by the ARB library.
#'
#' @section Composition Strategy:
#' To integrate complex expressions, decompose them into sums/differences of
#' built-in functions and integrate each term separately:
#'
#' For ∫(x² + sin(x))dx from 0 to π:
#' ```
#' term1 <- integrate_rigorous("quadratic", 0, pi)
#' term2 <- integrate_rigorous("sin", 0, pi)
#' total <- term1$value + term2$value  # π³/3 + 2
#' ```
#'
#' @examples
#' \dontrun{
#' # Example 1: Arctangent integral ∫₀¹ 1/(1+x²) dx = π/4
#' result1 <- integrate_rigorous("arctangent", 0, 1)
#' print(result1)
#' # Expected: 0.785398163397448309615660845819875721049...
#'
#' # Example 2: Gaussian integral ∫₀² exp(-x²) dx
#' result2 <- integrate_rigorous("gaussian", 0, 2)
#' print(result2)
#' # Expected: 0.882081390762421679967481035914054037224...
#'
#' # Example 3: Sine integral ∫₀π sin(x) dx = 2
#' result3 <- integrate_rigorous("sin", 0, pi)
#' print(result3)
#' # Expected: exactly 2.0
#'
#' # Example 4: High precision integration
#' result4 <- integrate_rigorous("gaussian", 0, 1, precision = 256)
#' cat("256-bit result:", result4$value_str, "\n")
#'
#' # Example 5: Gamma function integral ∫₀∞ x²exp(-x) dx = Γ(3) = 2
#' result5 <- integrate_rigorous("gamma_2", 0, 10)  # Truncate at 10
#' print(result5)
#' # Expected: ≈ 2.0 (with small truncation error)
#' }
#'
#' @references
#' Johansson, F. (2017). Arb: Efficient Arbitrary-Precision Midpoint-Radius
#' Interval Arithmetic. IEEE Transactions on Computers.
#'
#' @seealso [list_builtin_functions()], [print.rigorous_result()]
#' @useDynLib integralARB, .registration=TRUE
#' @export
integrate_rigorous <- function(func_name, a, b, precision = 128L,
                               abs_tol = NULL, rel_tol = NULL,
                               max_eval = 1000000L, max_depth = 100L,
                               verbose = FALSE, params = NULL) {

  # Input validation
  stopifnot(
    is.character(func_name), length(func_name) == 1,
    is.numeric(a), is.numeric(b), length(a) == 1, length(b) == 1,
    is.numeric(precision), precision >= 64, precision <= 1024,
    is.logical(verbose)
  )

  # Set default tolerances for high accuracy (1e-20 or better)
  if (is.null(abs_tol)) abs_tol <- 2^(-precision-10)
  if (is.null(rel_tol)) rel_tol <- 2^(-precision-10)

  # Validate function name
  valid_functions <- c(
    "constant", "linear", "quadratic", "cubic",
    "reciprocal", "arctangent",
    "exp", "exp_minus", "gaussian",
    "log",
    "sin", "cos", "sinc",
    "sinh", "cosh", "sech",
    "gamma_0", "gamma_1", "gamma_2"
  )

  if (!func_name %in% valid_functions) {
    stop("Unknown function '", func_name, "'. Available functions:\n  ",
         paste(valid_functions, collapse = ", "))
  }

  # Call C interface
  result <- .Call("integrate_rigorous_builtin_c", func_name, a, b, precision,
                  abs_tol, rel_tol, max_eval, max_depth, verbose, params,
                  PACKAGE = "integralARB")

  # Create S3 object
  class(result) <- "rigorous_result"
  result$call <- match.call()
  result$precision <- precision
  result$function_name <- func_name
  result$limits <- c(a, b)

  return(result)
}

#' List all available built-in functions
#'
#' Returns a data frame with information about all built-in functions
#' available for rigorous integration.
#'
#' @return Data frame with columns: name, description, latex_form, exact_antiderivative
#' @export
list_builtin_functions <- function() {
  data.frame(
    name = c(
      "constant", "linear", "quadratic", "cubic",
      "reciprocal", "arctangent",
      "exp", "exp_minus", "gaussian",
      "log",
      "sin", "cos", "sinc",
      "sinh", "cosh", "sech",
      "gamma_0", "gamma_1", "gamma_2"
    ),
    description = c(
      "f(x) = 1", "f(x) = x", "f(x) = x²", "f(x) = x³",
      "f(x) = 1/x", "f(x) = 1/(1+x²)",
      "f(x) = exp(x)", "f(x) = exp(-x)", "f(x) = exp(-x²)",
      "f(x) = ln(x)",
      "f(x) = sin(x)", "f(x) = cos(x)", "f(x) = sin(x)/x",
      "f(x) = sinh(x)", "f(x) = cosh(x)", "f(x) = sech(x)",
      "f(x) = exp(-x)", "f(x) = x*exp(-x)", "f(x) = x²*exp(-x)"
    ),
    exact_antiderivative = c(
      "x", "x²/2", "x³/3", "x⁴/4",
      "ln|x|", "arctan(x)",
      "exp(x)", "-exp(-x)", "√π/2 * erf(x)",
      "x*ln(x) - x",
      "-cos(x)", "sin(x)", "Si(x)",
      "cosh(x)", "sinh(x)", "arctan(sinh(x))",
      "Γ(1) = 1", "Γ(2) = 1", "Γ(3) = 2"
    ),
    stringsAsFactors = FALSE
  )
}

#' Print method for rigorous integration results
#' @param x A rigorous_result object
#' @param digits Number of digits to display for double values
#' @param show_string Show high-precision string representation
#' @param ... Additional arguments
#' @export
print.rigorous_result <- function(x, digits = 15, show_string = TRUE, ...) {
  cat("Rigorous Integration Result\n")
  cat("===========================\n")
  cat("Function: ", x$function_name %||% "unknown", "\n")
  cat("Limits:   [", x$limits[1], ", ", x$limits[2], "]\n")
  cat("Precision:", x$precision, "bits\n\n")

  if (x$status == 0) {
    if (show_string && !is.null(x$value_str) && nchar(x$value_str) > 0) {
      cat("Value (high-precision):\n")
      cat("  ", x$value_str, "\n")
      cat("Error bound:\n")
      cat("  ±", x$error_str, "\n")
      cat("Double approximation:", format(x$value, digits = digits), "\n")
    } else {
      cat("Value:     ", format(x$value, digits = digits), "\n")
      cat("Error:     ±", format(x$error_bound, digits = 5), "\n")
    }

    if (x$rel_error < Inf && x$rel_error > 0) {
      cat("Rel. Err:  ", format(x$rel_error, digits = 3), "\n")
    }

    # Show accuracy achieved
    if (!is.null(x$error_bound) && x$error_bound > 0) {
      accuracy_digits <- -log10(x$error_bound)
      if (accuracy_digits > 0) {
        cat("Accuracy:  ~", floor(accuracy_digits), "decimal digits\n")
      }
    }

  } else {
    cat("Status:    FAILED (", x$status, ")\n")
    cat("Message:   ", x$message %||% "Integration failed", "\n")
  }

  cat("Evaluations:", x$evaluations %||% 0, "\n")

  invisible(x)
}

#' Get exact analytical result for known integrals
#'
#' Returns the exact analytical value for definite integrals of built-in
#' functions over common intervals, when known.
#'
#' @param func_name Function name
#' @param a Lower limit
#' @param b Upper limit
#' @return Exact value if known, NA otherwise
#' @export
get_exact_result <- function(func_name, a, b) {

  # Known exact results for common intervals
  exact_results <- list(
    # π/4 integrals
    "arctangent" = function(a, b) {
      if (a == 0 && b == 1) return(pi/4)
      if (a == 0 && b == Inf) return(pi/2)
      return(atan(b) - atan(a))
    },

    # Polynomial integrals
    "constant" = function(a, b) b - a,
    "linear" = function(a, b) (b^2 - a^2)/2,
    "quadratic" = function(a, b) (b^3 - a^3)/3,
    "cubic" = function(a, b) (b^4 - a^4)/4,

    # Trigonometric integrals
    "sin" = function(a, b) {
      if (a == 0 && b == pi) return(2)
      return(cos(a) - cos(b))
    },
    "cos" = function(a, b) sin(b) - sin(a),

    # Exponential integrals
    "exp" = function(a, b) exp(b) - exp(a),
    "exp_minus" = function(a, b) exp(-a) - exp(-b),

    # Gaussian integral ∫₀^∞ exp(-x²) dx = √π/2
    "gaussian" = function(a, b) {
      if (a == 0 && b == Inf) return(sqrt(pi)/2)
      # For finite limits, no simple closed form
      return(NA)
    },

    # Hyperbolic functions
    "sinh" = function(a, b) cosh(b) - cosh(a),
    "cosh" = function(a, b) sinh(b) - sinh(a),

    # Gamma integrals (over [0,∞])
    "gamma_0" = function(a, b) {
      if (a == 0 && b == Inf) return(gamma(1))  # = 1
      return(NA)
    },
    "gamma_1" = function(a, b) {
      if (a == 0 && b == Inf) return(gamma(2))  # = 1
      return(NA)
    },
    "gamma_2" = function(a, b) {
      if (a == 0 && b == Inf) return(gamma(3))  # = 2
      return(NA)
    }
  )

  if (func_name %in% names(exact_results)) {
    tryCatch({
      exact_results[[func_name]](a, b)
    }, error = function(e) NA)
  } else {
    NA
  }
}

#' Compare numerical result with exact analytical result
#'
#' @param result A rigorous_result object
#' @return List with comparison information
#' @export
compare_with_exact <- function(result) {
  exact <- get_exact_result(result$function_name, result$limits[1], result$limits[2])

  if (is.na(exact)) {
    return(list(
      exact_available = FALSE,
      message = "No exact result available for this function and limits"
    ))
  }

  numerical_error <- abs(result$value - exact)
  relative_error <- numerical_error / abs(exact)

  list(
    exact_available = TRUE,
    exact_value = exact,
    numerical_value = result$value,
    absolute_error = numerical_error,
    relative_error = relative_error,
    error_bound = result$error_bound,
    bound_is_correct = numerical_error <= result$error_bound,
    accuracy_digits = -log10(relative_error)
  )
}
