#' Rigorous Numerical Integration with Guaranteed Error Bounds
#'
#' Computes definite integrals using ball arithmetic with mathematically
#' rigorous error bounds. Unlike heuristic methods, this guarantees that
#' the true integral value lies within the returned interval.
#'
#' @param f An R function, or the string "arctangent" for built-in test function
#' @param a Lower limit of integration (numeric)
#' @param b Upper limit of integration (numeric)
#' @param precision Working precision in bits (default: 64)
#' @param abs_tol Absolute error tolerance (default: 2^(-precision))
#' @param rel_tol Relative error tolerance (default: 2^(-precision))
#' @param max_eval Maximum number of function evaluations (default: 100000)
#' @param max_depth Maximum subdivision depth (default: 50)
#' @param verbose Print integration progress (default: FALSE)
#'
#' @return An object of class `rigorous_result` containing:
#' \itemize{
#'   \item `value` - The computed integral value
#'   \item `error_bound` - Rigorous error bound
#'   \item `abs_error` - Absolute error estimate
#'   \item `rel_error` - Relative error estimate
#'   \item `evaluations` - Number of function evaluations performed
#'   \item `subdivisions` - Number of interval subdivisions
#'   \item `status` - Integration status (0 = success)
#'   \item `message` - Status message
#' }
#'
#' @details
#' This function uses adaptive Gauss-Legendre quadrature with rigorous
#' error bounds based on complex analysis. It automatically detects
#' discontinuities and singularities, ensuring reliable results even
#' for pathological functions that break other integration methods.
#'
#' The integration is performed using the FLINT/ARB library, which provides
#' mathematically guaranteed error bounds rather than heuristic estimates.
#'
#' @examples
#' \dontrun{
#' # Example 1: Built-in arctangent function ∫[0,1] 1/(1+x^2) dx = π/4
#' result1 <- integrate_rigorous("arctangent", 0, 1)
#' print(result1)
#' # Expected: value = 0.785398163397448309615660845819875721048822
#'
#' # Example 2: Sine function ∫[0,π] sin(x) dx = 2
#' result2 <- integrate_rigorous(sin, 0, pi)
#' print(result2)
#' # Expected: value ≈ 2.0
#'
#' # Example 3: Gaussian function ∫[0,2] exp(-x²) dx
#' gaussian_func <- function(x) exp(-x^2)
#' result3 <- integrate_rigorous(gaussian_func, 0, 2)
#' print(result3)
#' # Expected: value ≈ 0.88208139076242167996748103591405403722405...
#' }
#'
#' @references
#' Johansson, F. (2017). Arb: Efficient Arbitrary-Precision Midpoint-Radius
#' Interval Arithmetic. IEEE Transactions on Computers.
#'
#' @seealso [integrate()], [print.rigorous_result()]
#' @export
integrate_rigorous <- function(f, a, b, precision = 128L,
                               abs_tol = NULL, rel_tol = NULL,
                               max_eval = 100000L, max_depth = 50L,
                               verbose = FALSE) {

  # Input validation
  stopifnot(
    is.numeric(a), is.numeric(b), length(a) == 1, length(b) == 1,
    is.numeric(precision), precision >= 32, precision <= 4096,
    is.logical(verbose)
  )

  # Set default tolerances
  if (is.null(abs_tol)) abs_tol <- 2^(-precision-10)
  if (is.null(rel_tol)) rel_tol <- 2^(-precision-10)

  # Call C interface
  result <- .Call("integrate_rigorous_c", f, a, b, precision,
                  abs_tol, rel_tol, max_eval, max_depth, verbose,
                  PACKAGE = "integralARB")

  # Create S3 object
  class(result) <- "rigorous_result"
  result$call <- match.call()
  result$precision <- precision
  result$expression <- f
  result$limits <- c(a, b)

  return(result)
}

#' @export
print.rigorous_result <- function(x, digits = 10, ...) {
  cat("Rigorous Integration Result\n")
  cat("===========================\n")

  if (is.function(x$expression)) {
    cat("Function:  [R function]\n")
  } else {
    cat("Function: ", deparse(x$expression), "\n")
  }

  cat("Limits:   [", x$limits[1], ", ", x$limits[2], "]\n")
  cat("Precision:", x$precision, "bits\n\n")

  if (x$status == 0) {
    # Use high-precision string if available, otherwise fall back to double
    if (!is.null(x$value_str) && nchar(x$value_str) > 0) {
      cat("Value:     ", x$value_str, "\n")
      cat("Error:     ±", x$error_str, "\n")
    } else {
      cat("Value:     ", format(x$value, digits = digits), "\n")
      cat("Error:     ±", format(x$error_bound, digits = 3), "\n")
    }

    if (x$rel_error < Inf) {
      cat("Rel. Err:  ", format(x$rel_error, digits = 3), "\n")
    }
  } else {
    cat("Status:    FAILED (", x$status, ")\n")
    cat("Message:   ", x$message, "\n")
  }

  cat("Evals:     ", x$evaluations, "\n")

  invisible(x)
}

#' Test if FLINT/ARB is working correctly
#'
#' @return The value of π computed using ARB
#' @export
test_flint <- function() {
  .Call("test_flint_basic", PACKAGE = "integralARB")
}
