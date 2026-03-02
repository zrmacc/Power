# Purpose: Input validation helpers for power and sample-size calculations.

#' Validate design matrices and variance for linear regression power
#'
#' Checks that \code{beta_x}, \code{var_x}, \code{var_z}, and \code{cov_xz}
#' have consistent dimensions and that \code{var_resid} is positive. Stops with
#' an error if any check fails; otherwise returns invisibly.
#'
#' @param beta_x Regression coefficient for X (matrix, e.g. column vector).
#' @param cov_xz Covariance of X and Z (matrix).
#' @param var_resid Residual variance of Y|(X, Z), scalar.
#' @param var_x Variance of X (matrix).
#' @param var_z Variance of Z (matrix).
#' @return \code{invisible(NULL)}. Stops with an error if validation fails.
#' @noRd
.ValidateLinRegDesign <- function(beta_x, cov_xz, var_resid, var_x, var_z) {
  if (any(var_resid <= 0)) {
    stop("var_resid must be positive.")
  }
  if (nrow(var_x) != nrow(beta_x) || ncol(cov_xz) != ncol(var_z) ||
      nrow(cov_xz) != nrow(var_x)) {
    stop("Dimension mismatch among beta_x, var_x, var_z, and cov_xz.")
  }
  return(invisible(NULL))
}

#' Validate sample size is positive
#'
#' Stops with an error if any element of \code{n} is non-positive. Used for
#' power calculations that require a valid sample size.
#'
#' @param n Sample size (numeric, possibly vector).
#' @return \code{invisible(NULL)}. Stops with an error if \code{n <= 0}.
#' @noRd
.ValidatePositiveN <- function(n) {
  if (any(n <= 0)) {
    stop("n must be positive.")
  }
  return(invisible(NULL))
}

#' Validate sample size search bounds
#'
#' Ensures \code{max_n >= 1} and \code{min_n} is in \code{[1, max_n]}.
#' Stops with an error if the bounds are invalid.
#'
#' @param min_n Lower bound for sample size search.
#' @param max_n Upper bound for sample size search.
#' @return \code{invisible(NULL)}. Stops with an error if bounds are invalid.
#' @noRd
.ValidateSampleSizeBounds <- function(min_n, max_n) {
  if (max_n < 1) {
    stop("max_n must be at least 1.")
  }
  if (min_n < 1 || min_n > max_n) {
    stop("min_n must be at least 1 and not greater than max_n.")
  }
  return(invisible(NULL))
}

#' Validate target power is in (0, 1)
#'
#' Stops with an error if \code{power} is not a single number in the open
#' interval (0, 1).
#'
#' @param power Target power (numeric, length 1).
#' @return \code{invisible(NULL)}. Stops with an error if invalid.
#' @noRd
.ValidateTargetPower <- function(power) {
  if (power <= 0 || power >= 1) {
    stop("power must be in (0, 1).")
  }
  return(invisible(NULL))
}

#' Validate that beta_x is not the null (zero) effect
#'
#' For sample size calculation, the target power is undefined when the
#' effect is null. Stops with an error if all elements of \code{beta_x}
#' are zero.
#'
#' @param beta_x Regression coefficient for X (matrix).
#' @return \code{invisible(NULL)}. Stops with an error if \code{all(beta_x == 0)}.
#' @noRd
.ValidateNonZeroBeta <- function(beta_x) {
  if (all(beta_x == 0)) {
    stop("Target power is undefined when beta_x is zero (null effect).")
  }
  return(invisible(NULL))
}

#' Validate design for logistic regression power
#'
#' Checks that \code{beta_0} is scalar, \code{beta_x}, \code{beta_z},
#' \code{var_x}, \code{var_z}, and \code{cov_xz} have consistent dimensions,
#' and that the combined covariance of \eqn{(X, Z)} is positive definite.
#'
#' @param beta_0 Intercept (scalar).
#' @param beta_x Coefficient for X (matrix).
#' @param beta_z Coefficient for Z (matrix).
#' @param var_x Variance of X (matrix).
#' @param var_z Variance of Z (matrix).
#' @param cov_xz Covariance of X and Z (matrix).
#' @return \code{invisible(NULL)}. Stops with an error if validation fails.
#' @noRd
.ValidateLogisticDesign <- function(
    beta_0,
    beta_x,
    beta_z,
    var_x,
    var_z,
    cov_xz
) {
  if (length(beta_0) != 1L) {
    stop("beta_0 must be a scalar.")
  }
  p <- nrow(as.matrix(beta_x))
  q <- nrow(as.matrix(beta_z))
  var_x <- as.matrix(var_x)
  var_z <- as.matrix(var_z)
  cov_xz <- as.matrix(cov_xz)
  if (nrow(var_x) != p || ncol(var_x) != p) {
    stop("var_x must be a ", p, " x ", p, " matrix.")
  }
  if (nrow(var_z) != q || ncol(var_z) != q) {
    stop("var_z must be a ", q, " x ", q, " matrix.")
  }
  if (nrow(cov_xz) != p || ncol(cov_xz) != q) {
    stop("cov_xz must be a ", p, " x ", q, " matrix.")
  }
  Sigma <- rbind(
    cbind(var_x, cov_xz),
    cbind(t(cov_xz), var_z)
  )
  tryCatch(
    chol(Sigma),
    error = function(e) {
      stop(
        "The combined covariance of (X, Z) must be positive definite. ",
        conditionMessage(e)
      )
    }
  )
  return(invisible(NULL))
}
