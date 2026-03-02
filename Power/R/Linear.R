# Purpose: Power calculation for linear regression.

# -----------------------------------------------------------------------------
# Power computation
# -----------------------------------------------------------------------------

#' Power from sample size and precomputed efficient information
#'
#' Internal helper: computes Wald-test power for \eqn{H_0: \beta_X = 0} given
#' sample size \code{n} and the efficient information matrix for \eqn{\beta_X}
#' (Schur complement). Used by \code{\link{SampleSizeLinReg}} to avoid
#' recomputing the information in each root-finding iteration.
#'
#' @param n Sample size (may be a vector).
#' @param beta_x True regression coefficient for X, matrix (e.g. column vector).
#' @param var_resid Residual variance of Y|(X, Z), scalar.
#' @param ixz Efficient information for \eqn{\beta_X}, matrix (output of
#'   \code{SchurC(var_x, var_z, cov_xz)}).
#' @param t1e Type I error level.
#' @return Numeric power (same length as \code{n}).
#' @noRd
.PowerFromNCP <- function(n, beta_x, var_resid, ixz, t1e) {
  ncp <- n * (1 / var_resid) * as.numeric(matQF(beta_x, ixz))
  k <- nrow(ixz)
  crit <- stats::qchisq(p = 1 - t1e, df = k)
  out <- stats::pchisq(q = crit, df = k, ncp = ncp, lower.tail = FALSE)
  return(out)
}

#' Power Calculation
#'
#' Power for testing \eqn{H_{0}:\beta_{X} = 0} in the linear model
#' \eqn{Y \sim X + Z}. All variables (X, Y, Z) are assumed centered.
#'
#' @param beta_x True regression coefficient for X, vector.
#' @param cov_xz Covariance of X and Z, matrix.
#' @param n Sample size.
#' @param t1e Type 1 error level.
#' @param var_resid Residual variance of Y|(X, Z), scalar.
#' @param var_x Variance of X, matrix.
#' @param var_z Variance of Z, matrix.
#' @export
#' @return Numeric power.
#'
#' @details Power is computed using the non-central chi-squared distribution
#'   under the Wald test. For multivariate X, the test has \code{nrow(var_x)}
#'   degrees of freedom.
#'
#' @seealso \code{\link{SampleSizeLinReg}} for sample size given a target power.
PowerLinReg <- function(
  beta_x,
  cov_xz,
  n,
  t1e = 0.05,
  var_resid,
  var_x,
  var_z
) {
  beta_x <- as.matrix(beta_x)
  cov_xz <- as.matrix(cov_xz)
  var_x <- as.matrix(var_x)
  var_z <- as.matrix(var_z)

  .ValidatePositiveN(n)
  .ValidateLinRegDesign(beta_x, cov_xz, var_resid, var_x, var_z)

  # Null effect: power equals type I error.
  if (all(beta_x == 0)) {
    return(t1e)
  }

  # Non-centrality parameter.
  ixz <- SchurC(var_x, var_z, cov_xz)
  ncp <- n * (1 / var_resid) * as.numeric(matQF(beta_x, ixz))

  # Critical value and power.
  k <- nrow(var_x)
  crit <- stats::qchisq(p = 1 - t1e, df = k)
  power <- stats::pchisq(q = crit, df = k, ncp = ncp, lower.tail = FALSE)
  return(power)
}


# -----------------------------------------------------------------------------

#' Sample Size Calculation
#'
#' Sample size for testing \eqn{H_{0}:\beta_{X} = 0} in the linear model
#' \eqn{Y \sim X + Z}. All variables (X, Y, Z) are assumed centered.
#'
#' @param beta_x True regression coefficient for X, vector.
#' @param cov_xz Covariance of X and Z, matrix.
#' @param max_n Maximum allowable sample size (upper bound for search).
#' @param power Target power.
#' @param t1e Type 1 error level.
#' @param var_resid Residual variance of Y|(X, Z), scalar.
#' @param var_x Variance of X, matrix.
#' @param var_z Variance of Z, matrix.
#' @param min_n Minimum allowable sample size (lower bound for search).
#'   Default 1.
#' @param tol Numeric tolerance passed to \code{\link[stats]{uniroot}} for
#'   the root-finding search. Default matches \code{uniroot}'s default
#'   (\code{.Machine$double.eps^0.25}).
#' @export
#' @return Data.frame with columns \code{n} (sample size) and \code{power}
#'   (achieved power at that n).
#'
#' @details Uses \code{\link[stats]{uniroot}} to find the smallest integer
#'   sample size in \code{[min_n, max_n]} at which power reaches the target.
#'   Stops with an error if target power is not attainable at \code{max_n}.
#'
#' @seealso \code{\link{PowerLinReg}} for power at a given sample size.
SampleSizeLinReg <- function(
  beta_x,
  cov_xz,
  max_n,
  power = 0.80,
  t1e = 0.05,
  var_resid,
  var_x,
  var_z,
  min_n = 1,
  tol = .Machine$double.eps^0.25
) {
  beta_x <- as.matrix(beta_x)
  cov_xz <- as.matrix(cov_xz)
  var_x <- as.matrix(var_x)
  var_z <- as.matrix(var_z)

  .ValidateSampleSizeBounds(min_n, max_n)
  .ValidateTargetPower(power)
  .ValidateLinRegDesign(beta_x, cov_xz, var_resid, var_x, var_z)
  .ValidateNonZeroBeta(beta_x)

  # Efficient information (computed once).
  ixz <- SchurC(var_x, var_z, cov_xz)

  # Objective: power(n) - target.
  obj <- function(n) {
    delta <- .PowerFromNCP(n, beta_x, var_resid, ixz, t1e) - power
    return(delta)
  }

  min_power <- .PowerFromNCP(min_n, beta_x, var_resid, ixz, t1e)
  max_power <- .PowerFromNCP(max_n, beta_x, var_resid, ixz, t1e)
  if (max_power < power) {
    stop(
      "Target power ", round(power, 3), " is unattainable with max_n = ", max_n,
      ". Achieved power at max_n is ", round(max_power, 3), ". ",
      "Increase max_n and try again."
    )
  }
  if (min_power >= power) {
    final_n <- min_n
  } else {
    zero <- stats::uniroot(f = obj, lower = min_n, upper = max_n, tol = tol)
    final_n <- ceiling(zero$root)
  }
  final_power <- .PowerFromNCP(final_n, beta_x, var_resid, ixz, t1e)

  out <- data.frame(
    n = final_n,
    power = final_power
  )
  return(out)
}
