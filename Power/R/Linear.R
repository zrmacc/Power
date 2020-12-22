# Purpose: Power calculation for linear regression.

#' Power Calculation
#' 
#' Power for testing \eqn{H_{0}:\beta_{X} = 0} in the linear model
#' \eqn{Y ~ X + Z}. All variables (X, Y, Z) are assumed centered. 
#' 
#' @param beta_x True regression coefficient for X, vector.
#' @param cov_xz Covariance of of X and Z, matrix.
#' @param n Sample size.
#' @param t1e Type 1 error level.
#' @param var_resid Residual variance of Y|(X, Z), scalar.
#' @param var_x Variance of X, matrix.
#' @param var_z Variance of Z, matrix.
#' @importFrom stats pchisq qchisq
#' @export 
#' @return Numeric power.

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
  
  # Non-centrality parameter.
  ixz <- SchurC(var_x, var_z, cov_xz)
  ncp <- n * (1 / var_resid) * as.numeric(matQF(beta_x, ixz))
  
  # Critical value.
  k <- nrow(var_x)
  
  # Power.
  crit <- qchisq(p = 1 - t1e, df = k)
  power <- pchisq(q = crit, df = k, ncp = ncp, lower.tail = FALSE)
  return(power)
}


# -----------------------------------------------------------------------------

#' Sample Size Calculation
#' 
#' Sample size for testing \eqn{H_{0}:\beta_{X} = 0} in the linear model
#' \eqn{Y ~ X + Z}. All variables (X, Y, Z) are assumed centered. 
#' 
#' @param beta_x True regression coefficient for X, vector.
#' @param cov_xz Covariance of of X and Z, matrix.
#' @param max_n Maximum allowable sample size.
#' @param power Target power. 
#' @param t1e Type 1 error level.
#' @param var_resid Residual variance of Y|(X, Z), scalar.
#' @param var_x Variance of X, matrix.
#' @param var_z Variance of Z, matrix.
#' @importFrom stats uniroot
#' @export 
#' @return Data.frame with final sampe size and power.

SampleSizeLinReg <- function(
  beta_x,
  cov_xz,
  max_n,
  power = 0.80,
  t1e = 0.05,
  var_resid,
  var_x,
  var_z
) {
  
  # Power wrapper.
  WrapPower <- function(n) {
    out <- PowerLinReg(
      beta_x = beta_x,
      cov_xz = cov_xz,
      n = n, 
      t1e = t1e,
      var_resid,
      var_x,
      var_z
    )
    return(out)
  }
  
  # Maximum possible power.
  max_power <- WrapPower(max_n)
  if (max_power < power) {
    stop("Target power is unattainable with current max_n.\n")
  }
  
  # Objective function.
  Obj <- function(n) {return(WrapPower(n) - power)}
  
  # Find sample size.
  zero <- uniroot(f = Obj, lower = 1, upper = max_n)
  
  # Sample size.
  final_n <- ceiling(zero$root)
  
  # Power at final sample size.
  final_power <- WrapPower(final_n)
  
  # Output.
  out <- data.frame(
    n = final_n,
    power = final_power
  )
  return(out)
}