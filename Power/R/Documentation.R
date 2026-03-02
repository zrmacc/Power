# Purpose: Package documentation
# Updated: 2020-12-22

#' @useDynLib Power, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Power
#'
#' Power and sample-size calculations for regression models. The package
#' assumes all outcome and predictor variables are centered (mean zero).
#'
#' **Linear regression:** For the model \eqn{Y \sim X + Z}, the functions
#' \code{\link{PowerLinReg}} and \code{\link{SampleSizeLinReg}} compute power
#' for the Wald test of \eqn{H_0: \beta_X = 0} and the sample size required
#' to achieve a target power, respectively. Both support univariate or
#' multivariate X and Z via matrix arguments (\code{var_x}, \code{var_z},
#' \code{cov_xz}).
#'
#' **Logistic regression:** \code{\link{PowerLogisticReg}} and
#' \code{\link{SampleSizeLogisticReg}} provide power and sample size for the
#' Wald test of \eqn{H_0: \beta_X = 0} in the logistic model.
#'
#' @author Zachary R. McCaw
#' @keywords internal
"_PACKAGE"
