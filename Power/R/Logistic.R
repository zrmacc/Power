# Purpose: Power and sample-size calculation for logistic regression.
# Updated: 2026-03-02

# -----------------------------------------------------------------------------
# Expected information matrix (internal)
# -----------------------------------------------------------------------------

#' Gauss-Hermite quadrature rule for standard normal
#'
#' Returns nodes and weights such that \eqn{E[g(Z)] \approx \sum_k w_k g(z_k)}
#' for \eqn{Z \sim N(0, 1)}. Uses \code{statmod::gauss.quad} (weight
#' \eqn{\exp(-x^2)}) and converts to the N(0,1) rule via \eqn{z = \sqrt{2}\,x}
#' and \eqn{w_{\mathrm{norm}} = w_{\mathrm{hermite}} / \sqrt{\pi}}.
#'
#' @param n Number of quadrature points.
#' @return List with \code{nodes} (length \code{n}) and \code{weights}
#'   (length \code{n}).
#' @noRd
.GaussHermiteNormal <- function(n) {
  # Standard Hermite rule is for weight exp(-x^2). Convert to N(0,1): z = sqrt(2)*x,
  # and weight for N(0,1) integral is w_hermite / sqrt(pi).
  rule <- statmod::gauss.quad(n, kind = "hermite")
  nodes <- sqrt(2) * rule$nodes
  weights <- rule$weights / sqrt(pi)
  out <- list(nodes = nodes, weights = weights)
  return(out)
}


#' Expected information matrix for logistic regression (Gauss-Hermite)
#'
#' Approximates \eqn{i_{\beta\beta'} = E[W W' \pi(1-\pi)]} using product
#' Gauss-Hermite quadrature over the distribution of \eqn{(X, Z) \sim N(0,
#' \Sigma)}. The number of quadrature points is \code{n_quad^d} where
#' \eqn{d = \dim(X) + \dim(Z)}. Accurate and fast for small \eqn{d} (e.g. 1–4).
#'
#' @param beta_0 Intercept (scalar).
#' @param beta_x Coefficient for X (matrix).
#' @param beta_z Coefficient for Z (matrix).
#' @param var_x Variance of X (matrix).
#' @param var_z Variance of Z (matrix).
#' @param cov_xz Covariance of X and Z (matrix).
#' @param n_quad Number of quadrature points per dimension.
#' @return The expected information matrix (numeric matrix).
#' @noRd
.LogisticInformationMatrixQuadrature <- function(
    beta_0,
    beta_x,
    beta_z,
    var_x,
    var_z,
    cov_xz,
    n_quad = 15L
) {
  # Dimensions: p = dim(X), q = dim(Z), d = dim(X, Z) for the quadrature.
  beta_x <- as.matrix(beta_x)
  beta_z <- as.matrix(beta_z)
  p <- nrow(beta_x)
  q <- nrow(beta_z)
  d <- p + q
  var_x <- as.matrix(var_x)
  var_z <- as.matrix(var_z)
  cov_xz <- as.matrix(cov_xz)

  # One-dimensional Gauss-Hermite rule for N(0, 1): nodes and weights so that
  # E[g(Z)] \approx sum(weights * g(nodes)) for Z ~ N(0, 1).
  rule_1d <- .GaussHermiteNormal(n_quad)
  nodes_1d <- rule_1d$nodes
  weights_1d <- rule_1d$weights

  # Covariance of (X, Z) and its Cholesky factor so that (X, Z)' = l_chol z_std
  # with z_std ~ N(0, I_d). We integrate over z_std and transform.
  sigma <- rbind(
    cbind(var_x, cov_xz),
    cbind(t(cov_xz), var_z)
  )
  l_chol <- t(chol(sigma))

  # Accumulator for the expected information: (1 + p + q) x (1 + p + q) with
  # rows/cols ordered as (intercept, X, Z).
  n_par <- 1L + p + q
  i_sum <- matrix(0, n_par, n_par)

  # Tensor-product grid: each point is a combination of d 1D quadrature
  # indices. grid[i, j] = which 1D node to use for dimension j at point i.
  grid <- as.matrix(expand.grid(rep(list(seq_len(n_quad)), d)))
  n_pts <- nrow(grid)

  for (idx in seq_len(n_pts)) {
    # Standard-normal quadrature node z_std in R^d and product of 1D weights.
    z_std <- vapply(seq_len(d), function(j) nodes_1d[grid[idx, j]], 0)
    wgt_prod <- prod(vapply(seq_len(d), function(j) weights_1d[grid[idx, j]], 0))
    
    # Transform to (X, Z): (X, Z)' = l_chol z_std.
    xz <- l_chol %*% z_std
    x_part <- xz[seq_len(p), , drop = FALSE]
    z_part <- xz[p + seq_len(q), , drop = FALSE]

    # Linear predictor and probability pi = P(Y=1 | X, Z).
    eta <- as.numeric(beta_0 + crossprod(x_part, beta_x) + crossprod(z_part, beta_z))
    pi_i <- stats::plogis(eta)

    # Contribution to E[W W' pi(1-pi)]: weight = (product of 1D weights) * pi(1-pi).
    wgt <- wgt_prod * pi_i * (1 - pi_i)
    w <- matrix(c(1, as.numeric(x_part), as.numeric(z_part)), ncol = 1L)
    i_sum <- i_sum + wgt * (w %*% t(w))
  }

  return(i_sum)
}


#' Expected information matrix for logistic regression (Monte Carlo)
#'
#' Approximates \eqn{i_{\beta\beta'} = E[W W' \pi(1-\pi)]} by Monte Carlo.
#' Assumes \eqn{(X, Z)} is multivariate normal with mean zero and covariance
#' given by \code{var_x}, \code{var_z}, and \code{cov_xz}.
#'
#' @param beta_0 Intercept (scalar).
#' @param beta_x Coefficient for X (vector or column matrix).
#' @param beta_z Coefficient for Z (vector or column matrix).
#' @param var_x Variance of X (matrix).
#' @param var_z Variance of Z (matrix).
#' @param cov_xz Covariance of X and Z (matrix).
#' @param n_sim Number of Monte Carlo draws.
#' @return The expected information matrix (numeric matrix).
#' @noRd
.LogisticInformationMatrixMonteCarlo <- function(
    beta_0,
    beta_x,
    beta_z,
    var_x,
    var_z,
    cov_xz,
    n_sim = 100000L
) {
  # Dimensions and covariance of (X, Z); Cholesky so (X, Z)' = l_chol z_std.
  beta_x <- as.matrix(beta_x)
  beta_z <- as.matrix(beta_z)
  p <- nrow(beta_x)
  q <- nrow(beta_z)
  var_x <- as.matrix(var_x)
  var_z <- as.matrix(var_z)
  cov_xz <- as.matrix(cov_xz)

  sigma <- rbind(
    cbind(var_x, cov_xz),
    cbind(t(cov_xz), var_z)
  )
  l_chol <- t(chol(sigma))

  # Accumulator for sum of W W' pi(1-pi) over draws; order (intercept, X, Z).
  n_par <- 1L + p + q
  i_sum <- matrix(0, n_par, n_par)

  for (i in seq_len(n_sim)) {
    z_std <- stats::rnorm(p + q)
    xz <- l_chol %*% z_std
    x_part <- xz[seq_len(p), , drop = FALSE]
    z_part <- xz[p + seq_len(q), , drop = FALSE]
    eta <- as.numeric(beta_0 + crossprod(x_part, beta_x) + crossprod(z_part, beta_z))
    pi_i <- stats::plogis(eta)
    wgt <- pi_i * (1 - pi_i)
    w <- matrix(c(1, as.numeric(x_part), as.numeric(z_part)), ncol = 1L)
    i_sum <- i_sum + wgt * (w %*% t(w))
  }

  i_expected <- i_sum / n_sim
  return(i_expected)
}

#' Expected information matrix for logistic regression
#'
#' Dispatches to Gauss-Hermite quadrature or Monte Carlo depending on
#' \code{method}. Quadrature is used when \code{method = "quadrature"} (or
#' \code{"auto"} and dimension \eqn{d \leq 5}); otherwise Monte Carlo.
#'
#' @param beta_0 Intercept (scalar).
#' @param beta_x Coefficient for X (matrix).
#' @param beta_z Coefficient for Z (matrix).
#' @param var_x Variance of X (matrix).
#' @param var_z Variance of Z (matrix).
#' @param cov_xz Covariance of X and Z (matrix).
#' @param method \code{"quadrature"}, \code{"montecarlo"}, or \code{"auto"}.
#' @param n_quad Number of quadrature points per dimension (quadrature only).
#' @param n_sim Number of Monte Carlo draws (Monte Carlo only).
#' @return The expected information matrix (numeric matrix).
#' @noRd
.LogisticInformationMatrix <- function(
    beta_0,
    beta_x,
    beta_z,
    var_x,
    var_z,
    cov_xz,
    method = c("auto", "quadrature", "montecarlo"),
    n_quad = 15L,
    n_sim = 100000L
) {
  method <- match.arg(method)
  p <- nrow(as.matrix(beta_x))
  q <- nrow(as.matrix(beta_z))
  d <- p + q

  # Choose quadrature when explicitly requested or when auto and d is small.
  use_quadrature <- switch(
    method,
    auto = d <= 5L,
    quadrature = TRUE,
    montecarlo = FALSE
  )

  if (use_quadrature) {
    if (d > 6L) {
      warning(
        "Dimension of (X, Z) is ", d, "; quadrature may be slow. ",
        "Consider method = \"montecarlo\"."
      )
    }
    i_full <- .LogisticInformationMatrixQuadrature(
      beta_0, beta_x, beta_z, var_x, var_z, cov_xz, n_quad
    )
    return(i_full)
  }
  i_full <- .LogisticInformationMatrixMonteCarlo(
    beta_0, beta_x, beta_z, var_x, var_z, cov_xz, n_sim
  )
  return(i_full)
}

#' Efficient information for beta_X in logistic regression
#'
#' From the full expected information matrix, extracts the block for
#' \eqn{(\beta_X, \gamma)} with \eqn{\gamma = (\beta_0, \beta_Z)} and returns
#' the Schur complement (efficient information for \eqn{\beta_X}).
#'
#' @param i_full Full expected information matrix from
#'   \code{.LogisticInformationMatrix}, with order (intercept, X, Z).
#' @param p Number of X covariates (length of \code{beta_x}).
#' @param q Number of Z covariates (length of \code{beta_z}).
#' @return The efficient information matrix for \eqn{\beta_X} (p x p).
#' @noRd
.LogisticEfficientInformation <- function(i_full, p, q) {
  # Block indices: intercept (1), X (2:(1+p)), Z ((2+p):(1+p+q)).
  # Efficient information for beta_X is Schur complement of gamma = (intercept, Z) block.
  idx_x <- seq.int(2L, 1L + p)
  idx_gamma <- c(1L, seq.int(2L + p, 1L + p + q))
  i_bb <- i_full[idx_x, idx_x, drop = FALSE]
  i_aa <- i_full[idx_gamma, idx_gamma, drop = FALSE]
  i_ba <- i_full[idx_x, idx_gamma, drop = FALSE]
  i_eff <- SchurC(i_bb, i_aa, i_ba)
  return(i_eff)
}

#' Power from n and precomputed efficient information (logistic)
#'
#' Internal helper: Wald-test power for \eqn{H_0: \beta_X = 0} in the
#' logistic model, given sample size \code{n} and the efficient information
#' matrix for \eqn{\beta_X}.
#'
#' @param n Sample size (may be a vector).
#' @param beta_x True coefficient for X (matrix, e.g. column vector).
#' @param i_eff Efficient information for \eqn{\beta_X} (p x p matrix).
#' @param t1e Type I error level.
#' @return Numeric power (same length as \code{n}).
#' @noRd
.PowerLogisticFromInfo <- function(n, beta_x, i_eff, t1e) {
  # Non-centrality parameter Delta = n * beta_x' i_eff beta_x; Wald ~ chi^2_k(Delta).
  ncp <- n * as.numeric(matQF(beta_x, i_eff))
  k <- nrow(i_eff)
  crit <- stats::qchisq(p = 1 - t1e, df = k)
  power <- stats::pchisq(q = crit, df = k, ncp = ncp, lower.tail = FALSE)
  return(power)
}

# -----------------------------------------------------------------------------
# Public API
# -----------------------------------------------------------------------------

#' Power Calculation (Logistic Regression)
#'
#' Power for testing \eqn{H_0: \beta_X = 0} in the logistic regression model
#' \eqn{\mathrm{logit}(\pi) = \beta_0 + X'\beta_X + Z'\beta_Z}, where
#' \eqn{\pi = P(Y=1 \mid X, Z)}. The distribution of \eqn{(X, Z)} is assumed
#' multivariate normal with mean zero and covariance given by \code{var_x},
#' \code{var_z}, and \code{cov_xz}. The expected information
#' \eqn{E[W W' \pi(1-\pi)]} can be approximated by Gauss-Hermite quadrature
#' (default for small dimension) or Monte Carlo.
#'
#' @param beta_0 Intercept (scalar).
#' @param beta_x True regression coefficient for X (vector or matrix).
#' @param beta_z True regression coefficient for Z (vector or matrix).
#' @param var_x Variance of X (matrix).
#' @param var_z Variance of Z (matrix).
#' @param cov_xz Covariance of X and Z (matrix).
#' @param n Sample size.
#' @param t1e Type I error level.
#' @param method \code{"auto"} (quadrature when \eqn{\dim(X)+\dim(Z) \leq 5},
#'   else Monte Carlo), \code{"quadrature"}, or \code{"montecarlo"}.
#' @param n_quad Number of quadrature points per dimension (quadrature only).
#' @param n_sim Number of Monte Carlo draws (Monte Carlo only).
#' @export
#' @return Numeric power.
#'
#' @details Power is computed using the non-central chi-squared distribution
#'   under the asymptotic Wald test. With \code{method = "quadrature"},
#'   the expectation is approximated by product Gauss-Hermite quadrature
#'   over the normal distribution of \eqn{(X, Z)}, which is typically more
#'   accurate and faster than Monte Carlo when the dimension is small.
#'
#' @seealso \code{\link{SampleSizeLogisticReg}} for sample size given a target
#'   power.
PowerLogisticReg <- function(
  beta_0,
  beta_x,
  beta_z,
  var_x,
  var_z,
  cov_xz,
  n,
  t1e = 0.05,
  method = c("auto", "quadrature", "montecarlo"),
  n_quad = 15L,
  n_sim = 100000L
) {
  beta_x <- as.matrix(beta_x)
  beta_z <- as.matrix(beta_z)
  var_x <- as.matrix(var_x)
  var_z <- as.matrix(var_z)
  cov_xz <- as.matrix(cov_xz)

  .ValidatePositiveN(n)
  .ValidateLogisticDesign(beta_0, beta_x, beta_z, var_x, var_z, cov_xz)

  # Under the null, power equals type I error.
  if (all(beta_x == 0)) {
    return(t1e)
  }

  # Expected information, then efficient information for beta_X.
  i_full <- .LogisticInformationMatrix(
    beta_0, beta_x, beta_z, var_x, var_z, cov_xz,
    method = method, n_quad = n_quad, n_sim = n_sim
  )
  p <- nrow(beta_x)
  q <- nrow(beta_z)
  i_eff <- .LogisticEfficientInformation(i_full, p, q)

  power <- .PowerLogisticFromInfo(n, beta_x, i_eff, t1e)
  return(power)
}

#' Sample Size Calculation (Logistic Regression)
#'
#' Sample size for testing \eqn{H_0: \beta_X = 0} in the logistic regression
#' model \eqn{\mathrm{logit}(\pi) = \beta_0 + X'\beta_X + Z'\beta_Z}. Finds the
#' smallest integer \eqn{n} in \eqn{[\code{min_n}, \code{max_n}]} such that
#' power reaches the target. The expected information is approximated once
#' (by quadrature or Monte Carlo) and reused in the root-finding search.
#'
#' @param beta_0 Intercept (scalar).
#' @param beta_x True regression coefficient for X (vector or matrix).
#' @param beta_z True regression coefficient for Z (vector or matrix).
#' @param var_x Variance of X (matrix).
#' @param var_z Variance of Z (matrix).
#' @param cov_xz Covariance of X and Z (matrix).
#' @param max_n Maximum allowable sample size (upper bound for search).
#' @param power Target power.
#' @param t1e Type I error level.
#' @param min_n Minimum allowable sample size (lower bound for search).
#' @param tol Numeric tolerance for \code{\link[stats]{uniroot}}.
#' @param method \code{"auto"}, \code{"quadrature"}, or \code{"montecarlo"}.
#' @param n_quad Number of quadrature points per dimension (quadrature only).
#' @param n_sim Number of Monte Carlo draws (Monte Carlo only).
#' @export
#' @return A data frame with columns \code{n} (sample size) and \code{power}
#'   (achieved power at that \code{n}).
#'
#' @details Uses \code{\link[stats]{uniroot}} to find the smallest sample size
#' at which power reaches the target. The efficient information is computed
#' once and held fixed during the search.
#'
#' @seealso \code{\link{PowerLogisticReg}} for power at a given sample size.
SampleSizeLogisticReg <- function(
  beta_0,
  beta_x,
  beta_z,
  var_x,
  var_z,
  cov_xz,
  max_n,
  power = 0.80,
  t1e = 0.05,
  min_n = 1,
  tol = .Machine$double.eps^0.25,
  method = c("auto", "quadrature", "montecarlo"),
  n_quad = 15L,
  n_sim = 100000L
) {
  beta_x <- as.matrix(beta_x)
  beta_z <- as.matrix(beta_z)
  var_x <- as.matrix(var_x)
  var_z <- as.matrix(var_z)
  cov_xz <- as.matrix(cov_xz)

  .ValidateSampleSizeBounds(min_n, max_n)
  .ValidateTargetPower(power)
  .ValidateLogisticDesign(beta_0, beta_x, beta_z, var_x, var_z, cov_xz)
  .ValidateNonZeroBeta(beta_x)

  # Expected information and efficient information for beta_X (computed once).
  i_full <- .LogisticInformationMatrix(
    beta_0, beta_x, beta_z, var_x, var_z, cov_xz,
    method = method, n_quad = n_quad, n_sim = n_sim
  )
  p <- nrow(beta_x)
  q <- nrow(beta_z)
  i_eff <- .LogisticEfficientInformation(i_full, p, q)

  # Objective for uniroot: power(n) - target.
  obj <- function(n) {
    delta <- .PowerLogisticFromInfo(n, beta_x, i_eff, t1e) - power
    return(delta)
  }

  min_power <- .PowerLogisticFromInfo(min_n, beta_x, i_eff, t1e)
  max_power <- .PowerLogisticFromInfo(max_n, beta_x, i_eff, t1e)
  if (max_power < power) {
    stop(
      "Target power ", round(power, 3), " is unattainable with max_n = ", max_n,
      ". Achieved power at max_n is ", round(max_power, 3), ". ",
      "Increase max_n and try again."
    )
  }

  # Smallest n in [min_n, max_n] achieving target power.
  if (min_power >= power) {
    final_n <- min_n
  } else {
    zero <- stats::uniroot(f = obj, lower = min_n, upper = max_n, tol = tol)
    final_n <- ceiling(zero$root)
  }
  final_power <- .PowerLogisticFromInfo(final_n, beta_x, i_eff, t1e)

  out <- data.frame(
    n = final_n,
    power = final_power
  )
  return(out)
}
