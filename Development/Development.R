# Purpose: Testing power calculations.
library(Power)

bx <- sqrt(0.005)
bz <- sqrt(0.050)
rho <- 0.50
pve <- bx^2 + bz^2 + 2 * bx * bz * rho

PowerLinReg(
  beta_x = bx,
  cov_xz = rho,
  n = 1000,
  t1e = 0.05,
  var_resid = 1 - pve,
  var_x = 1,
  var_z = 1
)

SampleSizeLinReg(
  beta_x = beta_x,
  cov_xz = cov_xz,
  max_n = 1000,
  t1e = 0.05,
  power = 0.80,
  var_resid = var_resid,
  var_x = var_x,
  var_z = var_z
)
