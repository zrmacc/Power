# Purpose: Testing power calculations.
library(Power)

beta_x <- sqrt(0.05)
cov_xz <- 0.5
n <- 10
max_n <- 1e3
t1e <- 0.05
var_resid <- 1
var_x <- 1
var_z <- 1

PowerLinReg(
  beta_x = beta_x,
  cov_xz = cov_xz,
  n = n,
  t1e = 0.05,
  var_resid = var_resid,
  var_x = var_x,
  var_z = var_z
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
