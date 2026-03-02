test_that("PowerLinReg returns a numeric in [0, 1]", {
  p <- PowerLinReg(
    beta_x = 1,
    cov_xz = 0.5,
    n = 10,
    t1e = 0.05,
    var_resid = 1,
    var_x = 1,
    var_z = 1
  )
  expect_type(p, "double")
  expect_length(p, 1L)
  expect_gte(p, 0)
  expect_lte(p, 1)
})

test_that("PowerLinReg matches known result from README", {
  p <- PowerLinReg(
    beta_x = 1,
    cov_xz = 0.5,
    n = 10,
    t1e = 0.05,
    var_resid = 1,
    var_x = 1,
    var_z = 1
  )
  expect_equal(p, 0.781908, tolerance = 1e-5)
})

test_that("PowerLinReg under null gives power equal to type I error", {
  p <- PowerLinReg(
    beta_x = 0,
    cov_xz = 0,
    n = 100,
    t1e = 0.05,
    var_resid = 1,
    var_x = 1,
    var_z = 1
  )
  expect_equal(p, 0.05, tolerance = 1e-6)
})

test_that("PowerLinReg increases with n", {
  p_small <- PowerLinReg(
    beta_x = 1,
    cov_xz = 0,
    n = 20,
    t1e = 0.05,
    var_resid = 1,
    var_x = 1,
    var_z = 1
  )
  p_large <- PowerLinReg(
    beta_x = 1,
    cov_xz = 0,
    n = 200,
    t1e = 0.05,
    var_resid = 1,
    var_x = 1,
    var_z = 1
  )
  expect_gt(p_large, p_small)
})

test_that("PowerLinReg accepts matrix inputs for multivariate X, Z", {
  beta_x <- matrix(c(0.5, 0.5), ncol = 1)
  var_x <- diag(2)
  var_z <- matrix(1)
  cov_xz <- matrix(c(0.2, 0.2), nrow = 2, ncol = 1)
  p <- PowerLinReg(
    beta_x = beta_x,
    cov_xz = cov_xz,
    n = 50,
    t1e = 0.05,
    var_resid = 1,
    var_x = var_x,
    var_z = var_z
  )
  expect_type(p, "double")
  expect_gte(p, 0)
  expect_lte(p, 1)
})
