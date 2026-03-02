test_that("PowerLogisticReg returns numeric in [0, 1]", {
  p <- PowerLogisticReg(
    beta_0 = 0,
    beta_x = 0.5,
    beta_z = 0.2,
    var_x = 1,
    var_z = 1,
    cov_xz = 0.3,
    n = 100,
    t1e = 0.05,
    n_sim = 5000L
  )
  expect_type(p, "double")
  expect_length(p, 1L)
  expect_gte(p, 0)
  expect_lte(p, 1)
})

test_that("PowerLogisticReg under null returns type I error", {
  p <- PowerLogisticReg(
    beta_0 = 0,
    beta_x = 0,
    beta_z = 0.2,
    var_x = 1,
    var_z = 1,
    cov_xz = 0,
    n = 200,
    t1e = 0.05,
    n_sim = 5000L
  )
  expect_equal(p, 0.05, tolerance = 0.01)
})

test_that("PowerLogisticReg increases with n", {
  p_small <- PowerLogisticReg(
    beta_0 = 0,
    beta_x = 0.4,
    beta_z = 0,
    var_x = 1,
    var_z = 1,
    cov_xz = 0,
    n = 50,
    t1e = 0.05,
    n_sim = 5000L
  )
  p_large <- PowerLogisticReg(
    beta_0 = 0,
    beta_x = 0.4,
    beta_z = 0,
    var_x = 1,
    var_z = 1,
    cov_xz = 0,
    n = 300,
    t1e = 0.05,
    n_sim = 5000L
  )
  expect_gt(p_large, p_small)
})

test_that("SampleSizeLogisticReg returns data.frame with n and power", {
  out <- SampleSizeLogisticReg(
    beta_0 = 0,
    beta_x = 0.5,
    beta_z = 0.2,
    var_x = 1,
    var_z = 1,
    cov_xz = 0.2,
    max_n = 500,
    power = 0.80,
    t1e = 0.05,
    n_sim = 5000L
  )
  expect_s3_class(out, "data.frame")
  expect_named(out, c("n", "power"))
  expect_gte(out$power, 0.75)
  expect_lte(out$power, 1)
  expect_true(out$n >= 1)
  expect_true(out$n <= 500)
})

test_that("SampleSizeLogisticReg is consistent with PowerLogisticReg", {
  out <- SampleSizeLogisticReg(
    beta_0 = 0,
    beta_x = 0.5,
    beta_z = 0,
    var_x = 1,
    var_z = 1,
    cov_xz = 0,
    max_n = 400,
    power = 0.85,
    t1e = 0.05,
    n_sim = 5000L
  )
  p_check <- PowerLogisticReg(
    beta_0 = 0,
    beta_x = 0.5,
    beta_z = 0,
    var_x = 1,
    var_z = 1,
    cov_xz = 0,
    n = out$n,
    t1e = 0.05,
    n_sim = 5000L
  )
  expect_equal(out$power, p_check, tolerance = 0.02)
})

test_that("SampleSizeLogisticReg errors when target unattainable", {
  expect_error(
    SampleSizeLogisticReg(
      beta_0 = 0,
      beta_x = 0.1,
      beta_z = 0,
      var_x = 1,
      var_z = 1,
      cov_xz = 0,
      max_n = 20,
      power = 0.90,
      n_sim = 5000L
    ),
    "unattainable"
  )
})

test_that("Gauss-Hermite quadrature and Monte Carlo give similar power", {
  p_quad <- PowerLogisticReg(
    beta_0 = 0,
    beta_x = 0.5,
    beta_z = 0.2,
    var_x = 1,
    var_z = 1,
    cov_xz = 0.2,
    n = 150,
    t1e = 0.05,
    method = "quadrature",
    n_quad = 20L
  )
  p_mc <- PowerLogisticReg(
    beta_0 = 0,
    beta_x = 0.5,
    beta_z = 0.2,
    var_x = 1,
    var_z = 1,
    cov_xz = 0.2,
    n = 150,
    t1e = 0.05,
    method = "montecarlo",
    n_sim = 50000L
  )
  expect_equal(p_quad, p_mc, tolerance = 0.03)
})
