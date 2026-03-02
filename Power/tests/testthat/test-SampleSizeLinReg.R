test_that("SampleSizeLinReg returns data.frame with n and power", {
  out <- SampleSizeLinReg(
    beta_x = 1,
    cov_xz = 0.5,
    max_n = 100,
    power = 0.90,
    t1e = 0.05,
    var_resid = 1,
    var_x = 1,
    var_z = 1
  )
  expect_s3_class(out, "data.frame")
  expect_named(out, c("n", "power"))
  expect_type(out$n, "double")
  expect_type(out$power, "double")
  expect_equal(nrow(out), 1L)
})

test_that("SampleSizeLinReg n is integer and power in [0, 1]", {
  out <- SampleSizeLinReg(
    beta_x = 1,
    cov_xz = 0.5,
    max_n = 100,
    power = 0.90,
    t1e = 0.05,
    var_resid = 1,
    var_x = 1,
    var_z = 1
  )
  expect_equal(out$n, round(out$n))
  expect_gte(out$power, 0)
  expect_lte(out$power, 1)
})

test_that("SampleSizeLinReg achieved power meets target", {
  out <- SampleSizeLinReg(
    beta_x = 1,
    cov_xz = 0.5,
    max_n = 100,
    power = 0.90,
    t1e = 0.05,
    var_resid = 1,
    var_x = 1,
    var_z = 1
  )
  expect_gte(out$power, 0.90)
})

test_that("SampleSizeLinReg is consistent with PowerLinReg", {
  out <- SampleSizeLinReg(
    beta_x = 1,
    cov_xz = 0.5,
    max_n = 100,
    power = 0.90,
    t1e = 0.05,
    var_resid = 1,
    var_x = 1,
    var_z = 1
  )
  p_check <- PowerLinReg(
    beta_x = 1,
    cov_xz = 0.5,
    n = out$n,
    t1e = 0.05,
    var_resid = 1,
    var_x = 1,
    var_z = 1
  )
  expect_equal(out$power, p_check, tolerance = 1e-6)
})

test_that("SampleSizeLinReg matches known result from README", {
  out <- SampleSizeLinReg(
    beta_x = 1,
    cov_xz = 0.5,
    max_n = 100,
    power = 0.90,
    t1e = 0.05,
    var_resid = 1,
    var_x = 1,
    var_z = 1
  )
  expect_equal(out$n, 15)
  expect_equal(out$power, 0.9183621, tolerance = 1e-5)
})

test_that("SampleSizeLinReg errors when target power not attainable", {
  expect_error(
    SampleSizeLinReg(
      beta_x = 0.01,
      cov_xz = 0,
      max_n = 10,
      power = 0.90,
      t1e = 0.05,
      var_resid = 1,
      var_x = 1,
      var_z = 1
    ),
    "unattainable"
  )
})
