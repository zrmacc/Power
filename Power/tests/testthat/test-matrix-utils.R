test_that("matQF computes x' A x correctly", {
  x <- matrix(c(1, 2), ncol = 1)
  A <- matrix(c(2, 0, 0, 3), nrow = 2, ncol = 2)
  out <- matQF(x, A)
  expect_equal(as.numeric(out), 2 * 1^2 + 3 * 2^2)
})

test_that("matQF is scalar for 1x1 A and vector x", {
  x <- matrix(1, ncol = 1)
  A <- matrix(2, nrow = 1, ncol = 1)
  out <- matQF(x, A)
  expect_equal(dim(out), c(1L, 1L))
  expect_equal(as.numeric(out), 2)
})

test_that("SchurC matches manual block formula for 2x2 case", {
  Ibb <- matrix(2, nrow = 1, ncol = 1)
  Iaa <- matrix(1, nrow = 1, ncol = 1)
  Iba <- matrix(0.5, nrow = 1, ncol = 1)
  out <- SchurC(Ibb, Iaa, Iba)
  # I_bb - I_ba I_aa^{-1} I_ab = 2 - 0.5 * 1 * 0.5 = 1.75
  expect_equal(as.numeric(out), 2 - 0.5 * 0.5, tolerance = 1e-10)
})

test_that("SchurC returns symmetric matrix for symmetric inputs", {
  Ibb <- matrix(c(2, 0.5, 0.5, 2), nrow = 2, ncol = 2)
  Iaa <- matrix(1, nrow = 1, ncol = 1)
  Iba <- matrix(c(0.3, 0.3), nrow = 2, ncol = 1)
  out <- SchurC(Ibb, Iaa, Iba)
  expect_equal(out, t(out), tolerance = 1e-10)
})
