test_that("package exports required user-facing functions", {
  expect_true(is.function(glmm_rr))
  expect_true(is.function(summary.glmm_rr))
  expect_true(is.function(plot.glmm_rr))
})

test_that("static helpers compute safe probit values", {
  expect_true(is.finite(glmmRRdiag:::.safe_qnorm(0.5)))
})
