library(testthat)
library(tlmixture)

context("standardize()")

test_that("standardize() scales all but the skipped columns, keeping order", {
  df = data.frame(a = c(1, 2, 3, 4), y = c(10, 20, 30, 40), b = c(2, 4, 6, 8))
  out = tlmixture:::standardize(df, skip_vars = "y")
  expect_identical(names(out), names(df))
  expect_identical(out$y, df$y)
  expect_equal(mean(out$a), 0)
  expect_equal(sd(out$a), 1)
  expect_equal(out$b, as.numeric(scale(df$b)))
})
