library(testthat)
library(tlmixture)

context("cluster-aware fold helpers")

set.seed(21, "L'Ecuyer-CMRG")
n_fam = 150
fam_size = sample(1:3, n_fam, replace = TRUE)
family = rep(seq_len(n_fam), fam_size)
n = length(family)
# Rare-ish binary outcome that clusters within families.
y = rbinom(n, 1, plogis(-1.5 + rnorm(n_fam, 0, 1.5)[family]))
dat = data.frame(x = rnorm(n), y = y)

intact = function(fold, g) {
  all(tapply(fold, g, function(x) length(unique(x))) == 1L)
}

fold_vector = function(rset) {
  fold = integer(nrow(dat))
  for (k in seq_along(rset$splits)) {
    fold[rsample::complement(rset$splits[[k]])] = k
  }
  fold
}

test_that("cluster_vfold_cv keeps clusters whole and never exposes the grouping", {
  set.seed(1)
  rset = tlmixture:::cluster_vfold_cv(dat, family, v = 5L, strata_values = dat$y)
  expect_length(rset$splits, 5L)
  fold = fold_vector(rset)
  expect_true(all(fold > 0L))
  expect_true(intact(fold, family))
  expect_identical(names(rsample::analysis(rset$splits[[1]])), names(dat))
})

test_that("stratifying on the cluster mean outcome keeps folds balanced", {
  set.seed(2)
  rset = tlmixture:::cluster_vfold_cv(dat, family, v = 5L, strata_values = dat$y)
  rates = tapply(dat$y, fold_vector(rset), mean)
  # Every fold has events, and no fold's event rate is far from the overall.
  expect_true(all(rates > 0))
  expect_lt(max(abs(rates - mean(dat$y))), 0.15)
})

test_that("cluster_vfold_cv falls back to unstratified folds when needed", {
  set.seed(3)
  expect_silent(rset <- tlmixture:::cluster_vfold_cv(dat, family, v = 2L,
                                                     strata_values = rep(1, n)))
  expect_true(intact(fold_vector(rset), family))
  expect_error(tlmixture:::cluster_vfold_cv(dat, family[-1], v = 2L),
               "cluster grouping has")
  expect_error(tlmixture:::cluster_vfold_cv(dat, rep(1:2, length.out = n), v = 3L),
               "fewer than folds_cvtmle")
})

test_that("crossfit_folds is unchanged without ids and keeps clusters with ids", {
  set.seed(4)
  old = sample(rep(seq(5L), length.out = n))
  set.seed(4)
  expect_identical(tlmixture:::crossfit_folds(n, 5L)$folds, old)

  cf = tlmixture:::crossfit_folds(n, 5L, family)
  expect_identical(cf$cv_folds, 5L)
  expect_true(intact(cf$folds, family))
  # Fewer clusters than folds: one fold per cluster.
  few = tlmixture:::crossfit_folds(9L, 5L, rep(1:3, each = 3))
  expect_identical(few$cv_folds, 3L)
})

test_that("sl_cluster_args caps folds and drops ids that cannot be used", {
  f = tlmixture:::sl_cluster_args
  expect_identical(f(NULL, 5L), list(id = NULL, V = 5L))
  expect_identical(f(family, 5L)$V, 5L)
  expect_identical(f(rep(1:3, 10), 5L)$V, 3L)
  expect_identical(f(rep(1, 10), 5L), list(id = NULL, V = 2L))
})

test_that("accepts_arg checks for the argument or dots", {
  expect_true(tlmixture:::accepts_arg(function(id) NULL, "id"))
  expect_true(tlmixture:::accepts_arg(function(...) NULL, "id"))
  expect_false(tlmixture:::accepts_arg(function(data) NULL, "id"))
})
