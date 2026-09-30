library(testthat)
library(tlmixture)

context("tlmixture() with cluster-robust standard errors")

# Families (1-3 children) nested in sites, with a shared family effect.
set.seed(5, "L'Ecuyer-CMRG")
n_fam = 120
fam_size = sample(1:3, n_fam, replace = TRUE)
family = rep(seq_len(n_fam), fam_size)
n = length(family)
site = rep(sample(1:10, n_fam, replace = TRUE), fam_size)
u = rnorm(n_fam, 0, 2)[family]
X = data.frame(w1 = rnorm(n), w2 = rbinom(n, 1, 0.5),
               a1 = rnorm(n_fam)[family] + rnorm(n, 0, 0.5),
               a2 = rnorm(n), a3 = rnorm(n))
X$y = 50 + 2 * X$a1 + X$a2 + X$w1 + u + rnorm(n, 0, 3)
exposures = c("a1", "a2", "a3")

run = function(data, ...) {
  suppressWarnings(
    tlmixture(data, "y", exposures,
              estimator_outcome = c("SL.mean", "SL.glm"),
              mixture_fn = mixture_glm,
              refit_mixtures = FALSE, ...))
}

fold_of = function(result) {
  folds = integer(nrow(X))
  for (k in seq_along(result$fold_sets$splits)) {
    folds[rsample::complement(result$fold_sets$splits[[k]])] = k
  }
  folds
}

test_that("an id with one observation per cluster changes nothing", {
  set.seed(1)
  base = run(X)
  set.seed(1)
  expect_message(singleton <- run(X, id = seq_len(n)),
                 "one observation per cluster")
  expect_identical(singleton$combined, base$combined)
  expect_null(singleton$cluster_id)
  expect_identical(singleton$cluster_df, Inf)
})

test_that("clustering by family keeps families within folds and uses t df", {
  set.seed(2)
  res = run(X, id = family)
  expect_equal(res$cluster_df, n_fam - 1)
  folds = fold_of(res)
  expect_true(all(folds > 0))
  expect_true(all(tapply(folds, family, function(x) length(unique(x))) == 1L))
  expect_true(is.data.frame(res$combined$groups_df))
  expect_true(is.finite(res$combined$groups_df$rd_se))
})

test_that("id can name columns of data, which are not used as confounders", {
  Xf = cbind(X, family = family)
  set.seed(3)
  by_name = run(Xf, id = "family")
  set.seed(3)
  by_vector = run(X, id = family)
  expect_identical(by_name$combined, by_vector$combined)
  split_data = rsample::analysis(by_name$fold_sets$splits[[1]])
  expect_false("family" %in% names(split_data))
})

test_that("nested family and site reduce to site", {
  Xfs = cbind(X, family = family, site = site)
  set.seed(4)
  expect_message(res <- run(Xfs, id = c("family", "site")),
                 "nested in 'site'")
  expect_identical(names(res$cluster_id), "site")
  expect_equal(res$cluster_df, length(unique(site)) - 1)
  folds = fold_of(res)
  expect_true(all(tapply(folds, site, function(x) length(unique(x))) == 1L))
})

test_that("mixture_fn receives the training ids when it accepts them", {
  got = new.env()
  spy = function(data, ...) {
    args = list(...)
    got$id = args$id
    got$n = nrow(data)
    mixture_glm(data, ...)
  }
  set.seed(5)
  suppressWarnings(
    tlmixture(X, "y", exposures, estimator_outcome = c("SL.mean", "SL.glm"),
              mixture_fn = spy, refit_mixtures = FALSE, id = family))
  expect_false(is.null(got$id))
  expect_length(got$id, got$n)

  # Without clustering, no id is passed.
  got$id = "unset"
  set.seed(5)
  suppressWarnings(
    tlmixture(X, "y", exposures, estimator_outcome = c("SL.mean", "SL.glm"),
              mixture_fn = spy, refit_mixtures = FALSE))
  expect_null(got$id)
})

test_that("a binary outcome runs with clustering", {
  Xb = X
  Xb$y = rbinom(n, 1, plogis(X$a1 + u / 2))
  set.seed(6)
  res = run(Xb, id = family)
  expect_true(is.finite(res$combined$groups_df$rd_se))
})

test_that("mixture_score_sl3 accepts id and keeps it out of the confounders", {
  skip_if_not_installed("sl3")
  lrnr = sl3::Lrnr_sl$new(
    learners = sl3::Stack$new(sl3::Lrnr_glm$new(), sl3::Lrnr_mean$new()),
    metalearner = sl3::Lrnr_nnls$new(convex = TRUE))
  set.seed(7)
  mix = mixture_score_sl3(X, "y", exposures, id = family,
                          estimator_mixture = lrnr)
  expect_setequal(mix$confounders, c("w1", "w2"))
  expect_length(suppressWarnings(predict(mix, X)), n)
  expect_error(mixture_score_sl3(X, "y", exposures, id = family[-1],
                                 estimator_mixture = lrnr),
               "id has")
})

test_that("crossed ids run multiway, fall back with a warning, and obey folds_by", {
  school = sample(1:15, n, replace = TRUE)
  Xc = cbind(X, site = site, school = school)
  set.seed(8)
  expect_warning(
    res <- tlmixture(Xc, "y", exposures,
                     estimator_outcome = c("SL.mean", "SL.glm"),
                     mixture_fn = mixture_glm, refit_mixtures = FALSE,
                     id = c("site", "school")),
    "crossed")
  expect_identical(names(res$cluster_id), c("site", "school"))
  expect_equal(res$cluster_df,
               min(length(unique(site)), length(unique(school))) - 1)
  # Fallback: the coarsest variable with enough clusters (site) is intact.
  expect_true(all(tapply(fold_of(res), site,
                         function(x) length(unique(x))) == 1L))

  set.seed(9)
  res2 = run(Xc, id = c("site", "school"), folds_by = "school")
  expect_true(all(tapply(fold_of(res2), school,
                         function(x) length(unique(x))) == 1L))
  expect_error(run(Xc, id = c("site", "school"), folds_by = "family"),
               "not one of the cluster variables")
})

test_that("the outcome and propensity SuperLearners receive the cluster ids", {
  calls = list()
  spy = function(...) {
    args = list(...)
    calls[[length(calls) + 1L]] <<- list(id = args$id, cv = args$cvControl)
    SuperLearner::SuperLearner(...)
  }
  local({
    local_mocked_bindings(SuperLearner = spy, .package = "tlmixture")
    set.seed(10)
    run(X, id = family)
  })
  expect_gt(length(calls), 0L)
  for (call in calls) {
    expect_false(is.null(call$id))
    expect_false(call$cv$stratifyCV)
    # Each call's ids are training-fold family ids.
    expect_true(all(call$id %in% family))
  }

  # Without clustering, no ids are passed and stratification is as before.
  calls = list()
  local({
    local_mocked_bindings(SuperLearner = spy, .package = "tlmixture")
    set.seed(10)
    run(X)
  })
  expect_gt(length(calls), 0L)
  for (call in calls) expect_null(call$id)
})

test_that("the full-data refit passes the ids to mixture_fn", {
  got = new.env()
  got$refit_id = "unset"
  spy = function(data, ...) {
    args = list(...)
    if (nrow(data) == n) got$refit_id = args$id
    mixture_glm(data, ...)
  }
  set.seed(11)
  suppressWarnings(
    tlmixture(X, "y", exposures, estimator_outcome = c("SL.mean", "SL.glm"),
              mixture_fn = spy, refit_mixtures = TRUE, id = family))
  expect_length(got$refit_id, n)
  expect_true(all(tapply(got$refit_id, family,
                         function(x) length(unique(x))) == 1L))
})

test_that("mixture functions with internal CV keep clusters together", {
  # mixture_sl: SuperLearner's internal folds keep families whole.
  mix_sl = mixture_sl(X, "y", exposures, family = "gaussian", id = family)
  splits_family = vapply(mix_sl$sl$validRows, function(rows) {
    length(intersect(family[rows], family[-rows])) > 0L
  }, logical(1))
  expect_false(any(splits_family))

  # mixture_score_glm and mixture_backfit_sl accept id and run.
  set.seed(12)
  expect_s3_class(suppressWarnings(
    mixture_score_glm(X, "y", exposures, id = family)), "mixture_score_glm")
  expect_error(mixture_score_glm(X, "y", exposures, id = family[-1]), "id has")
  skip_if_not_installed("quadprog")
  set.seed(13)
  expect_s3_class(suppressWarnings(
    mixture_backfit_sl(X, "y", exposures, id = family, max_iterations = 2L,
                       estimator_mixture = c("SL.mean", "SL.glm"))),
    "mixture_backfit_sl")
})

test_that("sl3 mixture functions pass the ids to their sl3 tasks", {
  skip_if_not_installed("sl3")
  lrnr = function() {
    sl3::Lrnr_sl$new(
      learners = sl3::Stack$new(sl3::Lrnr_glm$new(), sl3::Lrnr_mean$new()),
      metalearner = sl3::Lrnr_nnls$new(convex = TRUE))
  }
  set.seed(14)
  score = suppressWarnings(
    mixture_score_sl3(X, "y", exposures, id = family,
                      estimator_mixture = lrnr()))
  task = score$reg_mixture$training_task
  expect_identical(task$nodes$id, ".tlmixture_cluster_id")
  split_any = vapply(task$folds, function(fold) {
    length(intersect(family[fold$training_set],
                     family[fold$validation_set])) > 0L
  }, logical(1))
  expect_false(any(split_any))

  set.seed(15)
  backfit = suppressWarnings(
    mixture_backfit_sl3(X, "y", exposures, id = family, max_iterations = 2L,
                        estimator_mixture = lrnr()))
  expect_identical(backfit$reg_mixture$training_task$nodes$id,
                   ".tlmixture_cluster_id")
  expect_identical(backfit$reg_adjust$training_task$nodes$id,
                   ".tlmixture_cluster_id")
})

test_that("bad id input is reported clearly", {
  expect_error(tlmixture(X, "y", exposures, id = "familly"),
               "not found in data: familly")
  expect_error(tlmixture(X, "y", exposures, id = "a1"),
               "outcome or an exposure")
  expect_error(tlmixture(X, "y", exposures, id = family[-1]),
               "rows but there are")
  # A supplied id that is NULL (e.g. a misspelled data$column) warns.
  expect_warning(
    tryCatch(tlmixture(X, "y", exposures, id = X$famly,
                       folds_cvtmle = 0L),
             error = function(e) NULL),
    "NOT clustered")
})
