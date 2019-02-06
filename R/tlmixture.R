#' Targeted learning for exposure mixtures
#'
#' This is our main function.
#'
#' @param data Data frame with outcome, exposure, and adjustment variables.
#' @param outcome Name of the outcome variable.
#' @param exposures A vector of exposure names, or (not yet supported) a list where each element is
#' a vector of pre-clustered exposures.
#' @param quantiles_mixtures Number of quantiles to use for discretizing mixture
#' (default 3 - low, medium, high).
#' @param quantiles_exposures Number of quantiles to use for discretizing continuous exposures
#' (default 4).
#'
#' @param folds_cvtmle Number of CV-TMLE folds (default 2).
# TODO: support arbitrary estimators (sl3, mlr, etc.)
#' @param estimator_outcome SuperLearner library for outcome estimation.
# TODO: support arbitrary estimators (sl3, mlr, etc.)
#' @param estimator_propensity SuperLearner library for propensity estimation.
#' @param cluster_exposures Whether to automatically cluster a vector of exposures into
#' sub-groups (default FALSE, TRUE not yet supported).
#' @param verbose If TRUE, display more detailed info during execution.
#'
#' @export
tlmixture =
  function(data, outcome, exposures,
           quantiles_mixtures = 3L,
           quantiles_exposures = 4L,
           folds_cvtmle = 2L,
           estimator_outcome = c("SL.mean", "SL.glmnet"),
           estimator_propensity = estimator_outcome,
           cluster_exposures = FALSE,
           verbose = FALSE
           ) {

  ##################
  # Initial setup

  # Detect outcome type: "binary" if outcome is {0, 1}, "continuous" otherwise.
  family = outcome_family(data, outcome)

  if (verbose) {
    cat("Outcome family:", family, "\n")
  }

  # Setup cross-validation folds, stratified on the outcome.
  folds = rsample::vfold_cv(data, v = folds_cvtmle, strata = outcome)


  ##################
  # Loop over folds, analyzing training sets and then applying to test.

  fold_results =
    purrr::map(folds$splits, analyze_folds, outcome = outcome, exposures = exposures,
               family = family,
               quantiles_mixtures = quantiles_mixtures,
               quantiles_exposures = quantiles_exposures,
               estimator_outcome = estimator_outcome,
               estimator_propensity = estimator_propensity,
               cluster_exposures = cluster_exposures,
               verbose = verbose)


  ##################
  # Combine test set results (CV-TMLE etc.)

  combined_results =
    combine_test_results(fold_results, family = family,
                         quantiles_mixtures = quantiles_mixtures,
                         verbose = verbose)


  ##################
  # Compile everything we want to return.
  results =
    list(folds = fold_results,
         combined = combined_results,
         folds = folds,
         family = family)

  return(results)
}
