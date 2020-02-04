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
#' @param folds_sl Number of SL folds during outcome and propensity estimation.
# TODO: support arbitrary estimators (sl3, mlr, etc.)
#' @param estimator_outcome SuperLearner library for outcome estimation.
# TODO: support arbitrary estimators (sl3, mlr, etc.)
#' @param estimator_propensity SuperLearner library for propensity estimation.
#' @param cluster_exposures Whether to automatically cluster a vector of exposures into
#' sub-groups (default FALSE; TRUE not yet supported).
#' @param mixture_fn Current options: mixture_glm, mixture_pls, or mixture_sl
#' @param refit_mixtures After CV-TMEL, refit mixture functions to full dataset.
#' @param verbose If TRUE, display more detailed info during execution.
#'
#' @export
tlmixture =
  function(data, outcome, exposures,
           quantiles_mixtures = 3L,
           quantiles_exposures = 4L,
           folds_cvtmle = 2L,
           folds_sl = 2L,
           estimator_outcome = c("SL.mean", "SL.glmnet"),
           estimator_propensity = estimator_outcome,
           cluster_exposures = FALSE,
           mixture_fn = mixture_glm,
           refit_mixtures = TRUE,
           verbose = FALSE
           ) {

  ##################
  # Initial setup

  # Detect outcome type: "binary" if outcome is {0, 1}, "continuous" otherwise.
  family = outcome_family(data, outcome)

  if (verbose) {
    cat("Outcome family:", family, "\n")
  }


  # TODO: move this rescaling into its own function.

  # Save bounds on the full Y variables for later transformation if Y is not binary.
  # TODO: review tmle3 to see how it handles this rescaling.
  if (family == "binomial" || length(unique(data[[outcome]])) == 2L) {
    q_bounds = c(0, 1)
    needs_rescale = FALSE
  } else {
    # This part is duplicated from the TMLE code in tmle_init_stage1.

    # Define Qbounds just for continuous (non-binary) outcomes.
    q_bounds = range(data[[outcome]], na.rm = TRUE)
    # Extend bounds 10% beyond the observed range.
    # NOTE: if one of the bounds is zero then it won't be extended.
    q_bounds = q_bounds + 0.1 * c(-abs(q_bounds[1]), abs(q_bounds[2]))
    needs_rescale = TRUE
  }

  # Bound and transform outcome.
  if (needs_rescale) {
    #oldY = data[[outcome]]
    y_bounded = bound(data[[outcome]], q_bounds)

    outcome_range = range(y_bounded, na.rm = TRUE)
    # Ystar[is.na(Ystar)] <- 0

    # This rescales the outcome to be \in [0, 1]
    y_rescaled = (y_bounded - outcome_range[1]) / diff(outcome_range)

    # TODO: confirm that max(y_rescaled) <= 1 and min(y_rescaled) >= 0

    # Remove old version just to be safe.
    data[[outcome]] = NULL


    outcome_orig = outcome

    # Use the transformed version of the outcome.
    data[["Y_"]] = y_rescaled
    outcome = "Y_"

    if (verbose && family == "gaussian") {
      cat("Rescaled outcome:\n")
      print(summary(data[[outcome]]))
    }

  } else {
    outcome_orig = NULL
  }

  # Setup cross-validation folds, stratified on the outcome.
  folds = rsample::vfold_cv(data, v = folds_cvtmle, strata = outcome)

  ##################
  # Loop over folds, analyzing training sets and then applying to test.

  fold_results =
    purrr::map(folds$splits, analyze_folds,
               outcome = outcome,
               exposures = exposures,
               family = family,
               quantiles_mixtures = quantiles_mixtures,
               quantiles_exposures = quantiles_exposures,
               folds_sl = folds_sl,
               estimator_outcome = estimator_outcome,
               estimator_propensity = estimator_propensity,
               mixture_fn = mixture_fn,
               cluster_exposures = cluster_exposures,
               verbose = verbose)


  ##################
  # Combine test set results (CV-TMLE etc.)

  combined_results =
    combine_test_results(fold_results, family = family,
                         quantiles_mixtures = quantiles_mixtures,
                         verbose = verbose)

  ##################
  # Analyze results for exposure groups (FDR adjustment, etc.)
  group_results =
    analyze_exposure_groups(combined_results,
                            data = data,
                            outcome = outcome,
                            exposures = exposures,
                            mixture_fn = mixture_fn,
                            refit_mixtures = refit_mixtures,
                            verbose = verbose)

  ##################
  # Compile everything we want to return.
  results =
    list(folds = fold_results,
         combined = combined_results,
         groups = group_results,
         fold_sets = folds,
         outcome = outcome,
         outcome_orig = outcome_orig,
         outcome_bounds = q_bounds,
         outcome_rescaled = data[[outcome]],
         exposures = exposures,
         mixture_fn = mixture_fn,
         quantiles_exposures = quantiles_exposures,
         quantiles_mixtures = quantiles_mixtures,
         rescaled = needs_rescale,
         family = family)

  class(results) = "tlmixture"

  if (verbose) {
    cat("tlmixture analysis complete.\n")
  }

  return(results)
}
