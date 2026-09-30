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
#' @param mixture_fn Function that estimates the mixture on training data, e.g.
#'   mixture_glm, mixture_pls, mixture_sl, mixture_score_glm, or
#'   mixture_score_sl3. To change a mixture function's settings, wrap it and
#'   forward \code{...} (see Examples); forwarding \code{...} also passes
#'   \code{id} through when clustering, so the mixture's internal
#'   cross-validation keeps clusters together.
#' @param refit_mixtures After CV-TMEL, refit mixture functions to full dataset.
#' @param verbose If TRUE, display more detailed info during execution.
#' @param id Optional cluster identifier(s) for cluster-robust standard errors,
#'   for data in which observations are correlated within groups (e.g. siblings
#'   in a family, children in a site). The simplest form is the name(s) of
#'   columns in \code{data}, e.g. \code{id = "family"}, or
#'   \code{id = c("school", "neighborhood")} for multiway clustering (Cameron,
#'   Gelbach & Miller 2011); those columns are used as clusters and are not
#'   used as confounders. Alternatively pass the ids directly: a vector with one
#'   value per row of \code{data}, or a list or data frame of such vectors.
#'   When given: CV-TMLE folds keep clusters together
#'   (\code{rsample::group_vfold_cv()}, stratified on each cluster's mean
#'   outcome); the
#'   outcome and propensity SuperLearners group their internal CV folds by
#'   cluster; the ids are passed to \code{mixture_fn} as \code{id} when it
#'   accepts an \code{id} or \code{...} argument, so the built-in mixture
#'   functions with internal cross-validation (\code{mixture_sl},
#'   \code{mixture_score_glm}, \code{mixture_score_sl3},
#'   \code{mixture_backfit_sl}, \code{mixture_backfit_sl3}) keep clusters
#'   together too, including in the full-data refit; influence-curve variances sum within clusters;
#'   and confidence intervals and p-values use a t distribution with (fewest
#'   clusters - 1) degrees of freedom. A variable nested in another (family
#'   within site) is covered by clustering on the outer one and is dropped with
#'   a message. NULL (the default) assumes independent observations.
#' @param folds_by With several crossed cluster variables, the name of the
#'   cluster variable to form CV-TMLE folds on. By default folds keep every
#'   cluster variable intact when possible, and otherwise use the coarsest
#'   variable with at least \code{folds_cvtmle} clusters.
#'
#' @examples
#' \donttest{
#' # Simulated data: 150 families of 1-3 children, with a shared family effect.
#' set.seed(1)
#' family = rep(1:150, times = sample(1:3, 150, replace = TRUE))
#' n = length(family)
#' data = data.frame(w1 = rnorm(n), w2 = rnorm(n),
#'                   a1 = rnorm(n), a2 = rnorm(n), a3 = rnorm(n),
#'                   family = family)
#' data$y = data$a1 + 0.5 * data$a2 + data$w1 + rnorm(150)[family] + rnorm(n)
#'
#' # Cluster-robust standard errors: name the cluster column in data. It is
#' # used for clustering only, never as a confounder.
#' result = tlmixture(data, outcome = "y", exposures = c("a1", "a2", "a3"),
#'                    estimator_outcome = c("SL.mean", "SL.glm"),
#'                    id = "family")
#' result$combined$groups_df
#' result$cluster_df   # number of families - 1
#'
#' # Nested clusters (families within sites): pass both; clustering on the
#' # outer variable covers the inner one.
#' # tlmixture(data, ..., id = c("family", "site"))
#' }
#'
#' \dontrun{
#' # A customized mixture function: wrap it and forward the dots. Forwarding
#' # the dots also passes the cluster ids through, so the sl3 learners'
#' # internal cross-validation keeps families together.
#' library(sl3)
#' my_mixture = function(...) {
#'   mixture_score_sl3(
#'     ...,
#'     estimator_mixture = Lrnr_sl$new(
#'       learners = Stack$new(Lrnr_glm$new(), Lrnr_xgboost$new(nrounds = 100)),
#'       metalearner = Lrnr_nnls$new(convex = TRUE)),
#'     estimator_confounders = Lrnr_sl$new(
#'       learners = Stack$new(Lrnr_glm$new(), Lrnr_ranger$new(num.trees = 500)),
#'       metalearner = Lrnr_nnls$new(convex = TRUE)))
#' }
#' result = tlmixture(data, outcome = "y", exposures = c("a1", "a2", "a3"),
#'                    estimator_outcome = c("SL.mean", "SL.ranger"),
#'                    mixture_fn = my_mixture,
#'                    id = "family")
#' }
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
           verbose = FALSE,
           id = NULL,
           folds_by = NULL
           ) {

  ##################
  # Initial setup

  # An id that was supplied but evaluates to NULL is almost always a typo such
  # as id = data$famly (a missing column gives NULL), which would otherwise
  # silently turn clustering off.
  if (!missing(id) && is.null(id)) {
    warning("id was supplied but is NULL (for example data$column with a ",
            "misspelled column name), so standard errors are NOT clustered.",
            call. = FALSE)
  }

  # id can name columns of data. Those columns become the cluster ids and are
  # removed from data, because every column other than the outcome and the
  # exposures is otherwise used as a confounder. A character vector as long as
  # the data is taken as the ids themselves, not as column names.
  if (is.character(id) && length(id) < nrow(data)) {
    missing_cols = setdiff(id, names(data))
    if (length(missing_cols) > 0L) {
      stop("id names column(s) not found in data: ",
           paste(missing_cols, collapse = ", "), ".")
    }
    id_cols = unique(id)
    if (any(id_cols %in% c(outcome, unlist(exposures)))) {
      stop("id cannot name the outcome or an exposure.")
    }
    id = data[, id_cols, drop = FALSE]
    data = data[, setdiff(names(data), id_cols), drop = FALSE]
    if (verbose) {
      cat("Clustering on column(s)", paste(id_cols, collapse = ", "),
          "of data; they are not used as confounders.\n")
    }
  }

  # Validate the cluster ids once, up front. NULL means no clustering, which
  # includes an id in which every observation is its own cluster.
  cluster_id = check_cluster_id(id, nrow(data), verbose = TRUE)
  df_cluster = cluster_df(cluster_id)
  if (!is.null(cluster_id) && verbose) {
    cat("Cluster-robust standard errors using:",
        paste0(names(cluster_id), " (",
               vapply(cluster_id, function(x) length(unique(x)), numeric(1)),
               " clusters)", collapse = ", "),
        "\nt critical values with", df_cluster, "degrees of freedom.\n")
  }

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
  # With clustering, whole clusters are assigned to folds instead, stratified
  # on each cluster's mean outcome (group_vfold_cv() needs strata that are
  # constant within a cluster, so the outcome itself cannot be used).
  fold_groups = cluster_fold_groups(cluster_id, folds_cvtmle, folds_by)
  if (is.null(fold_groups)) {
    folds = rsample::vfold_cv(data, v = folds_cvtmle, strata = outcome)
  } else {
    folds = cluster_vfold_cv(data, fold_groups, v = folds_cvtmle,
                             strata_values = data[[outcome]])
  }

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
               fold_groups = fold_groups,
               verbose = verbose)


  ##################
  # Combine test set results (CV-TMLE etc.)

  combined_results =
    combine_test_results(fold_results, family = family,
                         quantiles_mixtures = quantiles_mixtures,
                         verbose = verbose,
                         cluster_id = cluster_id,
                         df = df_cluster)

  ##################
  # Analyze results for exposure groups (FDR adjustment, etc.)
  group_results =
    analyze_exposure_groups(combined_results,
                            data = data,
                            outcome = outcome,
                            exposures = exposures,
                            mixture_fn = mixture_fn,
                            refit_mixtures = refit_mixtures,
                            id = fold_groups,
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
         family = family,
         # Cluster ids actually used (after dropping nested and singleton
         # variables; NULL without clustering) and the degrees of freedom
         # behind the t critical values.
         cluster_id = cluster_id,
         cluster_df = df_cluster)

  class(results) = "tlmixture"

  if (verbose) {
    cat("tlmixture analysis complete.\n")
  }

  return(results)
}
