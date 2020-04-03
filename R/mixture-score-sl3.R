#' Create score-based mixture via sl3
#'
#' @param data tbd
#' @param outcome outcome column name
#' @param exposures tbd
#' @param estimator_mixture sl3::Lrnr_sl with offset-supporting learners
#' @param estimator_confounders sl3::Lrnr_sl with offset-supporting learners
#' @param folds_mixture SL folds to use during mixture estimation; default 5. (Not yet implemented).
#' @param folds_confounders SL folds to use during confounder estimation; default 5. (Not yet implemented).
#' @param exposure_groups List of all exposure groups
#' @param max_iterations tbd
#' @param tolerance tbd
#' @param adjust_other_exposures Whether or not to include other exposure groups as confounders to adjust for. Default TRUE.
#' @param force_offset Whether to use offset style estimation for gaussian outcomes; default True.
#' @param quantiles tbd
#' @param family Not used for this mixture estimator.
#' @param verbose tbd
#' @param debug If True show very detailed output.
#' @param ... tbd
#'
#' @importFrom stats as.formula as.formula binomial coef lm
#'
#' @export
mixture_score_sl3 =
  function(data, outcome, exposures,
           estimator_mixture =
             sl3::Lrnr_sl$new(learners =
                                sl3::make_learner(sl3::Stack,
                                                  sl3::make_learner(sl3::Lrnr_mean),
                                                  sl3::make_learner(sl3::Lrnr_glm)),
                         metalearner = sl3::make_learner(sl3::Lrnr_nnls)),
           # This will be cloned from estimator_mixture if it is not specified.
           estimator_confounders = NULL,
           folds_mixture = 5L,
           folds_confounders = 5L,
           exposure_groups = NULL,
           max_iterations = 5L,
           tolerance = 0.00001,
           adjust_other_exposures = TRUE,
           force_offset = TRUE,
           quantiles = NULL, family = "continuous", verbose = FALSE,
           debug = FALSE,
           ...) {

  if (verbose) {
    cat("Create score-based mixture via sl3.\n")
  }

  if (is.null(estimator_confounders)) {
    estimator_confounders = estimator_mixture$clone()
  }

  if (family == "binary") {
    family = "binomial"
    # For binomial we want the default offset to be a 50% probability,
    # so that it becomes a 0 when tranformed to the linear logit scale.
  }
    
  # Fix family for glm()
  #if (family == "continuous") {
  #  family = "gaussian"
  #}

  # This is a separate if statement because we want family to be standardized to
  # either "binomial" or "continuous" at this point.
  if (family == "binomial") {
    offset_null = 0.5
  } else {
    offset_null = 0
  }

  confounders = names(data)[!names(data) %in% c(outcome, exposures)]

  # Remove other exposure variables from the confounder list.
  if (!adjust_other_exposures) {
    confounders = setdiff(confounders, unlist(exposure_groups))
  }

  if (debug) {
    cat("Family:", family, "\n")
    cat("Exposures:", exposures, "\n")
    cat("Confounders:", confounders, "\n")
  }

  # Only use offset for binomial outcomes.
  #if (family == "binomial") {
  # Initialize mixture offset
  # sl3 expects offset to be on the probability scale rather than logit scale.
  # For continuous variable, this will be subtracted from the outcome before prediction.
  data$offset_mixture = offset_null
  #}
  
  # 1. Cross-validated exposure-naive outcome regression E[ Y | W ]
  # Use estimator_confounders
  
  # Create cross-validation folds.
  cv_folds = 5L
  folds = sample(rep(seq(cv_folds), length.out = nrow(data)))
 
  # Create a standalone vector for the exposure-naive outcome prediction. 
  # This could alternatively be a column in the dataframe.
  exposure_naive_outcome = vector(mode = "numeric", nrow(data))
  
  # Loop over training data and predict on test.
  for (test_fold in seq(cv_folds)) {
    if (verbose) {
      cat(paste("Exposure-naive outcome regression:", test_fold, "of", cv_folds, "\n"))
    }
    train_df = data[folds != test_fold, ]
    test_df = data[folds == test_fold, ]
    
    # Fit E[ Y | W] on training data.
    task_confounders =
      sl3::make_sl3_Task(data = train_df,
                         covariates = confounders,
                         outcome = outcome,
                         #offset = "offset_confounders",
                         outcome_type = family)
    
    fit_confounders = estimator_confounders$train(task_confounders)
    
    if (verbose) {
      print(fit_confounders)
    }
    
    # Apply to test
    task_confounders_pred =
      sl3::make_sl3_Task(data = test_df,
                         covariates = confounders,
                         outcome = outcome,
                        # offset = "offset_confounders",
                         outcome_type = family)
   
    pred_test = fit_confounders$predict(task_confounders_pred) 
    #pred_test = predict(train_model, test_df, type = "response")
    
    # Save the result.
    exposure_naive_outcome[folds == test_fold] = pred_test
    
  }
  
  # 2. Create meta-level dataset
  meta_df = cbind(data[, c(outcome, exposures, confounders)],
                  `_offset` = exposure_naive_outcome)
  
  # Create interaction-term covariates
  
  # 3. Mixture regression: E[ Y | A', W, E[ Y | W]
  task_mixture = sl3::make_sl3_Task(data = meta_df,
                                    covariates = c(exposures, confounders),
                                    offset = "_offset",
                                    outcome = outcome,
                                    outcome_type = family)
  # Include two-way interactions of A and W, including A:A
  #interaction_grid = expand.grid(c(exposures, confounders), exposures)
  #glm_formula =
  #  as.formula(paste0(outcome, "~ ",
  #                    paste0(interaction_grid$Var1, ":",
  #                           interaction_grid$Var2, collapse = " + ")))
  
  if (verbose) {
    cat("\nEstimating mixture function.\n")
  }
  
  fit_mix = estimator_mixture$train(task_mixture)
  
  if (verbose) {
    cat("Mixture fit:\n")
    print(fit_mix)
  }

  results = list(reg_mixture = fit_mix,
                 #reg_adjust = fit_confounders,
                 coefs_mixture = NULL,
                 outcome = outcome,
                 exposures = exposures,
                 confounders = confounders,
                 #iterations = iteration,
                 force_offset = force_offset,
                 #max_iterations = max_iterations,
                 #converged = iteration < max_iterations,
                 family = family,
                 # Set weights to 0 for now.
                 weights = rep(0, length(exposures)))

  class(results) = "mixture_score_sl3"

  return(results)
}

# TODO: document with roxygen
#' tbd
#'
#' @param object tbd
#' @param data tbd
#' @param type "mixture" (default) or "confounders"
#' @param ... tbd
#'
#' @export
predict.mixture_score_sl3 = function(object, data, ...) {
 
  if (object$family == "binomial") {
    # This becomes a 0 when transformed to the logit scale.
    offset_null = 0.5
  } else {
    offset_null = 0
  }

  # We need to create a blank outcome to avoid an sl3 error.
  data$`_outcome` = 0

  if (object$family == "binomial" || object$force_offset) {
    # We need to create a blank offset to avoid an sl3 error.
    #data$offset_mixture = 0
    data$offset_mixture = offset_null

    # Setup mixture estimation task
    task_mixture = sl3::make_sl3_Task(data = data,
                                      covariates = c(object$exposures, object$confounders),
                                      offset = "offset_mixture",
                                      outcome = "_outcome",
                                      # Continuous or binomial
                                      outcome_type = object$family)
    } else {
      # Setup mixture estimation task
      task_mixture = sl3::make_sl3_Task(data = data,
                                        covariates = c(object$exposures, object$confounders),
                                        outcome = "_outcome",
                                        # Continuous or binomial
                                        outcome_type = object$family)
    }

  preds = object$reg_mixture$predict(task_mixture)

  if (object$family == "binomial" || object$force_offset) {
    # Clear this variable in case we're modifying a data.table by reference.
    data$offset_mixture = NULL
  }

  # Clean up blank outcome
  data$`_outcome` = NULL 
  return(preds)
}
