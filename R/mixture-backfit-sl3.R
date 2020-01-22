#' Create exposure weights via SL with backfitting
#'
#' @param data tbd
#' @param outcome outcome column name
#' @param exposures tbd
#' @param estimator_mixture sl3::Lrnr_sl with offset-supporting learners
#' @param estimator_confounders sl3::Lrnr_sl with offset-supporting learners
#' @param exposure_groups List of all exposure groups
#' @param max_iterations tbd
#' @param tolerance tbd
#' @param quantiles tbd
#' @param family Not used for this mixture estimator.
#' @param verbose tbd
#' @param ... tbd
#'
#' @importFrom stats as.formula as.formula binomial coef lm
#'
#' @export
mixture_backfit_sl3 =
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
           max_iterations = 10L,
           tolerance = 0.00001,
           quantiles = NULL, family = "continuous", verbose = FALSE,
           ...) {

  if (verbose) {
    cat("Create mixture via backfit SL3.\n")
  }
    
  if (is.null(estimator_confounders)) {
    estimator_confounders = estimator_mixture$clone()
  }
    
  if (family == "binary") {
    family = "binomial"
  }
    
  confounders = names(data)[!names(data) %in% c(outcome, exposures)]
 
  # Initialize mixture offset
  data$offset_mixture = 0
  
  # Track coefficients over iterations (not used).
  coefs_mixture = data.frame(matrix(nrow = max_iterations,
                    ncol = length(exposures) + 1))
  colnames(coefs_mixture) = c("Intercept", exposures)
  
  #browser()

  for (iteration in seq(max_iterations)) {
    
    # Setup mixture estimation task
    task_mixture = sl3::make_sl3_Task(data = data,
                                      covariates = exposures,
                                      outcome = outcome,
                                      offset = "offset_mixture",
                                      # Continuous or binomial
                                      outcome_type = family)
    
    # TODO: incorporate the number of SL CV folds specified somehow.
    # TODO: use metalearner SuperLearner::method.CC_LS
    
    fit_mix = estimator_mixture$train(task_mixture)

    if (verbose) {
      cat("Mixture regression:\n")
      print(fit_mix)
    }
    
    # Predicted mixture value
    f_a = fit_mix$predict()

    # Calculate correction
    correction = mean(f_a)

    # Residualize
    # y_star = data[[outcome]] - (f_a - correction)
    data$offset_confounders = f_a - correction
    
    # Setup confounder adjustment task
    # TODO: don't use an offset if family = gaussian, to support more learners.
    # Just use residual instead.
    # TODO: use folds_confounders hyperparam.
    task_confounders = sl3::make_sl3_Task(data = data,
                                      covariates = confounders,
                                      outcome = outcome,
                                      offset = "offset_confounders",
                                      # Continuous or binomial
                                      outcome_type = family)
    
    fit_confounders = estimator_confounders$train(task_confounders)

    if (verbose) {
      cat("Confounder adjustment regression:\n")
      print(fit_confounders)
    }

    g_w = fit_confounders$predict()

    # Residualize
    #y_star = data[[outcome]] - g_w
    # Update offset
    data$offset_mixture = g_w

    # Check for convergence and stop early
    # Sum of the absolute change in the coefficients
    if (iteration > 1) {

      # Track optimization progress: max absolute change.
      # Could instead track the mean(abs(change))
      # Normalize by standard deviation so that tolerance is independent of outcome scale.

      # SD will be 0 if the outcome mean was given all weight in the SL estimator.
      if (sd(f_a) > 0) {
        change_f_a = max(abs(f_a - old_f_a)) / sd(f_a)
      } else {
        # Set to 0 so that it won't stop convergence for this iteration.
        change_f_a = 0
      }

      # SD will be 0 if the outcome mean was given all weight in the SL estimator.
      if (sd(g_w) > 0) {
        change_g_w = max(abs(g_w - old_g_w)) / sd(g_w)
      } else {
        #change_g_w = max(abs(g_w - old_g_w))
        # Set to 0 so that it won't stop convergence for this iteration.
        change_g_w = 0
      }

      # TODO: track these two values over time.
      if (verbose) {
        cat(paste0(iteration, "."), "max change f_a:", change_f_a, "max change g_w:", change_g_w, "\n")
      }

      # Stop iteration if we're below tolerance
      if (max(change_f_a, change_g_w) < tolerance) {

        # Truncate coefficients df so that we don't have a bunch of empty rows.
        #coefs_mixture = coefs_mixture[1:iteration, ]

        break
      }

    } else {
      if (verbose) {
        cat("1. Correction:", correction, "\n")
      }
    }

    # Save these values for the next iteration
    old_f_a = f_a
    old_g_w = g_w
  }
  
  
  # Check if our predicted mixture is constant, meaning that SL.mean has 100% weight.
  if (sd(f_a) == 0) {
    browser()
  }

  # TODO: consider checking if sd(f_a) = 0, in which case SL.mean is the estimator
  # and we won't have any variance in the histogram.

  results = list(reg_mixture = fit_mix,
                 reg_adjust = fit_confounders,
                 coefs_mixture = coefs_mixture,
                 outcome = outcome,
                 exposures = exposures,
                 iterations = iteration,
                 max_iterations = max_iterations,
                 converged = iteration < max_iterations,
                 family = family,
                 # Set weights to 0 for now.
                 weights = rep(0, length(exposures)))

  class(results) = "mixture_backfit_sl3"

  return(results)
}

# TODO: document with roxygen
#' tbd
#'
#' @param object tbd
#' @param data tbd
#' @param ... tbd
predict.mixture_backfit_sl3 = function(object, data, ...) {
  
  # We need to create a blank offset to avoid an sl3 error.
  data$offset_mixture = 0
  
  # Setup mixture estimation task
  task_mixture = sl3::make_sl3_Task(data = data,
                                    covariates = object$exposures,
                                    offset = "offset_mixture",
                                    # Continuous or binomial
                                    outcome_type = object$family)
  
  
  preds = object$reg_mixture$predict(task_mixture)
  
  # Clear this variable in case we're modifying a data.table by reference.
  data$offset_mixture = NULL
  
  return(preds)
}
