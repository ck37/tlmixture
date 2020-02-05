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
#' @param adjust_other_exposures Whether or not to include other exposure groups as confounders to adjust for. Default TRUE.
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
           max_iterations = 5L,
           tolerance = 0.00001,
           adjust_other_exposures = TRUE,
           force_offset = TRUE,
           quantiles = NULL, family = "continuous", verbose = FALSE,
           debug = FALSE,
           ...) {

  if (verbose) {
    cat("Create mixture via backfit SL3.\n")
  }

  if (is.null(estimator_confounders)) {
    estimator_confounders = estimator_mixture$clone()
  }

  if (family == "binary") {
    family = "binomial"
    # For binomial we want the default offset to be a 50% probability,
    # so that it becomes a 0 when tranformed to the linear logit scale.
  }

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

  for (iteration in seq(max_iterations)) {

    if (debug) {
      cat("Mixture offset sd:", sd(data$offset_mixture), "\n")
      cat("Mixture offset summary:\n")
      print(summary(data$offset_mixture))
    }

    # Setup mixture estimation task
    if (family == "binomial" || force_offset) {
      task_mixture = sl3::make_sl3_Task(data = data,
                                        covariates = exposures,
                                        outcome = outcome,
                                        offset = "offset_mixture",
                                        outcome_type = family)
    } else {
      # Create the residualized outcome.
      data$`_outcome_ym` = data[[outcome]] - data$offset_mixture
      task_mixture = sl3::make_sl3_Task(data = data,
                                        covariates = exposures,
                                        outcome = "_outcome_ym",
                                        outcome_type = family)
    }

    # TODO: incorporate the number of SL CV folds specified somehow.
    # TODO: use metalearner SuperLearner::method.CC_LS

    fit_mix = estimator_mixture$train(task_mixture)

    if (debug) {
      cat("Mixture regression:\n")
      print(fit_mix)
    }

    if (family == "binomial" || force_offset) {
      # Update offset to be 0.
      data$offset_mixture = offset_null

      task_mixture = sl3::make_sl3_Task(data = data,
                                        covariates = exposures,
                                        outcome = outcome,
                                        offset = "offset_mixture",
                                        outcome_type = family)
    }
    # We don't need a new task if family is gaussian.

    # Predicted mixture value
    f_a = fit_mix$predict(task_mixture)

    if (debug) {
      print(qplot(f_a) + ggtitle(paste("mixture f_a iteration", iteration)) + theme_minimal())
    }

    # Calculate correction
    if (family == "binomial") {
      # If binomial calculate on the logit scale to ensure that we remain in bounds.
      correction = mean(qlogis(f_a))
    } else {
      correction = mean(f_a)
    }

    if (debug) {
      cat("Correction:", correction, "\n")
    }

    # Check if our predicted mixture is constant, meaning that SL.mean has 100% weight.
    if (sd(f_a) == 0) {
      cat("Error: mixture prediction has no variation.\n")
      #browser()
    }

    # Residualize
    # y_star = data[[outcome]] - (f_a - correction)
    if (family == "binomial") {
      # Correct predictions on the logit scale, then convert back to probability scale.
      data$offset_confounders = plogis(qlogis(f_a) - correction)
      # TODO: confirm that offset is within (0, 1)
      if (debug) {
        cat("Summary of confounder offsets:\n")
        print(summary(data$offset_confounders))
      }
    } else {
      data$offset_confounders = f_a - correction
    }
    #if (family == "binomial") {
    #  # sl3 expects offset to be on the probability scale rather than logit scale.
    #  data$offset_confounders = qlogis(data$offset_confounders)
    #}

    # Setup confounder adjustment task
    # TODO: don't use an offset if family = gaussian, to support more learners.
    # Just use residual instead.
    # TODO: use folds_confounders hyperparam.
    if (family == "binomial" || force_offset) {
      task_confounders =
        sl3::make_sl3_Task(data = data,
                           covariates = confounders,
                           outcome = outcome,
                           offset = "offset_confounders",
                           outcome_type = family)
    } else {
      # Create the residualized outcome.
      data$`_outcome_yc` = data[[outcome]] - data$offset_confounders
      task_confounders =
        sl3::make_sl3_Task(data = data,
                           covariates = confounders,
                           outcome = "_outcome_yc",
                           outcome_type = family)
    }

    fit_confounders = estimator_confounders$train(task_confounders)

    if (debug) {
      cat("Confounder adjustment regression:\n")
      print(fit_confounders)
    }

    # Update offset to be 0.
    #data$offset_confounders = 0
    data$offset_confounders = offset_null

    if (family == "binomial" || force_offset) {
      task_confounders =
        sl3::make_sl3_Task(data = data,
                           covariates = confounders,
                           outcome = outcome,
                           offset = "offset_confounders",
                           outcome_type = family)
    }
    # For family = gaussian we don't need to update the task.

    g_w = fit_confounders$predict(task_confounders)

    if (debug) {
      print(qplot(g_w) + ggtitle(paste("g_w iteration", iteration)) + theme_minimal())
    }

    # Residualize
    #y_star = data[[outcome]] - g_w
    # Update offset
    # sl3 expects offset to be on the probability scale rather than logit scale.
    data$offset_mixture = g_w

    if (debug) {
      cat("Correlation of g_w and f_a:", cor(f_a, g_w), "\n")
    }

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


  # TODO: consider checking if sd(f_a) = 0, in which case SL.mean is the estimator
  # and we won't have any variance in the histogram.

  results = list(reg_mixture = fit_mix,
                 reg_adjust = fit_confounders,
                 coefs_mixture = NULL,
                 outcome = outcome,
                 exposures = exposures,
                 confounders = confounders,
                 iterations = iteration,
                 force_offset = force_offset,
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
#' @param type "mixture" (default) or "confounders"
#' @param ... tbd
#'
#' @export
predict.mixture_backfit_sl3 = function(object, data, type = "mixture", ...) {

  if (object$family == "binomial") {
    # This becomes a 0 when transformed to the logit scale.
    offset_null = 0.5
  } else {
    offset_null = 0
  }

  # We need to create a blank outcome to avoid an sl3 error.
  data$`_outcome` = 0

  if (type == "mixture") {
    if (object$family == "binomial" || object$force_offset) {
      # We need to create a blank offset to avoid an sl3 error.
      #data$offset_mixture = 0
      data$offset_mixture = offset_null

      # Setup mixture estimation task
      task_mixture = sl3::make_sl3_Task(data = data,
                                        covariates = object$exposures,
                                        offset = "offset_mixture",
                                        outcome = "_outcome",
                                        # Continuous or binomial
                                        outcome_type = object$family)
    } else {
      # Setup mixture estimation task
      task_mixture = sl3::make_sl3_Task(data = data,
                                        covariates = object$exposures,
                                        outcome = "_outcome",
                                        # Continuous or binomial
                                        outcome_type = object$family)
    }

    preds = object$reg_mixture$predict(task_mixture)

    if (object$family == "binomial" || object$force_offset) {
      # Clear this variable in case we're modifying a data.table by reference.
      data$offset_mixture = NULL
    }
  } else {
    # We need to create a blank offset to avoid an sl3 error.
    # TODO: or should this be set to the predicted mixture value?
    if (object$family == "binomial" || object$force_offset) {
      #data$offset_confounders = 0
      data$offset_confounders = offset_null

      # Setup confounder estimation task
      task_confounders =
        sl3::make_sl3_Task(data = data,
                           covariates = object$confounders,
                           outcome = "_outcome",
                           offset = "offset_confounders",
                           outcome_type = object$family)
    } else {
      # Setup confounder estimation task
      task_confounders =
        sl3::make_sl3_Task(data = data,
                           covariates = object$confounders,
                           outcome = "_outcome",
                           outcome_type = object$family)
    }

    preds = object$reg_adjust$predict(task_confounders)

    if (object$family == "binomial" || object$force_offset) {
      # Clear this variable in case we're modifying a data.table by reference.
      data$offset_confounders = NULL
    }
  }

  # Clean up blank outcome
  data$`_outcome` = NULL

  return(preds)
}
