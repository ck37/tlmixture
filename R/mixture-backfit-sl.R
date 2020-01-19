#' Create exposure weights via SL with backfitting
#'
#' @param data tbd
#' @param outcome outcome column name
#' @param exposures tbd
#' @param estimator_mixture tbd
#' @param estimator_confounders tbd
#' @param exposure_groups List of all exposure groups
#' @param max_iterations tbd
#' @param tolerance tbd
#' @param quantiles tbd
#' @param family Not used for this mixture estimator.
#' @param verbose tbd
#' @param ... tbd
#'
#' @importFrom stats as.formula as.formula binomial coef lm
#' @importFrom SuperLearner SuperLearner
#'
#' @export
mixture_backfit_sl =
  function(data, outcome, exposures,
           estimator_mixture = c("SL.mean", "SL.lm"),
           estimator_confounders = estimator_confounders,
           folds_mixture = 5L,
           folds_confounders = 5L,
           exposure_groups = NULL,
           max_iterations = 10L,
           tolerance = 0.00001,
           quantiles = NULL, family = gaussian(), verbose = FALSE,
           ...) {
    
  if (verbose) {
    cat("Create mixture via backfit SL.\n")
  }
    
  # Setup predictor dataframes so that we don't have to do it repeatedly.
  df_confounders = data[, !names(data) %in% c(outcome, exposures), drop = FALSE]
  df_exposures = data[, exposures, drop = FALSE]
  
  # Track coefficients over iterations (not used).
  coefs_mixture = data.frame(matrix(nrow = max_iterations,
                    ncol = length(exposures) + 1))
  colnames(coefs_mixture) = c("Intercept", exposures)

  # Initialize  Y*
  y_star = data[[outcome]]

  for (iteration in seq(max_iterations)) {
    
    # Estimate mixture function
    # (Could use GLM with offset)
    reg_mixture =
      SuperLearner(Y = y_star,
                   X = df_exposures,
                   family = gaussian(),
                   SL.library = estimator_mixture,
                   cvControl = list(V = folds_mixture))
    
    #coefs_mixture[iteration, ] = coef(reg_mixture)
    
    # Predicted mixture value
    f_a = reg_mixture$SL.predict
    
    # Calculate correction
    correction = mean(f_a)
    
    # Residualize
    y_star = data[[outcome]] - (f_a - correction)
    
    # Estimate adjustment function 
    reg_adjust =
      SuperLearner(Y = y_star,
                   X = df_confounders,
                   family = gaussian(),
                   SL.library = estimator_mixture,
                   cvControl = list(V = folds_confounders))
    
    g_w = reg_adjust$SL.predict
    
    # Residualize
    y_star = data[[outcome]] - g_w
    
    # Check for convergence and stop early
    # Sum of the absolute change in the coefficients
    if (iteration > 1) {
      
      # Track optimization progress: mean absolute change.
      # Could instead track the max(abs(change))
      change_f_a = max(abs(f_a - old_f_a))
      change_g_w = max(abs(g_w - old_g_w))
      
      # TODO: track these two values over time.
      cat(iteration, ". max change f_a:", change_f_a, "max change g_w:", change_g_w, "\n")
      
      # Stop iteration if we're below tolerance
      if (max(change_f_a, change_g_w) < tolerance) {
        
        # Truncate coefficients df so that we don't have a bunch of empty rows.
        #coefs_mixture = coefs_mixture[1:iteration, ]
        
        break
      }
      
    }
    
    # Save these values for the next iteration
    old_f_a = f_a
    old_g_w = g_w
  }
  
  results = list(reg_mixture = reg_mixture,
                 reg_adjust = reg_adjust,
                 coefs_mixture = coefs_mixture,
                 outcome = outcome,
                 exposures = exposures,
                 iterations = iteration,
                 max_iterations = max_iterations,
                 converged = iteration < max_iterations,
                 # Set weights to 0 for now.
                 weights = rep(0, length(exposures)))

  class(results) = "mixture_backfit_sl"

  return(results)
}

# TODO: document with roxygen
#' tbd
#'
#' @param object tbd
#' @param data tbd
#' @param ... tbd
predict.mixture_backfit_sl = function(object, data, ...) {
  preds = predict(object$reg_mixture, data[, object$exposures], onlySL = TRUE)
  return(preds$pred[, 1])
}
