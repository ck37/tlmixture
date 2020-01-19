#' Create exposure weights via lm with backfitting
#'
#' @param data tbd
#' @param outcome outcome column name
#' @param exposures tbd
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
mixture_backfit_lm =
  function(data, outcome, exposures,
           exposure_groups = NULL,
           max_iterations = 50L,
           tolerance = 0.00001,
           quantiles = NULL, family = gaussian(), verbose = FALSE,
           ...) {
    
  if (verbose) {
    cat("Create mixture via backfit LM.\n")
  }
    
  # Setup predictor dataframes so that we don't have to do it repeatedly.
  df_confounders = data[, !names(data) %in% c(outcome, exposures), drop = FALSE]
  df_exposures = data[, exposures, drop = FALSE]
  
  # Track coefficients over iterations.
  coefs_mixture = data.frame(matrix(nrow = max_iterations,
                    ncol = length(exposures) + 1))
  colnames(coefs_mixture) = c("Intercept", exposures)

  # Initialize  Y*
  y_star = data[[outcome]]

  for (iteration in seq(max_iterations)) {
    
    # Estimate mixture function
    # (Could use GLM with offset)
    reg_mixture = lm(y_star ~ ., data = df_exposures)
    
    coefs_mixture[iteration, ] = coef(reg_mixture)
    
    # Predicted mixture value
    f_a = reg_mixture$fitted.values
    
    # Calculate correction
    correction = mean(f_a)
    
    # Residualize
    y_star = data[[outcome]] - (f_a - correction)
    
    # Estimate adjustment function 
    reg_adjust = lm(y_star ~ ., data = df_confounders)
    
    g_w = reg_adjust$fitted.values
    
    # Residualize
    y_star = data[[outcome]] - g_w
    
    # Check for convergence and stop early
    # Sum of the absolute change in the coefficients
    if (iteration > 1) {
      
      # Track optimization progress.
      current_coefs = coefs_mixture[iteration, ]
      prior_coefs = coefs_mixture[iteration - 1, ]
      coef_change = max(abs(current_coefs - prior_coefs))
      
      # Stop iteration if we're below tolerance
      if (coef_change < tolerance) {
        
        # Truncate coefficients df so that we don't have a bunch of empty rows.
        coefs_mixture = coefs_mixture[1:iteration, ]
        break
      }
    }
    
  }
  
  results = list(reg_mixture = reg_mixture,
                 reg_adjust = reg_adjust,
                 coefs_mixture = coefs_mixture,
                 outcome = outcome,
                 exposures = exposures,
                 iterations = iteration,
                 max_iterations = max_iterations,
                 converged = iteration < max_iterations,
                 weights = coef(reg_mixture)[exposures])

  class(results) = "mixture_backfit_lm"

  return(results)
}

# TODO: document with roxygen
#' tbd
#'
#' @param object tbd
#' @param data tbd
#' @param ... tbd
predict.mixture_backfit_lm = function(object, data, ...) {
  preds = predict(object$reg_mixture, data[, object$exposures])
  return(preds)
}
