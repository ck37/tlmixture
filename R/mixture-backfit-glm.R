#' Create exposure weights via glm with backfitting
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
#' @importFrom stats gaussian cor sd
#'
#' @export
mixture_backfit_glm =
  function(data, outcome, exposures,
           exposure_groups = NULL,
           max_iterations = 50L,
           tolerance = 0.00001,
           quantiles = NULL, family = gaussian(), verbose = FALSE,
           ...) {


  # Convert from tlmixture family values.
  if (family == "continuous") {
    family = "gaussian"
  }

  if (family == "binary") {
    family = "binomial"
  }

  if (verbose) {
    cat("Create mixture via backfit GLM.\n")
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

  # Initialize offset (only used if family = binomial)
  mixture_offset = rep(0, nrow(df_exposures))

  for (iteration in seq(max_iterations)) {

    # Estimate mixture function
    # (Could use GLM with offset)
    reg_mixture = glm(y_star ~ ., offset = mixture_offset,
                      family = family,
                      data = df_exposures)

    coefs_mixture[iteration, ] = coef(reg_mixture)

    # Predicted mixture value
    # For family = binomial() these are predicted probabilities.
    f_a = reg_mixture$fitted.values

    if (family == "binomial") {
      # Convert to linear log-odds scale.
      f_a = log(f_a / (1 - f_a))
    }

    # Calculate correction - should this be on the log-odds (linear) or probability scale (non-linear)?
    #correction = mean(f_a)
    correction = mean(f_a)

    # Residualize
    #y_star = data[[outcome]] - (f_a - correction)

    confounder_offset = f_a - correction
    # Convert to log-odds scale for linear prediction.
    # confounder_offset = log(confounder_offset / (1 - confounder_offset))

    # Estimate confounder adjustment function
    #reg_adjust = lm(y_star ~ ., data = df_confounders)
    reg_adjust = glm(y_star ~ ., offset = confounder_offset,
                     family = family,
                     data = df_confounders)

    g_w = reg_adjust$fitted.values

    # Residualize
    #y_star = data[[outcome]] - g_w

    # TODO: shouldn't this be on the logit scale, rather than probabilities?
    if (family == "binomial") {
      mixture_offset = log(g_w / (1 - g_w))
    }

    # Check for convergence and stop early
    # Sum of the absolute change in the coefficients
    if (iteration > 1) {

      # Track optimization progress.
      current_coefs = coefs_mixture[iteration, ]
      prior_coefs = coefs_mixture[iteration - 1, ]
      # Some coefficients may be NA due to a singular covariance matrix
      coef_change = max(abs(current_coefs - prior_coefs), na.rm = TRUE)

      # Stop iteration if we're below tolerance
      ##if (coef_change < tolerance) {

        # Truncate coefficients df so that we don't have a bunch of empty rows.
      ##  coefs_mixture = coefs_mixture[1:iteration, ]
      ## ## break
      ##}

      # Track optimization progress: max absolute change.
      # Could instead track the mean(abs(change))
      # Normalize by standard deviation
      change_f_a = max(abs(f_a - old_f_a)) / sd(f_a)
      change_g_w = max(abs(g_w - old_g_w)) / sd(g_w)

      # TODO: track these two values over time.
      if (verbose) {
        cat(paste0(iteration, "."), "max change f_a:", change_f_a, "max change g_w:", change_g_w, "\n")
      }

      # Stop iteration if we're below tolerance
      if (max(change_f_a, change_g_w) < tolerance) {

        # Truncate coefficients df so that we don't have a bunch of empty rows.
        coefs_mixture = coefs_mixture[1:iteration, ]

        break
      }

    }

    # Save these values for the next iteration
    old_f_a = f_a
    old_g_w = g_w

  }

  results = list(reg_mixture = reg_mixture,
                 reg_adjust = reg_adjust,
                 family = family,
                 coefs_mixture = coefs_mixture,
                 outcome = outcome,
                 mixture_offset = mixture_offset,
                 exposures = exposures,
                 iterations = iteration,
                 max_iterations = max_iterations,
                 converged = iteration < max_iterations,
                 weights = coef(reg_mixture)[exposures])

  class(results) = "mixture_backfit_glm"

  return(results)
}

# TODO: document with roxygen
#' tbd
#'
#' @param object tbd
#' @param data tbd
#' @param ... tbd
predict.mixture_backfit_glm = function(object, data, ...) {

  reg_obj = object$reg_mixture

  # We need to clear the offset to avoid an error unfortunately.
  # "Error in eval(object$call$offset, newdata) : object 'mixture_offset' not found"
  #if (object$family == "binomial") {
  #}

  # We are currently using an offset for both binomial and gaussian outcomes.
  reg_obj$call$offset = NULL

  preds = try({ predict(reg_obj, newdata = data[, object$exposures]) })

  if (class(preds) == "try-error") {
    browser()
  }

  return(preds)
}
