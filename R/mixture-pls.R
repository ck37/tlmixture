#' Create exposure weights using partial least squares
#'
#' @param data tbd
#' @param outcome outcome column name
#' @param exposures tbd
#' @param quantiles tbd
#' @param verbose tbd
#'
#' @importFrom stats as.formula as.formula binomial coef glm lm
#'
#' @export
mixture_pls =
  function(data, outcome, exposures, quantiles, verbose = FALSE) {
  if (verbose) {
    cat("Create exposure weights.\n")
  }

  # Let's use partial least squares for now.
  formula = as.formula(paste0(outcome, " ~ ."))

  # Only calculate the first component.
  # TODO: decide about incorporating adjustment variables and other exposures as well.
  # E.g. to residualize the outcome.
  reg = pls::plsr(formula,
                  data = data[, c(outcome, exposures)],
                  ncomp = 1L)

  # Extract the loadings for each component.
  # TODO: confirm we should use $loadings and not $loadings.weights (unclear what the difference is).
  loadings_mat = matrix(reg$loadings, nrow = nrow(reg$loadings))

  # Extract the first component.
  weights = loadings_mat[, 1]

  # Determine which is more common: positive or negative loadings.
  avg_positive = mean(weights >= 0, na.rm = TRUE)

  # If we have more negative than positive loadings, reverse the sign.
  # (This would happen with a protective mixture.)
  if (avg_positive < 0.5) {
    weights = weights * -1
  }

  # Set negative values to 0.
  weights[weights < 0] = 0

  # Could normalize to sum to 1.
  weights_norm = weights / sum(weights)

  if (verbose) {
    cat("Raw weights:", round(weights, 3),
        "Normalized:", round(weights_norm, 3), "\n")
  }

  # Return the normalized weights.
  # TODO: also consider saving the quantiling parameter, if it's implemented.
  results = list(weights = weights_norm, exposures = exposures)

  # Don't use this for now - could be a separate function.
  if (FALSE) {
    # We use various functions from wqs-helpers.R to estimate these weights.

    # Family is currently hard-coded to binomial but this should be detected
    # in our main function and passed through.
    family = "binomial"


    # If TRUE, exposure has a positive relationship to outcome,
    # otherwise it has a negative relationship.
    # TODO: can we determine this automatically, e.g. based on partial least squares or supervised PCA?
    b1_pos = TRUE

    # calculate parameters
    params = optim.f(data_b[, q_name, drop = FALSE],
                     # Outcome variable.
                     data_b[, y_name, drop = FALSE],

                     b1_pos,
                     b1_constr,
                     # Outcome family, either binomial or gaussian presumably.
                     family,
                     # Adjustment covariates.
                     data_b[, cov_name, drop = FALSE])
  }

  class(results) = "mixture_pls"

  return(results)
}

# TODO: document with roxygen
#' prediction for mixture_pls object
#' @param obj tbd
#' @param data tbd
#'
#' @export
predict.mixture_pls = function(obj, data) {
  # Return predicted mixture
  mixture = as.vector(as.matrix(data[, obj$exposures]) %*%
                      matrix(obj$weights, ncol = 1))

  return(mixture)
}


