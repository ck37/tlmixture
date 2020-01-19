#' Create exposure weights via lm
#'
#' @param data tbd
#' @param outcome outcome column name
#' @param exposures tbd
#' @param exposure_groups List of all exposure groups
#' @param quantiles tbd
#' @param family Not used for this mixture estimator.
#' @param verbose tbd
#' @param ... tbd
#'
#' @importFrom stats as.formula as.formula binomial coef lm
#'
#' @export
mixture_lm =
  function(data, outcome, exposures,
           exposure_groups,
           quantiles, family = gaussian(), verbose = FALSE,
           ...) {
    
  if (verbose) {
    cat("Create mixture via LM.\n")
  }

  formula = as.formula(paste(outcome, "~ ."))

  # TODO: standardize variables? Or do outside of the mixture estimation? Or option?

  estimator =
    stats::lm(formula, data = data)

  results = list(estimator = estimator,
                 outcome = outcome,
                 exposures = exposures,
                 weights = coef(estimator)[exposures])

  class(results) = "mixture_lm"

  return(results)
}

# TODO: document with roxygen
#' tbd
#'
#' @param object tbd
#' @param data tbd
#' @param ... tbd
predict.mixture_lm = function(object, data, ...) {
  # Extract just the exposures and apply the coefficients, then the link family.
  # TODO: should we generate mixture on the logit scale or the probability scale?
  weights = coef(object$estimator)[object$exposures]
  # Replace any NA weights with 0 (e.g. due to collinearity)
  weights[is.na(weights)] = 0

  preds = as.vector(weights %*% t(data[, object$exposures, drop = FALSE]))
  return(preds)
}
