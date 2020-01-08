#' Create exposure weights via glm
#'
#' @param data tbd
#' @param outcome outcome column name
#' @param exposures tbd
#' @param quantiles tbd
#' @param family tbd
#' @param verbose tbd
#' @param ... tbd
#'
#' @importFrom stats as.formula as.formula binomial coef glm lm
#'
#' @export
mixture_glm =
  function(data, outcome, exposures, quantiles, family = gaussian(), verbose = FALSE,
           ...) {
  if (verbose) {
    cat("Create mixture via GLM.\n")
  }

  data_x = data[, !names(data) %in% outcome, drop = FALSE]
  
  formula = as.formula(paste(outcome, "~ ."))
  
  estimator =
    stats::glm(formula, data = data, family = family)

  results = list(estimator = estimator,
                 outcome = outcome,
                 exposures = exposures,
                 weights = coef(estimator)[exposures])
  
  class(results) = "mixture_glm"

  return(results)
}

# TODO: document with roxygen
#' tbd
#'
#' @param object tbd
#' @param data tbd
#' @param ... tbd
predict.mixture_glm = function(object, data, ...) {
  # Extract just the exposures and apply the coefficients, then the link family.
  # TODO: should we generate mixture on the logit scale or the probability scale?
  weights = coef(object$estimator)[object$exposures]
  preds = as.vector(weights %*% t(data[, object$exposures]))
  return(preds)
}
