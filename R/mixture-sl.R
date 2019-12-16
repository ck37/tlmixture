#' Create exposure weights
#'
#' @param data tbd
#' @param outcome outcome column name
#' @param exposures tbd
#' @param quantiles tbd
#' @param family tbd
#' @param sl_library tbd
#' @param verbose tbd
#'
#' @importFrom stats as.formula as.formula binomial coef glm lm
#'
#' @export
mixture_sl =
  function(data, outcome, exposures, quantiles, family = "binomial", verbose = FALSE,
           sl_library = c("SL.mean", "SL.glm")) {
  if (verbose) {
    cat("Create mixture via SuperLearner.\n")
  }

  data_x = data[, !names(data) %in% outcome, drop = FALSE]

  sl = SuperLearner::SuperLearner(Y = data[[outcome]], X = data_x,
                                  family = family,
                                  SL.library = sl_library,
                                  cvControl = list(V = 2),
                                  verbose = FALSE)

  results = list(sl = sl, weights = NULL)
  class(results) = "mixture_sl"

  return(results)
}

# TODO: document with roxygen
#' tbd
#'
#' @param object bd
#' @param data tbd
#' @param ... tbd
predict.mixture_sl = function(object, data, ...) {
  preds = predict(object$sl, data, onlySL = TRUE)$pred[, 1]
  return(preds)
}
