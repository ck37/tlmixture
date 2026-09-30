#' Create exposure weights
#'
#' @param data tbd
#' @param outcome outcome column name
#' @param exposures tbd
#' @param exposure_groups tbd
#' @param quantiles tbd
#' @param family tbd
#' @param sl_library tbd
#' @param id Optional cluster identifier, one value per row of \code{data}
#'   (tlmixture() passes it when its \code{id} argument is used). When given,
#'   SuperLearner's internal cross-validation keeps clusters together. NULL
#'   (default) leaves the estimator unchanged.
#' @param verbose tbd
#'
#' @importFrom stats as.formula as.formula binomial coef glm lm
#'
#' @export
mixture_sl =
  function(data, outcome,
           exposures,
           exposure_groups,
           quantiles, family = "binomial", verbose = FALSE,
           sl_library = c("SL.mean", "SL.glm"),
           id = NULL) {
  if (verbose) {
    cat("Create mixture via SuperLearner.\n")
  }

  data_x = data[, !names(data) %in% outcome, drop = FALSE]

  # With cluster ids, keep clusters within SuperLearner's internal folds.
  sl_args = sl_cluster_args(id, 2)
  sl = SuperLearner::SuperLearner(Y = data[[outcome]], X = data_x,
                                  family = family,
                                  SL.library = sl_library,
                                  id = sl_args$id,
                                  cvControl = list(V = sl_args$V),
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
