#' Prints tlmixture result object
#' 
#' Prints tlmixture result object
#'
#' @param x tlmixture result object
#' @param ... Not used currently.
#' @export
print.tlmixture = function(x, ...) {
  cat("Combined results:\n")
  print(x$combined$results)
}