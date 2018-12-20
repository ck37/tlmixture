#' Detect outcome family based on the outcome data vector
#' @export
outcome_family = function(data, outcome_field) {

  unique_vals = unique(data[[outcome]])
  if (length(unique_vals) == 2L &&
      length(setdiff(unique(data[[outcome]]), c(0, 1))) == 0L) {
    family = "binomial"
  } else {
    family = "gaussian"
  }

  return(family)

}
