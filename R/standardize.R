#' Rescale variables, possibly excluding some columns
#'
#' A simple extension of base::scale() that skips the named columns and keeps
#' the original column order. Copied from ck37r::standardize() (same author,
#' MIT license) so that tlmixture does not need ck37r and its dependencies,
#' several of which are no longer on CRAN.
#'
#' @param x Data frame or matrix with column names.
#' @param skip_vars Names of variables not to scale.
#' @param ... Extra arguments passed through to base::scale().
#'
#' @return Data frame with the appropriate variables scaled.
#' @noRd
standardize = function(x, skip_vars = NULL, ...) {
  original_order = colnames(x)
  new_df = cbind(data.frame(scale(x[, !colnames(x) %in% skip_vars, drop = FALSE], ...)),
                 x[, colnames(x) %in% skip_vars, drop = FALSE])
  # Return with columns in their original order, rather than re-ordered
  # due to the data.frame() call above.
  new_df[, original_order]
}
