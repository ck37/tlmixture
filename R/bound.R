# CK: adapted from tmle package.
#---------- function .bound ---------------
# set outliers to min/max allowable values
# assumes data_vec contains only numerical data
#-----------------------------------------
#' @export
bound =
  function(data_vec, bounds) {
  data_vec[data_vec > max(bounds)] = max(bounds)
  data_vec[data_vec < min(bounds)] = min(bounds)
  return(data_vec)
}
