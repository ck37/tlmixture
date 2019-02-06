#' Combine test results
#'
combine_test_results =
  function(result,
           family,
           verbose = FALSE) {


  # Hopefully we get the same number of groups for each training fold.
  # But if it's data-adaptive then we could get different results across folds.
  # TODO: talk about this scenario with Alan.
  # (Perhaps we can somehow ensure/constrain that we get the same # of groups?)
  exposure_groups = result[[1]]$exposure_groups

  weight_dfs = list()

  #browser()

  # Loop over exposure groups.
  for (group_i in seq(length(exposure_groups))) {

    # Create a dataframe with a row for each fold and each column in the weight.
    weights = data.frame(t(sapply(result, function(fold_i) {
      # Return the weights for this group, across all training folds.
      fold_i$weights[[group_i]]
    })))
    names(weights) = exposure_groups[[group_i]]

    # Save this weight dataframe.
    weight_dfs[[group_i]] = weights

    # Report on average weight across folds.
    if (verbose) {
      cat("\nWeight dataframe:\n")
      print(weights)
      cat("\nAverage normalized weights:\n")
      print(colMeans(weights))
    }

  }
  results = list(weight_dfs = weight_dfs)
  return(results)
}
