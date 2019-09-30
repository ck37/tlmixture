#' Combine test results
#'
#' @param result tbd
#' @param family tbd
#' @param quantiles_mixtures tbd
#' @param verbose tbd
#'
#' @importFrom stats glm qlogis plogis var
#'
combine_test_results =
  function(result,
           family,
           quantiles_mixtures,
           verbose = FALSE) {


  # Hopefully we get the same number of groups for each training fold.
  # But if it's data-adaptive then we could get different results across folds.
  # TODO: talk about this scenario with Alan.
  # (Perhaps we can somehow ensure/constrain that we get the same # of groups?)
  exposure_groups = result[[1]]$exposure_groups

  weight_dfs = list()

  #browser()

  # Dataframe to contain results.
  results_df = NULL

  # Loop over exposure groups.
  for (group_i in seq(length(exposure_groups))) {

    ###############
    #  Create a dataframe with a row for each fold and each column in the weight.

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

    ################

    # Loop over mixture quantiles
    for (quantile_i in seq(quantiles_mixtures)) {

      # Extract the dataframes for this quantile across all CV-folds.
      test_results =
        do.call(rbind, lapply(seq(length(result)), function(fold_i) {
          fold = result[[fold_i]]

          # Return the weights for this group, across all training folds.
          df = fold$test_results[[group_i]]

          # Restrict to the current quantile that we're analyzing
          df = df[df$quantile == quantile_i, ]

          # Save the fold for future reference.
          df$fold = fold_i

          # Return the dataframe.
          df
      }))

      # This dataframe should be the same size as the original dataset.
      #stopifnot(nrow(test_results) == nrow(data))

      ########################
      # TODO: implement CV-TMLE

      # From tmle::estimateQ.
      # TODO: fix warning "NaNs produced"
      test_results$logit_q_hat = qlogis(test_results$q_pred)
      # TODO: stop if all results are NAs.

      # Estimate epsilon
      # This is for a logistic fluctuation.
      # TODO: try alternative version at
      # https://github.com/ck37/varimpact/blob/master/R/estimate_pooled_results.R#L96-L106
      reg = try(glm(y ~ -1 + offset(logit_q_hat) + haw,
                    data = test_results, family = "binomial"))
      if ("try-error" %in% class(reg)) {
        cat("Error in epsilon rgression.\n")
        browser()
      }

      epsilon = try(coef(reg))

      # Make sure that epsilon isn't NA.
      stopifnot(!is.na(epsilon))

      # Fluctuate Q_star
      q_star = test_results$logit_q_hat + epsilon * test_results$h1w

      # Transform Q_star
      q_star = plogis(q_star)

      if (verbose) cat("Estimating per-fold thetas: ")
      # Estimate parameter on every validation fold.
      thetas = tapply(q_star, test_results$fold, mean, na.rm = TRUE)
      if (verbose) cat(thetas, "\n")

      # Move Q_star into the data so that it can be analyzed per-fold.
      test_results$q_star = q_star
      rm(q_star)

      # Calculate ICs (per fold)
      if (verbose) cat("Calculating per-fold influence curves\n")
      # Get influence curve per fold.
      # Influence_curves here is a list, where each element is a result.
      # We can't convert to a matrix (one IC per col) because lengths may be different.
      # But perhaps it could be a long dataframe instead, with an extra column to denote fold.
      # TODO: figure out why this can generate NaNs
      influence_curves =
        base::by(test_results, test_results$fold, function(fold_data) {
        if (FALSE && verbose) {
          with(fold_data,
               cat("A:", length(A), "g1W_hat:", length(g1W_hat), "Y_star:", length(Y_star),
                   "Q_star:", length(Q_star), "\n"))
        }
        # (Here in_quantile is A).
        result = with(fold_data, (in_quantile / g_pred) * (y - q_star) +
                        q_star - mean(q_star, na.rm = TRUE))
        result
      })

      # Calculate fold sizes.

      # Calculate IC variances.
      ic_vars = sapply(influence_curves, var)

      # Calculate SEs
      std_errs = lapply(influence_curves, function(ic) {
        var(ic) / length(ic)
      })

      # Calculate combined SE.
      std_err = sqrt(mean(ic_vars) / nrow(test_results))

      # Calculate exposure-specific mean.
      psi = mean(thetas)

      # Calculate CI
      ci_lower = psi - 1.96 * std_err
      ci_upper = psi + 1.96 * std_err

      # Calculate P-value (two-sided).

      # Compile results into a one-row dataframe.
      new_results =
        data.frame(exposure_group = group_i,
                   quantile = quantile_i,
                   psi = psi,
                   std_err = std_err,
                   # TODO: fill in these values.
                   ci_lower = ci_lower,
                   ci_upper = ci_upper,
                   p_value = NA)

      # Add results to dataframe.
      results_df = rbind(results_df, new_results)

    }


  }
  results = list(weight_dfs = weight_dfs,
                 results = results_df
                 # TODO: what other results?
                 )
  return(results)
}
