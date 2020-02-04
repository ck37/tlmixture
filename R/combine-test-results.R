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

  if (verbose) {
    cat("\nCombining test results.\n")
  }


  # Hopefully we get the same number of groups for each training fold.
  # But if it's data-adaptive then we could get different results across folds.
  # TODO: talk about this scenario with Alan.
  # (Perhaps we can somehow ensure/constrain that we get the same # of groups?)
  exposure_groups = result[[1]]$exposure_groups

  weight_dfs = list()

  #browser()

  # Dataframe to contain results.
  results_df = NULL

  # Dataframe to contain exposure-group results.
  groups_df = NULL

  # Loop over exposure groups.
  for (group_i in seq(length(exposure_groups))) {

    # Attempt to recover from any errors that occur during this analysis.
    # If an exposure group had a zero-variation mixture function then all
    # of this analysis will fail.
    tryCatch({
    ###############
    #  Create a dataframe with a row for each fold and each column in the weight.

    # This doesn't work for SL currently.
    if ("weights" %in% names(result[[1]])) {
      weights = data.frame(t(sapply(result, function(fold_i) {
        # Return the weights for this group, across all training folds.
        if (length(fold_i$weights) > 0) {
          fold_i$weights[[group_i]]
        } else {
          NULL
        }
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
    } else {
      # Estimators that don't generate parametric vector of weights or coefficients.
      weight_dfs[[group_i]] = NULL
    }

    ################

    # Non-dataframe list to save arbitrary objects for each quantile.
    quantile_results = list()

    # Loop over mixture quantiles
    for (quantile_i in seq(quantiles_mixtures)) {

      # Extract the dataframes for this quantile across all CV-folds.
      test_results =
        do.call(rbind, lapply(seq(length(result)), function(fold_i) {
          fold = result[[fold_i]]

          # Check if estimation failed for this group in this fold.
          if (!group_i %in% names(fold$test_results)) {
            return(NULL)
          }

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
      # For continuous variables this will yield a warning:
      # "In eval(family$initialize) : non-integer #successes in a binomial glm!"
      reg = try(suppressWarnings(glm(y ~ -1 + offset(logit_q_hat) + haw,
                    data = test_results, family = "binomial")))
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

      # Calculate quantile-specific mean.
      psi = mean(thetas)

      # Calculate CI
      ci_lower = psi - 1.96 * std_err
      ci_upper = psi + 1.96 * std_err

      # Calculate P-value (two-sided).

      # Compile results into a one-row dataframe.
      new_results =
        data.frame(exposure_group = group_i,
                   quantile = quantile_i,
                   # Quantile-specific mean, averaged over all CV-TMLE folds.
                   psi = psi,
                   std_err = std_err,
                   # TODO: fill in these values.
                   ci_lower = ci_lower,
                   ci_upper = ci_upper,
                   p_value = NA)

      # Add results to dataframe.
      results_df = rbind(results_df, new_results)

      # Also save objects needed for the RD and RR effect calculations.
      results = list(
        # Adjusted quantile-specific means, one per CV-TMLE fold
        means = thetas,
        # This is a list, with one element per CV-TMLE fold
        curves = influence_curves
      )
      quantile_results[[quantile_i]] = results


    } # Looping over quantiles of the current mixture

    #####################################################
    # Now that we have the quantile-specific, estimate group-specific parameters:
    # 1. Risk difference of high quantile - low quantile, with 95% CI
    # 2. Risk ratio of high quantile / low quantile, with 95% CI
    #   - But only if outcome is binary?

    # Extract the quantile results for convenience access.
    qmin = quantile_results[[1]]

    # quantile_i should already be set to the maximum quantile due to the FOR loop above.
    qmax = quantile_results[[quantile_i]]

    # 1. Risk difference
    # Calculate estimate within each fold
    # One element per fold.
    rd_vec = qmax$means - qmin$means

    # Calculate variance of difference of influence curves for each CV-TMLE fold.
    # One element per fold.
    rd_ic_var = sapply(seq(length(qmax$curves)), function(ic_i) {
      var(qmax$curves[[ic_i]] - qmin$curves[[ic_i]])
    })

    n_validation = sapply(seq(length(qmax$curves)), function(ic_i) {
      # Picking the maximum quantile arbitrarily, the length of the influence curve
      # tells us the sample size in that validation set.
      length(qmax$curves[[ic_i]])
    })

    # Calculate the standard error for each CV-TMLE fold, based on the size of the
    # validation dataset
    # This will be a vector of: sqrt(varIC / n_validation)
    rd_se = sqrt(rd_ic_var / n_validation)

    if (verbose) {
      cat("Risk differences:", rd_vec, "\n")
      cat("Variances:", rd_ic_var, "\n")
      cat("Standard errors:", rd_se, "\n")
    }

    total_validation_size = sum(n_validation)

    # Calculate overall RD, se, and p-value
    rd_overall = mean(rd_vec)

    # Overall standard error is sqrt(mean(IC variances) / n) where n = total size of all validation sets.
    rd_se_overall = sqrt(mean(rd_ic_var) / total_validation_size)

    # 1-sided p-value
    # rd_pval_overall = 1 - pnorm(rd_overall / rd_se_overall)
    # 2-sided p-value
    rd_pval_overall = 2 * pnorm(-abs(rd_overall / rd_se_overall))

    # Confidence interval
    z_1.96 = qnorm(0.975)
    rd_ci_lower = rd_overall - z_1.96 * rd_se_overall
    rd_ci_upper = rd_overall + z_1.96 * rd_se_overall

    # TODO: also do for RR

    if (verbose) {
      cat("Risk difference:", rd_overall,
          paste0(" (", round(rd_ci_lower, 3), "-", round(rd_ci_upper, 3), ")"),
          paste0(" p = ", round(rd_pval_overall, 5)),
          "\n")
    }

    # Include group_name if it's defined.
    if (!is.null(names(exposure_groups))) {
      group_name = names(exposure_groups)[group_i]
    } else {
      group_name = ""
    }

    # Compile results for this exposure group.
    group_result = list(
      i = group_i,
      # Not working if group name is null - needs to be an empty string.
      group = group_name,
      rd = rd_overall,
      rd_pval = rd_pval_overall,
      rd_ci_lower = rd_ci_lower,
      rd_ci_upper = rd_ci_upper,
      rd_se = rd_se_overall
    )

    if (verbose) {
      cat("Result for group", group_i, "\n")
      print(group_result)
    }

    # Append to the dataframe.
    groups_df = rbind.data.frame(groups_df, group_result,
                                 # We need strings to not be factors because we
                                 # are adding a unique group name at each iteration.
                                 stringsAsFactors = FALSE)
    }, error = function(error) {
      cat("Analysis failed for group", group_i, "\n")
      print(error)
    })

  } # Looping over exposure groups

  results = list(#weight_dfs = NULL, #weight_dfs,
                 weight_dfs = weight_dfs,
                 results = results_df,
                 groups_df = groups_df
                 # TODO: what other results?
                 )
  return(results)
}
