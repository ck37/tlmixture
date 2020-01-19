#' Analyze results for the overall exposure groups
#'
#' Not relevant if there is only a single mixture
#'
#' @param combined_results tbd
#' @param data tbd
#' @param outcome tbd
#' @param exposures tbd
#' @param mixture_fn tbd
#' @param refit_mixtures tbd
#' @param verbose tbd
#'
#' @importFrom dplyr arrange
#' @importFrom magrittr %>%
analyze_exposure_groups = function(combined_results,
                                   data,
                                   outcome,
                                   exposures,
                                   mixture_fn,
                                   refit_mixtures = FALSE,
                                   verbose = FALSE) {
  
  # Convert the exposures vector to a list.
  # TODO: do this earlier in tlmixture()
  if (!is.list(exposures)) {
    exposures = list(exposures)
    names(exposures) = "all"
  }

  groups_df = combined_results$groups_df
  # Set to the raw pval by default, in case there is only 1 group.
  groups_df$rd_pval_fdr = groups_df$rd_val

  if (nrow(groups_df) > 1L) {

    # Calculate Benjamini-Hochberg FDR-adjusted p-values.
    groups_df$rd_pval_fdr = stats::p.adjust(groups_df$rd_pval, method = "BH")

    # Sort by p-value (asc), then standard err (asc)

    groups_df = groups_df %>% arrange(rd_pval_fdr, rd_pval, rd_se) %>% as.data.frame()
  }
  
  cat("Refitting mixtures to full data.\n")
  
  mixture_df = NULL
  mixture_objs = list()
  
  if (refit_mixtures) {
  
    # Loop over exposure groups and fit mixture to full dataset.
    for (group_i in seq(length(exposures))) {
      # Fit mixture to full dataframe.
      # TODO: also save mix obj.
      mix_obj = mixture_fn(data, outcome, exposures[[group_i]], verbose = verbose)
      mixture_objs[[group_i]] = mix_obj
      
      # Predict on full df.
      mixture = predict(mix_obj, data)
      
      # TODO: cbind all at once, which will be more efficient.
      if (!is.null(mixture_df)) {
        mixture_df = cbind.data.frame(mixture_df, mixture)
      } else {
        mixture_df = data.frame(mixture)
      }
    }
    
    names(mixture_objs) = names(exposures)
    names(mixture_df) = names(exposures)
  }
  

  # Return the updated dataframe.
  results = list(
    groups_df,
    mixture_df = mixture_df,
    mixture_objs = mixture_objs)
}
