#' Analyze results for the overall exposure groups
#'
#' Not relevant if there is only a single mixture
#'
#' @param combined_results tbd
#' @param verbose tbd
#'
#' @importFrom dplyr arrange
#' @importFrom magrittr %>%
analyze_exposure_groups = function(combined_results, verbose = FALSE) {

  groups_df = combined_results$groups_df
  # Set to the raw pval by default, in case there is only 1 group.
  groups_df$rd_pval_fdr = groups_df$rd_val

  if (nrow(groups_df) > 1L) {

    # Calculate Benjamini-Hochberg FDR-adjusted p-values.
    groups_df$rd_pval_fdr = stats::p.adjust(groups_df$rd_pval, method = "BH")

    # Sort by p-value (asc), then standard err (asc)

    groups_df = groups_df %>% arrange(rd_pval_fdr, rd_se) %>% as.data.frame()
  }

  # Return the updated dataframe.
  groups_df
}
