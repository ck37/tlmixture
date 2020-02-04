#' Display mixture
#'
#' @param analysis tbd
#' @param digits tbd
#'
#' @export
display_mixture =
  function(analysis, digits = 3,
           caption = "Mixture analysis",
           label = NULL,
           booktabs = TRUE) {

  # Save current value of this option.
  old_value = getOption("knitr.kable.NA")
  # We need to do this to print NA's as empty strings.
  options(knitr.kable.NA = '')


  # TODO: add adjusted outcome
  tab2 = analysis$tab2
  rownames(tab2)[1:2] = c("Mixture", "Outcome")

  tab3 = analysis$tab3

  # Set weight column to NA.
  # tab3[, 1] = NA

  # TODO: labe mixture quantiles as the columns.

  kab_tab = kableExtra::kable(rbind(tab2, tab3), digits = digits, booktabs = booktabs,
                              caption = caption, label = label) %>%
          kableExtra::kable_styling(bootstrap_options = c("striped", "hover"),
                        full_width = FALSE) %>%
          kableExtra::add_header_above(c(" " = 1, "Mean by mixture quantile" = ncol(tab2))) %>%
          #kableExtra::pack_rows("Mixture", 1, 2) %>%
          kableExtra::pack_rows("Exposures", 3, nrow(tab2)) %>%
          kableExtra::pack_rows("Covariates", nrow(tab2) + 1, nrow(tab2) + nrow(tab3))


  # Restore old option value if we had to modify it.
  if (!is.null(old_value)) {
    options(knitr.kable.NA = old_value)
  }

  kab_tab
}
