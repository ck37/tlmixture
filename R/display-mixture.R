#' Display mixture
#'
#' @param analysis tbd
#' @param digits tbd
#'
#' @export
display_mixture = function(analysis, digits = 3) {

  # Save current value of this option.
  old_value = getOption("knitr.kable.NA")
  # We need to do this to print NA's as empty strings.
  options(knitr.kable.NA = '')

  tab2 = analysis$tab2
  rownames(tab2)[1:2] = c("Mixture", "Outcome")

  tab3 = analysis$tab3

  # Set weight column to NA.
  # tab3[, 1] = NA

  kab_tab = kableExtra::kable(rbind(tab2, tab3), digits = digits) %>%
          kableExtra::kable_styling(bootstrap_options = c("striped", "hover"),
                        full_width = FALSE) %>%
          kableExtra::pack_rows("Mixture", 1, 2) %>%
          kableExtra::pack_rows("Exposures", 3, nrow(tab2)) %>%
          kableExtra::pack_rows("Adjustment variables", nrow(tab2) + 1, nrow(tab2) + nrow(tab3))


  # Restore old option value if we had to modify it.
  if (!is.null(old_value)) {
    options(knitr.kable.NA = old_value)
  }

  kab_tab
}
