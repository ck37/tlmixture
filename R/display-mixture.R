#' Display mixture
#'
#' @export
display_mixture = function(analysis, digits = 4) {
  tab2 = analysis$tab2
  tab3 = analysis$tab3
  print(kableExtra::kable(rbind(tab2, tab3), digits = digits) %>%
          kableExtra::kable_styling(bootstrap_options = c("striped", "hover"),
                        full_width = FALSE) %>%
          kableExtra::pack_rows("Mixture", 1, 2) %>%
          kableExtra::pack_rows("Exposures", 3, nrow(tab2)) %>%
          kableExtra::pack_rows("Demographics", nrow(tab2) + 1, nrow(tab2) + nrow(tab3)))
}
