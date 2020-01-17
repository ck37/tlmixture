#' Plot analysis
#'
#' @param analysis tbd
#'
#' @importFrom dplyr mutate if_else
#' @importFrom ggplot2 ggplot aes geom_line theme_minimal theme element_blank expand_limits
#' @export
plot_analysis = function(analysis) {

  # Extract demographic data.
  if (!is.null(analysis$tab3)) {
    df = analysis$tab3
    df$Weights = NULL
    names(df) = paste0(seq(analysis$num_quantiles))
    df = cbind("var" = rownames(df), df)
    rownames(df) = NULL
    df

    # Convert to long format for plotting.
    df_long = df %>% tidyr::gather(quartile, mean, as.character(seq(analysis$num_quantiles)))

    # Put labels at the end of the lines, via:
    # https://stackoverflow.com/questions/29357612/plot-labels-at-ends-of-lines

    g = df_long %>%
      mutate(label = if_else(quartile == max(quartile), as.character(var), NA_character_)) %>%
      ggplot(aes(x = quartile, y = mean, group = var, color = var)) +
      geom_line() + theme_minimal() +
      labs(#title = paste("Adjustment variables"),
           y = "Mean within quantile",
           x = "Mixture quantile") +
      expand_limits(x = 5) +
      #scale_x_discrete(expand = c(0, 4)) +
      ggrepel::geom_label_repel(aes(label = label),
                                #ggrepel::geom_text_repel(aes(label = label),
                                nudge_x = 0.2,
                                direction = "both",
                                na.rm = TRUE) +
      theme(legend.title = element_blank(),
            legend.position = "none")#,
    #plot.margin = unit(c(1,6,1,1), "lines"))
    print(g)
  }
  
  result = list(
    demos = g
  )
  
  invisible(result)
}
