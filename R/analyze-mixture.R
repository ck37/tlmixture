#' Analyze mixture
#'
#' @param data tbd
#' @param tlmixture_result tbd
#' @param plots tbd
#' @param span tbd
#' @param vars_desc tbd
#' @param name tbd
#' @param weight_var tbd
#' @param round_digits tbd
#' @param rescale_adjustment tbd
#' @param reg_vars tbd
#'
#' @importFrom dplyr filter group_by one_of n mutate_at vars summarize_at select
#' @importFrom ggplot2 ggplot geom_point geom_smooth theme_minimal geom_errorbar aes labs theme qplot geom_vline theme_bw
#' @export
analyze_mixture =
  function(data, tlmixture_result,
           # TODO: make plot_mixture() a separate function.
           plots = TRUE,
           span = 0.9,
           # Descriptive variables, e.g. demographics.
           vars_desc = NULL,
           name = "exposures",
           weight_var = NULL,
           round_digits = 4,
           rescale_adjustment = TRUE,
           reg_vars = NULL) {

    # Extract key hyperparameters from the tlmixture result object.
    exposures = tlmixture_result$exposures
    num_quantiles = tlmixture_result$quantiles_mixtures
    outcome = tlmixture_result$outcome

    # Calculate overall weights.
    avg_wgts = colMeans(tlmixture_result$combined$weight_dfs[[1]])

    # Calculate estimated mixture from the average weights.
    mixture = as.vector(as.matrix(data[, exposures]) %*% matrix(avg_wgts, ncol = 1L))

    # Calculate mixture quantiles.
    # Extract the internal values - don't need the 0% and 100% marks.
    (quants_full = quantile(mixture,
                            probs = seq(0, 1, length.out = num_quantiles + 1),
                            na.rm = TRUE))
    (quants_plot = quants_full[2:(num_quantiles)])

    data$mixture = mixture
    data$quantile = Hmisc::cut2(mixture, cuts = quants_full)

    # First table examines the mixture and the outcome across the mixture quantiles.
    # This isn't actually used anymore.
    if (FALSE) {
    tab = data %>%
      filter(!is.na(mixture)) %>%
      group_by(quantile) %>%
      dplyr::select(one_of(c(exposures, outcome)), mixture, quantile) %>%
      dplyr::mutate(n = n()) %>%
      mutate_at(vars(one_of(c(exposures, outcome)), mixture),
                list(mean = ~ mean(., na.rm = TRUE)))  %>%
      unique#()
   # %>% as.data.frame()
    }

    # Rename outcome to just be "Outcome"
    #rownames(tab) = c("Mixture", "Outcome")

    # TODO: combine these into a data.frame so that we don't have to order
    # separate vectors.
    ordered_exposures = exposures[order(avg_wgts, decreasing = TRUE)]
    ordered_weights = avg_wgts[order(avg_wgts, decreasing = TRUE)]

    tab2 = data %>%
      filter(!is.na(mixture)) %>%
      group_by(quantile) %>%
      dplyr::select(one_of(c(ordered_exposures, outcome)), mixture, quantile) %>%
      dplyr::mutate(n = n()) %>%
      summarize_at(vars("mixture", outcome, ordered_exposures),
                   list(~mean(., na.rm = TRUE))) %>%
      select(-quantile) %>%
      #mutate(quantile = as.numeric(quantile)) %>%
      as.data.frame()

    #tab2$quantile = as.character(tab2$quantile)

    # Transpose
    tab2 = data.frame(t(tab2))

    # Add quantile number to the labels.
    names(tab2) = paste(paste0(1:length(levels(data$quantile)), "."),
                        levels(data$quantile))

    # TODO: make this a separate table.
    tab2 = cbind("Weights" = c(NA, NA, ordered_weights), tab2)

    tab2 = round(tab2, round_digits)

    tab3 = NULL

    # Examine adjustment vars over the quantiles of the mixture.
    if (!is.null(vars_desc)) {
      if (rescale_adjustment) {
        # Rescale to mean 0, std dev. 1
        data2 = ck37r::standardize(data,
                                   skip_vars = setdiff(names(data), vars_desc))
      } else {
        data2 = data
      }

      tab3 = data2 %>%
        # TODO: report on this missingness.
        filter(!is.na(mixture)) %>%
        group_by(quantile) %>%
        dplyr::select(one_of(vars_desc), quantile) %>%
        dplyr::mutate(n = n()) %>%
        summarize_at(vars(vars_desc),
                     list(~mean(., na.rm = TRUE))) %>%
        select(-quantile) %>%
        #mutate(quantile = as.numeric(quantile)) %>%
        as.data.frame()

      #tab2$quantile = as.character(tab2$quantile)

      # Transpose
      tab3 = data.frame(t(tab3))

      # Add quantile number to the labels.
      names(tab3) = paste(paste0(1:length(levels(data$quantile)), "."),
                          levels(data$quantile))

      #browser()

      # TODO: make this a separate table.
      tab3 = cbind("Weights" = rep(0, nrow(tab3)), tab3)
      #tab3$Weights = NA

      tab3 = round(tab3, round_digits)
    }

    ######################
    # Plots

    # Plot mixture distribution with quantiles.
    print({mixture_dist = qplot(mixture, fill = I("gray70"), show.legend = FALSE) +
            geom_vline(aes(xintercept = quants_plot)) +
            theme_minimal() +
            labs(#title = paste("Mixture distribution:", name),
                 y = "Frequency",
                 x = "Estimated mixture")})


    # General Lowess smooth
    print({unadjusted_smooth = ggplot(data = data.frame(mixture, y = data[[outcome]]),
                 aes(x = mixture, y = y)) + #, weight = weight_var)) +
            #aes_(x = "mixture", y = "y")) + #, weight = weight_var)) +
            geom_point(alpha = I(0.5), show.legend = FALSE, stroke = 0) +
            geom_smooth(se = TRUE, span = span) +
            geom_vline(aes(xintercept = quants_plot),
                       data = data.frame(quants_plot), alpha = I(0.8)) +
            theme_minimal() +
            labs(x = paste("Estimated mixture"),
                 y = "Outcome"#,
                 #title = paste("Unadjusted relationship")
    )})
    #print(unadjusted_smooth)

    plot_df = tlmixture_result$combined$results

    # Quantile plot.
    print({adjusted_effects = ggplot(data = plot_df, aes(x = quantile, y = psi)) +
            geom_point() +
            geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
            labs(x = "Mixture quantile",
                 y = "Adjusted outcome mean") +
            theme_minimal()})
    #print(adjusted_effects)

    reg_result = try({
      reg_str = paste(outcome, " ~ mixture + ", paste(reg_vars, collapse = " + "))
      cat(reg_str, "\n")
      reg_form = as.formula(reg_str)

      reg = glm(reg_form, data = data, family = binomial())
      print(summary(reg))
    })

    result =
      list(tab2 = tab2,
           tab3 = tab3,
           mixture = mixture,
           weights = avg_wgts,
           plots = list(
             mixture_dist,
             unadjusted_smoooth = unadjusted_smooth,
             adjusted_effects = adjusted_effects),
           name = name,
           num_quantiles = num_quantiles, # needed for plot_analysis()
           reg = reg)
    return(result)
  }
