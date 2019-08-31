#' Analyze mixture
#'
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

    library(dplyr)

    tab = data %>%
      filter(!is.na(mixture)) %>%
      group_by(quantile) %>%
      dplyr::select(one_of(c(exposures, outcome)), mixture, quantile) %>%
      dplyr::mutate(n = n()) %>%
      mutate_at(vars(one_of(c(exposures, outcome)), mixture),
                list(mean = ~ mean(., na.rm = TRUE)))  %>%
      unique

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

    if (!is.null(vars_desc)) {
      tab3 = data %>%
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

    # Plot mixture quantiles.
    print(qplot(mixture) +
            geom_vline(aes(xintercept = quants_plot)) +
            theme_bw() +
            labs(title = paste("Mixture distribution:", name)))


    # General Lowess smooth
    print(ggplot(data = data.frame(mixture, y = data[[outcome]]),
                 aes(x = mixture, y = y)) + #, weight = weight_var)) +
            #aes_(x = "mixture", y = "y")) + #, weight = weight_var)) +
            geom_point() + geom_smooth(se = TRUE, span = span) +
            geom_vline(aes(xintercept = quants_plot),
                       data = data.frame(quants_plot)) +
            theme_minimal() +
            labs(x = paste("Estimated mixture:", name),
                 title = paste("Unadjusted risk:", name)))

    plot_df = tlmixture_result$combined$results

    # Quantile plot.
    print(ggplot(data = plot_df, aes(x = quantile, y = psi)) +
            geom_point() +
            geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
            theme_minimal())

    reg_result = try({
      reg_str = paste(outcome, " ~ mixture + ", paste(reg_vars, collapse = " + "))
      cat(reg_str, "\n")
      reg_form = as.formula(reg_str)

      reg = glm(reg_form, data = data, family = binomial())
      print(summary(reg))
    })

    result = list(tab2 = tab2,
                  tab3 = tab3,
                  mixture = mixture,
                  weights = avg_wgts,
                  name = name,
                  reg = reg)
    return(result)
  }
