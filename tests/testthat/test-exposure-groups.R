library(testthat)
library(tlmixture)
library(ggplot2)

# Create a simple simulated dataset.

samples = 100L
means = 0:3
sigma = diag(1:4)

# Add correlation between 1 and 3
sigma[1, 3] = sigma[3, 1] = 0.2
# And between 2 and 4.
sigma[2, 4] = sigma[4, 2] = 0.1
sigma

set.seed(1, "L'Ecuyer-CMRG")
adjustment = data.frame(MASS::mvrnorm(samples, mu = means, Sigma = sigma))
names(adjustment) = tolower(names(adjustment))

colMeans(adjustment)

# Adjustment variables:
# 1: strong confounder through exposure 1
# 2: weak confounder through exposure 2
# 3: related to exposures but not outcome
# 4: related to outcome but not exposures

#########
# Create exposures

exposure1 = 0.6 * adjustment$x1 + 0.2 * adjustment$x3 + 0.2 * rnorm(samples)
exposure2 = 0.2 * adjustment$x2 + 0.2 * adjustment$x3 + 0.6 * rnorm(samples)
exposure3 = 0.1 * adjustment$x1 + 0.1 * adjustment$x3 + 0.8 * rnorm(samples)

exposure_data = data.frame(exposure1, exposure2, exposure3)
names(exposure_data) = paste0("e", 1:3)
names(exposure_data)

exposures = names(exposure_data)

################
# Exposure relations to outcome.
# 1: highly impact on outcome
# 2: low impact on outcome
# 3: no impact on outcome

outcome_signal =
  with(adjustment,
       with(exposure_data,
            e1 * 0.9 + e2 * 0.5 + x1 * 0.2 + x2 * 0.1 + x4 * 0.2))

# Center to mean 0.
outcome_signal = outcome_signal - mean(outcome_signal)

y = rbinom(samples, 1, plogis(outcome_signal))
summary(y)

################
# Combine into dataset.
data = data.frame(y, exposure_data, adjustment)
outcome = "y"


# tlmixture test run.
library(tlmixture)
#folds_cvtmle = 2L
folds_cvtmle = 5L
estimators = c("SL.mean", "SL.glmnet")
cluster_exposures = FALSE
verbose = TRUE
# Not currently used.
quantiles_exposures = 4L
quantiles_mixtures = 3L

# Two exposure groups.
exposure_groups = list(c("e1", "e2"), c("e2", "e3"))

result =
  tlmixture(data, outcome = "y",
            exposures = exposure_groups,
            # This isn't being used currently.
            quantiles_exposures = quantiles_exposures,
            #quantiles_mixtures = quantiles_mixtures,
            #quantiles_mixtures = 5,
            #quantiles_mixtures = 4,
            quantiles_mixtures = 3,
            estimator_outcome = estimators,
            cluster_exposures = cluster_exposures,
            #mixture_fn = tlmixture::mixture_sl,
            mixture_fn = tlmixture::mixture_pls,
            #folds_cvtmle = folds_cvtmle,
            #folds_cvtmle = 3,
            #folds_cvtmle = 20,
            folds_cvtmle = 5,
            verbose = FALSE)

names(result)
#result

result$combined$results


# Weight results for the first (and only) exposure group.
result$combined$weight_dfs[[1]]
(avg_wgts_1 = colMeans(result$combined$weight_dfs[[1]]))
(avg_wgts_2 = colMeans(result$combined$weight_dfs[[2]]))

library(ggplot2)

qplot(data$e1)
qplot(data$e2)

# Calculate estimated mixture from the average weights.
mixture_1 = as.vector(as.matrix(data[, exposure_groups[[1]]]) %*% matrix(avg_wgts_1, ncol = 1L))
mixture_2 = as.vector(as.matrix(data[, exposure_groups[[2]]]) %*% matrix(avg_wgts_2, ncol = 1L))

# Mixture distribution.
# TODO: different types of histogram-like plots.
qplot(mixture_1)
qplot(mixture_2)
qplot(mixture_2, data$e2)

# Mixture 2 is just exposure 2.
cor(mixture_2, data$e2)
cor.test(mixture_2, data$e2)

# Plot the mixture versus risk - this looks pretty good!
ggplot(data = data.frame(mixture_1, y = data$y),
       aes(x = mixture_1, y = y)) + geom_point() +
  geom_smooth(se = FALSE, size = 0.2) +
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 0.1) +
  theme_minimal() +
  labs(x = "Estimated mixture 1")

# Plot the mixture versus risk - this looks pretty good!
ggplot(data = data.frame(mixture_2, y = data$y),
       aes(x = mixture_2, y = y)) + geom_point() +
  geom_smooth(se = FALSE, size = 0.2) +
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 0.1) +
  theme_minimal() +
  labs(x = "Estimated mixture 2")

# TODO: add in quantile lines or something.

# TODO: Plot each mixture over the CV-TMLE folds.

result$combined$results
plot_df = result$combined$results

# Separate quantile plots for the two exposure groups.
ggplot(data = subset(plot_df, exposure_group == 1), aes(x = quantile, y = psi)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  theme_minimal()

ggplot(data = subset(plot_df, exposure_group == 2), aes(x = quantile, y = psi)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  theme_minimal()
