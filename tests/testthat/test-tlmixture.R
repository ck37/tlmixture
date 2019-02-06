library(testthat)
library(tlmixture)

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

# Examine logistic regression results.
reg = glm(y ~ ., data = data, family = binomial())
summary(reg)

# Compare to OLS.
reg = lm(y ~ ., data = data)
summary(reg)

# What does partial least squares give us?
library(pls)
reg = plsr(y ~ ., data = data[, c(outcome, exposures)], ncomp = 3)
summary(reg)
# Look at explained variances.
explvar(reg)
# Unclear how to read this one.
# loadingplot(reg)
# We would ideally enforce nonnegativity on these loadings.
reg$loading.weights
reg$loadings
# TODO: confirm we should use $loadings and not $loadings.weights (unclear what the difference is).
loadings_mat = matrix(reg$loadings, nrow = nrow(reg$loadings))

# Extract the first component.
weights = loadings_mat[, 1]
# Set negative values to 0.
weights[weights < 0] = 0
# These are very close to the correct values!
weights

# Could normalize to sum to 1.
# weights = weights / sum(weights)

# tlmixture test run.
library(tlmixture)
#folds_cvtmle = 2L
folds_cvtmle = 5L
estimators = c("SL.mean", "SL.glmnet")
cluster_exposures = FALSE
verbose = TRUE
quantiles_exposures = 4L
quantiles_mixtures = 3L

result = tlmixture(data, outcome = "y",
                   exposures = exposures,
                   quantiles_exposures = quantiles_exposures,
                   #quantiles_mixtures = quantiles_mixtures,
                   quantiles_mixtures = 5,
                   estimator_outcome = estimators,
                   cluster_exposures = cluster_exposures,
                   #folds_cvtmle = folds_cvtmle,
                   #folds_cvtmle = 3,
                   folds_cvtmle = 20,
                   verbose = FALSE)

# Weight results for the first (and only) exposure group.
result$combined$weight_dfs[[1]]
(avg_wgts = colMeans(result$combined$weight_dfs[[1]]))

library(ggplot2)

qplot(data$e1)
qplot(data$e2)

# Calculate estimated mixture from the average weights.
mixture = as.vector(as.matrix(data[, exposures]) %*% matrix(avg_wgts, ncol = 1L))
# Plot the mixture versus risk.
ggplot(data = data.frame(mixture, y = data$y),
       aes(x = mixture, y = y)) + geom_point() + geom_smooth(se = FALSE) + theme_minimal() +
  labs(x = "Estimated mixture")
# TODO: add in quantile lines or something.

result$combined$results
plot_df = result$combined$results

# Quantile plot.
ggplot(data = plot_df, aes(x = quantile, y = psi)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  theme_minimal()

# Compare to weighted quantile sum.
if (FALSE) {
  if (requireNamespace("gWQS", quietly = TRUE)) {
    library(gWQS)

    # This uses only 40% of the data for training by default, and 60% for validation.
    wqs = gWQS::gwqs(y ~ x1 + x2 + x3 + x4, mix_name = exposures,
                     data = data, family = "binomial")
    # Doesn't seem that great: e1 is 0.33, e2 is 0.41, e3 is 0.26.
    wqs$final_weights


    # What if we use only 10% of data for validation (as though we did 10-fold CV)?
    wqs = gWQS::gwqs(y ~ x1 + x2 + x3 + x4, mix_name = exposures,
                     validation = 0.1,
                     data = data, family = "binomial")
    # Much better! e1 is 0.69, e2 is 0.29, e3 is 0.02.
    # But is it reporting the wrong names for e1 vs e2? It seems to have swapped them.
    print(wqs$final_weights)
  }
}
