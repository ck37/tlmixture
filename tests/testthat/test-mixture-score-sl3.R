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
cluster_exposures = FALSE
verbose = TRUE
# Not currently used.
quantiles_exposures = 4L
quantiles_mixtures = 3L

result =
  tlmixture(data, outcome = "y",
            exposures = exposures,
            # This isn't being used currently.
            quantiles_exposures = quantiles_exposures,
            #quantiles_mixtures = quantiles_mixtures,
            #quantiles_mixtures = 5,
            #quantiles_mixtures = 4,
            quantiles_mixtures = 3,
            cluster_exposures = cluster_exposures,
            mixture_fn = mixture_score_sl3,
            #folds_cvtmle = folds_cvtmle,
            #folds_cvtmle = 3,
            #folds_cvtmle = 20,
            folds_cvtmle = 5,
            verbose = TRUE)

# Run manually
if (FALSE) {
  
  get_cores = function() getOption("sl.cores", RhpcBLASctl::get_num_cores())

  library(sl3)
  
  lrnr_glm <- make_learner(Lrnr_glm)
  #lrnr_glm_fast <- make_learner(Lrnr_glm_fast)
  lrnr_mean <- make_learner(Lrnr_mean)
  # Can't use ranger - doesn't support offsets (presumably).
  # lrnr_ranger100 <- make_learner(Lrnr_ranger, num.trees = 100)
  # TODO: add offset support to Lrnr_glmnet.
  lrnr_glmnet <- make_learner(Lrnr_glmnet)
  lrnr_xgb <- make_learner(Lrnr_xgboost, nrounds = 100, nthread = get_cores(), eta = 0.1)
  lrnr_xgb2 <- make_learner(Lrnr_xgboost, nrounds = 40, nthread = get_cores(), eta = 0.3, max_depth = 3L)
  lrnr_hal_d3 = make_learner(Lrnr_hal9001)
  lrnr_hal_d2 = make_learner(Lrnr_hal9001, max_degree = 2)
  lrnr_hal_d1 = make_learner(Lrnr_hal9001, max_degree = 1)
  
  # stack <- make_learner(Stack, lrnr_glm, lrnr_mean)
  sl3_screen_cor <- make_learner(Lrnr_screener_corP)
  sl3_cor_glm <- make_learner(Pipeline, sl3_screen_cor, lrnr_glm)
  
  stack <- make_learner(Stack,
                        lrnr_glm,
                        #lrnr_glmnet,
                        #sl3_cor_glm,
                        #lrnr_mean,
                        lrnr_hal_d2,
                        #lrnr_hal_d1,
                        #lrnr_glm_fast,
                        lrnr_xgb,
                        lrnr_xgb2)#, lrnr_glmnet)
  sl <- Lrnr_sl$new(learners = stack,
                    metalearner = make_learner(Lrnr_nnls, convex = TRUE))
  
  estimator_sl3 = function(...) {
    #tlmixture::
    mixture_score_sl3(...,
                                 debug = TRUE,
                                 estimator_mixture = sl$clone(),
                                 # TODO: use a different library for the confounders.
                                 estimator_confounders = sl$clone())
  }
  
  # Orthogonalize exposures.
  df2 = data
  (covariates = setdiff(names(data), c(exposures, outcome)))

  # Orthogonalize exposures per manuscript algorithm.
  for (exposure in exposures) {
    # TODO: use SuperLearner here.
    reg = glm(as.formula(paste0(exposure, " ~ .")), data = df2[, c(covariates, exposure)])
    # A' = A - E[A | W]
    df2[[exposure]] = df2[[exposure]] - reg$fitted.values
  }
  
  result =
    tlmixture(#data,
              df2,
              outcome = "y",
            exposures = exposures,
            # This isn't being used currently.
            quantiles_exposures = quantiles_exposures,
            #quantiles_mixtures = quantiles_mixtures,
            #quantiles_mixtures = 5,
            #quantiles_mixtures = 4,
            quantiles_mixtures = 3,
            cluster_exposures = cluster_exposures,
            mixture_fn = estimator_sl3,
            #folds_cvtmle = folds_cvtmle,
            #folds_cvtmle = 3,
            #folds_cvtmle = 20,
            folds_cvtmle = 5,
            verbose = TRUE)

}
result
# Weight results for the first (and only) exposure group.
result$combined$weight_dfs[[1]]
(avg_wgts = colMeans(result$combined$weight_dfs[[1]]))

library(ggplot2)

qplot(data$e1)
qplot(data$e2)

# Calculate estimated mixture from the average weights.
mixture = as.vector(as.matrix(data[, exposures]) %*% matrix(avg_wgts, ncol = 1L))

# Mixture distribution.
# TODO: different types of histogram-like plots.
qplot(mixture)

# Plot the mixture versus risk - this looks pretty good!
ggplot(data = data.frame(mixture, y = data$y),
       aes(x = mixture, y = y)) + geom_point() +
  geom_smooth(se = FALSE, size = 0.2) +
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 0.1) +
  theme_minimal() +
  labs(x = "Estimated mixture")
# TODO: add in quantile lines or something.

# TODO: Plot each mixture over the CV-TMLE folds.

result$combined$results
plot_df = result$combined$results

result$groups$mixture_objs$all$reg_mixture

# Quantile plot.
ggplot(data = plot_df, aes(x = quantile, y = psi)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  theme_minimal()

