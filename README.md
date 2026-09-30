
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tlmixture

**tlmixture** is an R package to construct mixtures of groups of
correlated exposures (treatments) and estimate the relationship between
the mixture and an outcome.

[![Travis build
Status](https://travis-ci.com/ck37/tlmixture.svg?token=dUHb6GuutEciSMYLsscV&branch=master)](https://travis-ci.com/ck37/tlmixture)

## Installation

You can install the development version of tlmixture from GitHub:

``` r
if (!requireNamespace("remotes")) install.packages("remotes")
remotes::install_github("ck37/tlmixture")
```

## Example

This is a simple example which shows how to use some basic function
arguments.

``` r
library(tlmixture)

# Basic example code
result =
  tlmixture(
            # Dataframe containing outcome, exposures, and adjustment variables.
            data,
            # Name of the outcome variable.
            outcome = "y",
            # Vector of exposure names (single group), or a list with separate vectors per group.
            exposures = c("exposure1", "exposure2", "exposure3")
            # This will evaluate mixtures at low/medium/high levels.
            quantiles_mixtures = 3,
            # This SuperLearner library will be used for propensity too.
            estimator_outcome = c("SL.mean", "SL.glmnet", "SL.ranger")
            # How many CV-TMLE folds to use; more is generally better, but slower to compute.
            folds_cvtmle = 3)

# Review parameter estimates and confidence intervals.
result$combined$results
```

## Clustered data

When observations are not independent (e.g. siblings within families, patients within hospitals), pass the cluster variable(s) as `id` to get cluster-robust standard errors. Columns named in `id` are used for clustering only, never as confounders.

```r
result =
  tlmixture(data,
            outcome = "y",
            exposures = c("exposure1", "exposure2", "exposure3"),
            # Name of the cluster column in data (or a vector of cluster ids).
            id = "family")

# Degrees of freedom for the t-based confidence intervals (clusters - 1).
result$cluster_df
```

With `id`, CV-TMLE folds and SuperLearner's internal folds keep clusters together, influence-curve variances are summed within clusters, and confidence intervals use a t distribution with (clusters - 1) degrees of freedom. Several clustering variables can be given, e.g. `id = c("school", "neighborhood")`; a variable nested in another (family within site) is handled by clustering on the outer one.

### Customizing the mixture function

To change a mixture function's settings, wrap it and forward `...`. Forwarding `...` also passes the cluster ids through, so the mixture's own cross-validation keeps clusters together:

```r
library(sl3)

my_mixture = function(...) {
  mixture_score_sl3(
    ...,
    estimator_mixture = Lrnr_sl$new(
      learners = Stack$new(Lrnr_glm$new(), Lrnr_xgboost$new(nrounds = 100)),
      metalearner = Lrnr_nnls$new(convex = TRUE)),
    estimator_confounders = Lrnr_sl$new(
      learners = Stack$new(Lrnr_glm$new(), Lrnr_ranger$new(num.trees = 500)),
      metalearner = Lrnr_nnls$new(convex = TRUE)))
}

result = tlmixture(data, outcome = "y", exposures = c("exposure1", "exposure2"),
                   mixture_fn = my_mixture, id = "family")
```