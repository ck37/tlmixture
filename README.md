
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tlmixture

**tlmixture** is an R package to construct mixtures of groups of
correlated exposures (treatments) and estimate the relationship between
the mixture and an outcome.

[![Travis build
Status](https://travis-ci.com/ck37/tlmixture.svg?token=dUHb6GuutEciSMYLsscV&branch=master)](https://travis-ci.com/ck37/tlmixture)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/ck37/tlmixture?branch=master&svg=true)](https://ci.appveyor.com/project/ck37/tlmixture)

## Installation

You can install the development version of tlmixture from GitHub:

``` r
install.packages("remotes")
remotes::install_github("ck37/tlmixture")
```

## Example

This is a simple example which shows how to use some basic function
arguments.

``` r
library(tlmixture)

# Basic example code
result =
            # dataframe containing outcome, exposures, and adjustment variables.
  tlmixture(data,
            # Name of the outcome variable.
            outcome = "y",
            # Vector of exposure names (single group), or a list with separate vectors per group.
            exposures = c("exposure1", "exposure2", "exposure3")
            # This will evaluate mixtures at low/medium/high levels.
            quantiles_mixtures = 3,
            # This SuperLearner library will be used for propensity too.
            estimator_outcome = c("SL.mean", "SL.glmnet", "SL.ranger")
            # How many CV-TMLE folds to use; more is better, but slower to compute.
            folds_cvtmle = 5)

# Review parameter estimates and confidence intervals.
result$combined$results
```
