# =============================================================================
# Clock association models
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 9, 2026
# Description: Linear models of age acceleration / DunedinPACE on a predictor
#              (diagnosis group or clinical trait) plus covariates, shared by
#              clocks/2_compare.R and clocks/3_associations.R
# =============================================================================

# Covariates that can be estimated in `data`: numeric ones, and categorical
# ones with more than one value (e.g. a single chip row among the samples)
usable_covariates <- function(data, covariates) {
  covariates[sapply(covariates, function(v) {
    x <- data[[v]]
    is.numeric(x) || length(unique(na.omit(x))) > 1
  })]
}

# Outcomes with at least one value (clocks below the coverage threshold are NA)
estimated_outcomes <- function(data, outcomes) {
  outcomes[sapply(outcomes, function(o) any(!is.na(data[[o]])))]
}

# Fit <outcome> ~ predictor + covariates for each outcome and return the
# predictor's term(s) with 95% CI, n and the model. Samples with a missing
# value in the model are dropped (lm default)
fit_predictor <- function(data, outcomes, predictor, covariates) {
  rhs <- paste(c(predictor, covariates), collapse = " + ")
  map_dfr(outcomes, function(outcome) {
    fit <- lm(as.formula(paste(outcome, "~", rhs)), data = data)
    broom::tidy(fit, conf.int = TRUE) %>%
      filter(str_starts(term, fixed(predictor))) %>%
      mutate(outcome = outcome, n = nobs(fit), covariates = rhs, .before = 1)
  })
}
