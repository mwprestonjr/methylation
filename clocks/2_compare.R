# =============================================================================
# Clocks 2: compare age acceleration between groups
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 8, 2026
# Description: Compares age acceleration and DunedinPACE between CLOCK_GROUPS
#              (first = reference), adjusted for CLOCK_COVARIATES
#              (clocks/config.R), using the per-sample table from
#              1_estimate.R. Saves the group effect per outcome for 3_plots.R.
#              Skipped if a group has fewer than CLOCK_MIN_GROUP_N samples
#              (e.g. a data source without controls)
# Usage:       Rscript clocks/2_compare.R <data source>   (see config.R)
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------
library(tidyverse)

# Load shared configuration
source("config.R")
source("clocks/config.R")

# combat_batch is a label (plate or chip number), not a number
clocks <- read.csv(CLOCK_RESULTS,
                   colClasses = c(GP2ID = "character", combat_batch = "character"))
cat("Samples:", nrow(clocks), "\n")

# --- 1. Samples in the comparison --------------------------------------------

# CLOCK_GROUPS, with the first as the reference level
clocks_cmp <- clocks %>%
  filter(phenotype %in% CLOCK_GROUPS) %>%
  mutate(phenotype = factor(phenotype, levels = CLOCK_GROUPS))
group_label <- paste(rev(CLOCK_GROUPS), collapse = " vs ")
group_sizes <- table(clocks_cmp$phenotype)
cat("\nSamples in", group_label, "comparison:\n")
print(group_sizes)

# --- 2. Adjusted group differences -------------------------------------------

if (all(group_sizes >= CLOCK_MIN_GROUP_N)) {
  # Every outcome gets the same covariates, including age: the AgeAccel
  # residuals are uncorrelated with age across all samples, but not necessarily
  # within this subset, and adjusting in the same model handles an age
  # difference between the groups correctly. Categorical covariates with a
  # single value here can't be estimated and are left out
  usable <- CLOCK_COVARIATES[sapply(CLOCK_COVARIATES, function(v) {
    x <- clocks_cmp[[v]]
    is.numeric(x) || length(unique(na.omit(x))) > 1
  })]
  covars <- paste(c("phenotype", usable), collapse = " + ")
  cat("\nAge acceleration,", group_label, "- model: <outcome> ~", covars, "\n")
  if (length(setdiff(CLOCK_COVARIATES, usable))) {
    cat("Covariates left out (one value only):", setdiff(CLOCK_COVARIATES, usable), "\n")
  }

  # Outcomes not estimated on this array (all NA, below the coverage
  # threshold) are left out
  outcomes <- c(paste0("AgeAccel_", ADULT_CLOCKS), "DunedinPACE")
  outcomes <- outcomes[sapply(outcomes, function(o) any(!is.na(clocks_cmp[[o]])))]

  group_effects <- map_dfr(outcomes, function(outcome) {
    fit <- lm(as.formula(paste(outcome, "~", covars)), data = clocks_cmp)
    broom::tidy(fit, conf.int = TRUE) %>%
      filter(str_starts(term, "phenotype")) %>%
      mutate(outcome = outcome, n = nobs(fit), covariates = covars, .before = 1)
  })
  print(group_effects %>% select(outcome, n, term, estimate, std.error, p.value))
  write.csv(group_effects, CLOCK_GROUP_EFFECTS, row.names = FALSE)
  cat("Group comparison saved to:", CLOCK_GROUP_EFFECTS, "\n")
} else {
  cat("\nFewer than", CLOCK_MIN_GROUP_N, "samples in a group; skipping",
      group_label, "comparison\n")
  # Remove results from an earlier run so 3_plots.R doesn't draw stale ones
  if (file.exists(CLOCK_GROUP_EFFECTS)) file.remove(CLOCK_GROUP_EFFECTS)
}

cat("\nNext: Rscript clocks/3_plots.R", DATA_SOURCE, "\n")
