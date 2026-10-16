# =============================================================================
# Clocks 2: compare age acceleration between groups
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 8, 2026
# Updated: Oct 13, 2026
# Description: Compares age acceleration and DunedinPACE between CLOCK_GROUPS
#              (first = reference), adjusted for CLOCK_COVARIATES
#              (clocks/config.R), using the per-sample table from
#              1_estimate.R, with FDR (BH) across the outcomes. Saves the
#              group effect per outcome for 4_plots.R. Repeats the comparison
#              without the cell proportions (CLOCK_COVARIATES_NO_CELLS) as a
#              sensitivity analysis, saved separately.
#              Skipped if a group has fewer than CLOCK_MIN_GROUP_N samples
#              (e.g. a data source without controls)
# Usage:       Rscript clocks/2_compare.R <data source>   (see config.R)
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------
library(tidyverse)

# Load shared configuration
source("config.R")
source("clocks/config.R")
source("R/clock_models.R")

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

# Main model (CLOCK_COVARIATES) and the sensitivity model without cell
# proportions (CLOCK_COVARIATES_NO_CELLS), each saved to its own file
models <- list(
  list(label = "main model",               covariates = CLOCK_COVARIATES,          file = CLOCK_GROUP_EFFECTS),
  list(label = "without cell proportions", covariates = CLOCK_COVARIATES_NO_CELLS, file = CLOCK_GROUP_EFFECTS_NO_CELLS)
)

if (all(group_sizes >= CLOCK_MIN_GROUP_N)) {
  # Outcomes not estimated on this array (all NA, below the coverage
  # threshold) are left out
  outcomes <- estimated_outcomes(clocks_cmp, c(paste0("AgeAccel_", ADULT_CLOCKS), "DunedinPACE"))

  for (m in models) {
    # Every outcome gets the same covariates, including age: the AgeAccel
    # residuals are uncorrelated with age across all samples, but not
    # necessarily within this subset, and adjusting in the same model handles
    # an age difference between the groups correctly. Categorical covariates
    # with a single value here can't be estimated and are left out
    usable <- usable_covariates(clocks_cmp, m$covariates)
    covars <- paste(c("phenotype", usable), collapse = " + ")
    cat("\nAge acceleration,", group_label, paste0("(", m$label, ")"),
        "- model: <outcome> ~", covars, "\n")
    if (length(setdiff(m$covariates, usable))) {
      cat("Covariates left out (one value only):", setdiff(m$covariates, usable), "\n")
    }

    # FDR across all outcomes (and group terms, if CLOCK_GROUPS has more than 2)
    group_effects <- fit_predictor(clocks_cmp, outcomes, "phenotype", usable) %>%
      mutate(fdr = p.adjust(p.value, method = "BH"))
    print(group_effects %>% select(outcome, n, term, estimate, std.error, p.value, fdr))
    cat("FDR <", CLOCK_FDR, ":", sum(group_effects$fdr < CLOCK_FDR), "of", nrow(group_effects), "\n")
    write.csv(group_effects, m$file, row.names = FALSE)
    cat("Group comparison saved to:", m$file, "\n")
  }
} else {
  cat("\nFewer than", CLOCK_MIN_GROUP_N, "samples in a group; skipping",
      group_label, "comparison\n")
  # Remove results from an earlier run so 4_plots.R doesn't draw stale ones
  for (m in models) if (file.exists(m$file)) file.remove(m$file)
}

cat("\nNext: Rscript clocks/3_associations.R", DATA_SOURCE, "\n")
