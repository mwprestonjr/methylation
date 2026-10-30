# =============================================================================
# EWAS 3: exposure
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 16, 2026
# Description: CpG-wise limma model of the ComBat-corrected M values
#              (preprocessing/3_combat.R) on an exposure in the final sample
#              sheet or smoking from the GP2 clinical data (EWAS_EXPOSURE_COL:
#              ever/never regular smoking for both data sources), adjusted
#              for EWAS_COVARIATES and phenotype
#              (ewas/config.R). Saves the per-sample smoking values
#              (EWAS_SMOKING_VALUES). A continuous
#              exposure gives one result set (exposure_<column>); a
#              categorical one gives one per level vs EWAS_EXPOSURE_REFERENCE
#              (exposure_<column>_<level>), all from the same model. Saved in
#              DIR_RESULTS_EWAS. Skipped if the exposure isn't available or
#              has too little data
# Usage:       Rscript ewas/3_exposure.R <data source>   (see config.R)
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

library(tidyverse)
library(limma)

# Load shared configuration and EWAS helpers
source("config.R")
source("ewas/config.R")
source("R/ewas_functions.R")

PREFIX <- paste0("exposure_", EWAS_EXPOSURE_COL)

# Exposure outputs from earlier runs (possibly of another exposure column or
# level) would otherwise be picked up by 4_plots.R
remove_ewas_outputs("exposure_.*")

inp <- load_ewas_inputs()

skip <- function(reason) {
  cat("\nSkipping the exposure EWAS:", reason, "\n")
  quit(save = "no", status = 0)
}

# --- 1. Smoking from the clinical data ---------------------------------------

# Joined by GP2ID (replicate samples of a participant get the same value);
# samples without a known status are NA and drop out of the model
smoking_cols <- c("smoking_source", "smoking_status", "smoking_ever")
if (EWAS_EXPOSURE_COL %in% smoking_cols) {
  cat("\nLoading smoking from", CLINICAL_DATA, "...\n")
  inp$targets <- inp$targets %>%
    select(-any_of(smoking_cols)) %>%
    left_join(load_smoking(inp$targets$GP2ID), by = "GP2ID")
  stopifnot(all(inp$targets$Sample == colnames(inp$mVals)))

  cat("Smoking by source questionnaire (NA: unknown or not recorded):\n")
  print(table(inp$targets$smoking_source, inp$targets$smoking_ever, useNA = "ifany"))
  cat("Ever smoking by phenotype:\n")
  print(table(inp$targets$phenotype, inp$targets$smoking_ever, useNA = "ifany"))
  cat("Smoking status by phenotype (Former/Current only where current smoking is asked):\n")
  print(table(inp$targets$phenotype, inp$targets$smoking_status, useNA = "ifany"))
  write.csv(inp$targets %>% select(Sample, GP2ID, phenotype, all_of(smoking_cols)),
            EWAS_SMOKING_VALUES, row.names = FALSE)
  cat("Smoking values saved to:", EWAS_SMOKING_VALUES, "\n")
}

# --- 2. Exposure -------------------------------------------------------------

if (!EWAS_EXPOSURE_COL %in% names(inp$targets)) {
  skip(paste0("'", EWAS_EXPOSURE_COL, "' is not a column of ", SAMPLE_SHEET_FINAL,
              " (EWAS_EXPOSURE_COL in ewas/config.R)"))
}
exposure_raw <- inp$targets[[EWAS_EXPOSURE_COL]]
exposure_raw[exposure_raw %in% ""] <- NA

if (EWAS_EXPOSURE_TYPE == "categorical") {
  counts <- table(exposure_raw)
  cat("\n", EWAS_EXPOSURE_COL, ":\n", sep = "")
  print(counts)
  # Levels too small to estimate are left out (NA); reference level first
  levels <- names(counts)[counts >= EWAS_MIN_GROUP_N]
  if (!EWAS_EXPOSURE_REFERENCE %in% levels) {
    skip(paste0("reference level '", EWAS_EXPOSURE_REFERENCE, "' missing or has fewer than ",
                EWAS_MIN_GROUP_N, " samples"))
  }
  if (length(levels) < 2) skip("fewer than 2 levels with enough samples")
  if (length(setdiff(names(counts), levels))) {
    cat("Levels left out (fewer than", EWAS_MIN_GROUP_N, "samples):",
        setdiff(names(counts), levels), "\n")
  }
  exposure <- factor(exposure_raw, levels = c(EWAS_EXPOSURE_REFERENCE,
                                              setdiff(levels, EWAS_EXPOSURE_REFERENCE)))
} else if (EWAS_EXPOSURE_TYPE == "continuous") {
  exposure <- suppressWarnings(as.numeric(exposure_raw))
  cat("\n", EWAS_EXPOSURE_COL, ": ", sum(!is.na(exposure)), " samples with a value\n", sep = "")
  if (sd(exposure, na.rm = TRUE) %in% c(0, NA)) skip("no variation")
} else {
  stop("EWAS_EXPOSURE_TYPE must be 'continuous' or 'categorical'")
}

# --- 3. Model ----------------------------------------------------------------

model_df <- tibble(exposure = exposure) %>%
  bind_cols(build_ewas_covariates(inp$targets, phenotype = EWAS_ADJUST_PHENOTYPE))

cat("\nFitting", EWAS_EXPOSURE_COL, "EWAS...\n")
fit <- fit_ewas(inp$mVals, model_df, "exposure")

# --- 4. Save, one result set per coefficient ---------------------------------

if (EWAS_EXPOSURE_TYPE == "categorical") {
  n_model <- table(fit$model_df$exposure)
  groups  <- paste(names(n_model), n_model, collapse = " / ")
  for (lvl in levels(fit$model_df$exposure)[-1]) {
    coef    <- paste0("exposure", lvl)
    prefix  <- paste0(PREFIX, "_", make.names(lvl))
    label   <- paste0(EWAS_EXPOSURE_COL, ": ", lvl, " vs ", EWAS_EXPOSURE_REFERENCE)
    cat("\n", label, "\n", sep = "")
    results <- ewas_table(fit$fit, coef)
    write_ewas_result_set(results, prefix, inp$cpg_annotation,
                          summarize_ewas(prefix, label, coef, fit, results, groups))
  }
} else {
  results <- ewas_table(fit$fit, "exposure")
  write_ewas_result_set(results, PREFIX, inp$cpg_annotation,
                        summarize_ewas(PREFIX, paste(EWAS_EXPOSURE_COL, "(per unit)"),
                                       "exposure", fit, results))
}

cat("\nNext: Rscript ewas/4_plots.R", DATA_SOURCE, "\n")
