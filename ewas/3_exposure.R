# =============================================================================
# EWAS 3: exposures
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 16, 2026
# Description: CpG-wise limma models of the ComBat-corrected M values
#              (preprocessing/3_combat.R) on each exposure in EWAS_EXPOSURES
#              (ewas/config.R), from smoking questionnaires (load_smoking()),
#              the GP2 clinical data or the final sample sheet, adjusted for
#              EWAS_COVARIATES and phenotype. A continuous exposure gives one
#              result set (exposure_<name>); a categorical one gives one per
#              level vs its reference (exposure_<name>_<level>), all from the
#              same model. Saved in DIR_RESULTS_EWAS, with the per-sample
#              values (EWAS_EXPOSURE_VALUES). An exposure is skipped if it has
#              too little data
# Usage:       Rscript ewas/3_exposure.R <data source>   (see config.R)
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

library(tidyverse)
library(limma)

# Load shared configuration and EWAS helpers
source("config.R")
source("ewas/config.R")
source("R/ewas_functions.R")

# Exposure outputs from earlier runs (possibly of another exposure or level)
# would otherwise be picked up by 5_plots.R
remove_ewas_outputs("exposure_.*")

inp <- load_ewas_inputs()

# --- 1. Exposure values ------------------------------------------------------

# Joined by GP2ID; samples without a known value are NA and drop out of a model
smoking_cols <- c("smoking_source", "smoking_status", "smoking_ever")
if (any(EWAS_EXPOSURES$column %in% smoking_cols)) {
  cat("\nLoading smoking from", CLINICAL_DATA, "...\n")
  inp$targets <- inp$targets %>%
    select(-any_of(smoking_cols)) %>%
    left_join(load_smoking(inp$targets$GP2ID), by = "GP2ID")
  stopifnot(identical(inp$targets$Sample, colnames(inp$mVals)))

  cat("Smoking by source questionnaire (NA: unknown or not recorded):\n")
  print(table(inp$targets$smoking_source, inp$targets$smoking_ever, useNA = "ifany"))
  cat("Ever smoking by phenotype:\n")
  print(table(inp$targets$phenotype, inp$targets$smoking_ever, useNA = "ifany"))
  cat("Smoking status by phenotype (Former/Current only where current smoking is asked):\n")
  print(table(inp$targets$phenotype, inp$targets$smoking_status, useNA = "ifany"))
}

# Other exposures: CLINICAL_DATA columns, unless in the sample sheet already
clinical_items <- setdiff(EWAS_EXPOSURES$column, c(smoking_cols, names(inp$targets)))
if (length(clinical_items)) {
  cat("\nLoading", paste(clinical_items, collapse = ", "), "from", CLINICAL_DATA, "...\n")
  inp$targets <- inp$targets %>%
    left_join(load_clinical_items(inp$targets$GP2ID, clinical_items), by = "GP2ID")
  stopifnot(identical(inp$targets$Sample, colnames(inp$mVals)))
}

write.csv(inp$targets %>% select(Sample, GP2ID, phenotype, any_of(smoking_cols),
                                 all_of(setdiff(EWAS_EXPOSURES$column, smoking_cols))),
          EWAS_EXPOSURE_VALUES, row.names = FALSE)
cat("Exposure values saved to:", EWAS_EXPOSURE_VALUES, "\n")

# --- 2. Models ---------------------------------------------------------------

for (i in seq_len(nrow(EWAS_EXPOSURES))) {
  spec <- EWAS_EXPOSURES[i, ]
  fit_trait_model(inp, spec, inp$targets[[spec$column]], "exposure",
                  phenotype = EWAS_ADJUST_PHENOTYPE)
}

cat("\nNext: Rscript ewas/4_clinical.R", DATA_SOURCE, "\n")
