# =============================================================================
# EWAS 2: age
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 16, 2026
# Description: CpG-wise limma model of the ComBat-corrected M values
#              (preprocessing/3_combat.R) on age at the blood draw
#              (EWAS_AGE_COL), across all samples, adjusted for the other
#              EWAS_COVARIATES and phenotype (EWAS_ADJUST_PHENOTYPE in
#              ewas/config.R). Effects are per year. Saves
#              age_{results.csv.gz,top_hits.csv,summary.csv} in
#              DIR_RESULTS_EWAS. Age wasn't protected in ComBat, so batches
#              that differ in age lose some of their age signal
# Usage:       Rscript ewas/2_age.R <data source>   (see config.R)
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

library(tidyverse)
library(limma)

# Load shared configuration and EWAS helpers
source("config.R")
source("ewas/config.R")
source("R/ewas_functions.R")

PREFIX <- "age"

inp <- load_ewas_inputs()

# --- 1. Model ----------------------------------------------------------------

# Age is the variable of interest, so it's taken out of the covariates
age <- suppressWarnings(as.numeric(inp$targets[[EWAS_AGE_COL]]))
cat("\nAge:", sum(!is.na(age)), "samples with a value, mean", round(mean(age, na.rm = TRUE), 1),
    "range", paste(range(age, na.rm = TRUE), collapse = "-"), "\n")

model_df <- tibble(age = age) %>%
  bind_cols(build_ewas_covariates(inp$targets, exclude = "age",
                                  phenotype = EWAS_ADJUST_PHENOTYPE))

cat("\nFitting age EWAS...\n")
fit     <- fit_ewas(inp$mVals, model_df, "age")
results <- ewas_table(fit$fit, "age")

# --- 2. Save -----------------------------------------------------------------

summary_row <- summarize_ewas(PREFIX, "Age at sample collection (per year)", "age", fit, results)
write_ewas_result_set(results, PREFIX, inp$cpg_annotation, summary_row)

cat("\nNext: Rscript ewas/3_exposure.R", DATA_SOURCE, "\n")
