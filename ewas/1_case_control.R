# =============================================================================
# EWAS 1: case-control
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 16, 2026
# Description: CpG-wise limma model of the ComBat-corrected M values
#              (preprocessing/3_combat.R) on diagnosis group, EWAS_GROUPS[2]
#              vs EWAS_GROUPS[1] (ewas/config.R), adjusted for EWAS_COVARIATES.
#              Samples and cell proportions come from SAMPLE_SHEET_FINAL; CpGs
#              are annotated with CPG_ANNOTATION (preprocessing/4_cpg_annotation.R).
#              Saves case_control_{results.csv.gz,top_hits.csv,summary.csv}
#              in DIR_RESULTS_EWAS. Skipped if a group has fewer than
#              EWAS_MIN_GROUP_N samples (e.g. Psomagen, no controls)
# Usage:       Rscript ewas/1_case_control.R <data source>   (see config.R)
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

library(tidyverse)
library(limma)

# Load shared configuration and EWAS helpers
source("config.R")
source("ewas/config.R")
source("R/ewas_functions.R")

PREFIX <- "case_control"

inp <- load_ewas_inputs()

# --- 1. Samples in the comparison --------------------------------------------

# Phenotypes outside EWAS_GROUPS become NA and are dropped by fit_ewas()
group       <- factor(inp$targets$phenotype, levels = EWAS_GROUPS)
group_sizes <- table(group)
group_label <- paste(rev(EWAS_GROUPS), collapse = " vs ")
cat("\nSamples in", group_label, "comparison:\n")
print(group_sizes)

if (any(group_sizes < EWAS_MIN_GROUP_N)) {
  cat("\nFewer than", EWAS_MIN_GROUP_N, "samples in a group; skipping the",
      group_label, "EWAS\n")
  remove_ewas_outputs(PREFIX)
  quit(save = "no", status = 0)
}

# --- 2. Model ----------------------------------------------------------------

model_df <- tibble(phenotype = group) %>%
  bind_cols(build_ewas_covariates(inp$targets))

cat("\nFitting", group_label, "EWAS...\n")
fit  <- fit_ewas(inp$mVals, model_df, "phenotype")
coef <- paste0("phenotype", EWAS_GROUPS[2])
results <- ewas_table(fit$fit, coef)

# Group sizes among the samples actually in the model (complete data)
n_model <- table(fit$model_df$phenotype)
groups  <- paste(names(n_model), n_model, collapse = " / ")

# --- 3. Save -----------------------------------------------------------------

summary_row <- summarize_ewas(PREFIX, group_label, coef, fit, results, groups)
write_ewas_result_set(results, PREFIX, inp$cpg_annotation, summary_row)

cat("\nNext: Rscript ewas/2_age.R", DATA_SOURCE, "\n")
