# =============================================================================
# EWAS 1: case-control
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 16, 2026
# Description: CpG-wise limma models of the ComBat-corrected M values
#              (preprocessing/3_combat.R) on diagnosis group, EWAS_GROUPS[2]
#              vs EWAS_GROUPS[1] (ewas/config.R), adjusted for EWAS_COVARIATES.
#              Where monogenic PD variants are known (R/pd_genetics.R, PPMI),
#              also the same comparison without carriers, and carriers of
#              each gene in EWAS_CARRIER_GENES vs non-carriers among cases.
#              Samples and cell proportions come from SAMPLE_SHEET_FINAL; CpGs
#              are annotated with CPG_ANNOTATION (preprocessing/4_cpg_annotation.R).
#              Saves <model>_{results.csv.gz,top_hits.csv,summary.csv} in
#              DIR_RESULTS_EWAS (models: case_control,
#              case_control_no_carriers, carrier_<gene>). A model is skipped
#              if a group has fewer than EWAS_MIN_GROUP_N samples (e.g.
#              Psomagen, no controls)
# Usage:       Rscript ewas/1_case_control.R <data source>   (see config.R)
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

library(tidyverse)
library(limma)

# Load shared configuration and EWAS helpers
source("config.R")
source("ewas/config.R")
source("R/ewas_functions.R")
source("R/pd_genetics.R")

# Outputs from earlier runs (possibly of other carrier genes) would otherwise
# be picked up by 4_plots.R
remove_ewas_outputs("(case_control|case_control_no_carriers|carrier_.*)")

inp <- load_ewas_inputs()

# Fit one two-group model: `group` is a factor (reference level first) with
# NA for samples left out. Skipped if a group is too small
fit_group_model <- function(prefix, label, group) {
  cat("\n---", label, "---\n")
  group_sizes <- table(group)
  print(group_sizes)
  if (any(group_sizes < EWAS_MIN_GROUP_N)) {
    cat("Fewer than", EWAS_MIN_GROUP_N, "samples in a group; skipping\n")
    return(invisible(NULL))
  }

  model_df <- tibble(group = group) %>%
    bind_cols(build_ewas_covariates(inp$targets))
  fit     <- fit_ewas(inp$mVals, model_df, "group")
  coef    <- paste0("group", levels(group)[2])
  results <- ewas_table(fit$fit, coef)

  # Group sizes among the samples actually in the model (complete data)
  n_model <- table(fit$model_df$group)
  groups  <- paste(names(n_model), n_model, collapse = " / ")
  write_ewas_result_set(results, prefix, inp$cpg_annotation,
                        summarize_ewas(prefix, label, coef, fit, results, groups))
}

# --- 1. Monogenic variant status ---------------------------------------------

# By clinical_id; NA where the data source has no variant data
inp$targets <- inp$targets %>% left_join(load_pd_gene(), by = "clinical_id")
stopifnot(identical(inp$targets$Sample, colnames(inp$mVals)))   # one variant row per participant
has_variants <- any(!is.na(inp$targets$pd_gene))
if (has_variants) {
  cat("\nMonogenic PD variants by phenotype (NA: no variant data):\n")
  print(table(inp$targets$phenotype, inp$targets$pd_gene, useNA = "ifany"))
}

# --- 2. Case-control ---------------------------------------------------------

# Phenotypes outside EWAS_GROUPS become NA and are dropped by fit_ewas()
group_label <- paste(rev(EWAS_GROUPS), collapse = " vs ")
diagnosis   <- factor(inp$targets$phenotype, levels = EWAS_GROUPS)
fit_group_model("case_control", group_label, diagnosis)

# --- 3. Case-control without carriers ----------------------------------------

if (has_variants) {
  carrier <- !is.na(inp$targets$pd_gene) & inp$targets$pd_gene != EWAS_NONCARRIER
  no_carriers <- replace(diagnosis, carrier, NA)
  fit_group_model("case_control_no_carriers",
                  paste0(group_label, ", without monogenic variant carriers"), no_carriers)
}

# --- 4. Carriers vs non-carriers among cases ---------------------------------

if (has_variants) {
  for (gene in EWAS_CARRIER_GENES) {
    in_model <- inp$targets$phenotype %in% EWAS_GROUPS[2] &
                inp$targets$pd_gene %in% c(EWAS_NONCARRIER, gene)
    carrier_group <- factor(ifelse(in_model, inp$targets$pd_gene, NA),
                            levels = c(EWAS_NONCARRIER, gene))
    fit_group_model(paste0("carrier_", gene),
                    paste0(gene, " vs non-carrier ", EWAS_GROUPS[2]), carrier_group)
  }
} else {
  cat("\nNo monogenic variant data for", DATA_SOURCE,
      "(PPMI_GENETICS in config.R); skipping the carrier models\n")
}

cat("\nNext: Rscript ewas/2_age.R", DATA_SOURCE, "\n")
