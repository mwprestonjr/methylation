# =============================================================================
# EWAS functions
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 16, 2026
# Description: Shared helpers for the ewas/ scripts: loading the preprocessing
#              outputs (ComBat-corrected M values, final sample sheet, CpG
#              annotation) and smoking status (GP2 clinical data), building
#              the model covariates, fitting limma
#              CpG-wise models, and writing each model's result set and summary
# =============================================================================

library(limma)
library(tidyverse)
library(data.table)   # fread/fwrite for the per-CpG tables

# --- Inputs ------------------------------------------------------------------

# ComBat-corrected M values (COMBAT_MVALS), the final sample sheet
# (SAMPLE_SHEET_FINAL: QC-passed samples + cell proportions + combat_batch) in
# M value column order, and the CpG annotation (CPG_ANNOTATION), all from
# preprocessing/
load_ewas_inputs <- function() {
  for (f in c(COMBAT_MVALS, SAMPLE_SHEET_FINAL, CPG_ANNOTATION)) {
    if (!file.exists(f)) stop(f, " not found - run make preprocessing SOURCE=", DATA_SOURCE)
  }

  cat("Loading ComBat-corrected M values...\n")
  mVals <- readRDS(COMBAT_MVALS)
  cat("M values dimensions:", dim(mVals), "\n")

  # combat_batch and the Sentrix fields are labels, not numbers
  targets <- read.csv(SAMPLE_SHEET_FINAL,
                      colClasses = c(GP2ID            = "character",
                                     clinical_id      = "character",
                                     Sentrix_ID       = "character",
                                     Sentrix_Position = "character",
                                     Basename         = "character",
                                     combat_batch     = "character")) %>%
    mutate(Sample   = basename(Basename),
           chip_row = substr(Sentrix_Position, 1, 3))

  # Sample sheet in M value column order (columns are idat basenames)
  targets <- targets[match(colnames(mVals), targets$Sample), ]
  stopifnot(all(colnames(mVals) == targets$Sample))
  cat("Samples:", nrow(targets), "\n")
  cat("Phenotypes:\n")
  print(table(targets$phenotype, useNA = "ifany"))

  # One row per probe on the array: cpg, chr (no "chr" prefix), pos_hg38,
  # island, gene annotation (preprocessing/4_cpg_annotation.R)
  cpg_annotation <- fread(CPG_ANNOTATION, colClasses = c(chr = "character"))
  cat("CpG annotation:", nrow(cpg_annotation), "probes,",
      sum(rownames(mVals) %in% cpg_annotation$cpg), "of", nrow(mVals), "M value CpGs annotated\n")

  list(mVals = mVals, targets = targets, cpg_annotation = cpg_annotation)
}

# Smoking at the blood draw from the GP2 extended clinical data
# (CLINICAL_DATA), one row per GP2ID with a known status: the harmonized
# smoking_status at the first of EWAS_CLINICAL_VISIT_MONTHS with one, as
#   smoking_status: Never, Former, Current
#   smoking_ever:   Never, Ever
# "Unknown" and missing values are left out (NA after joining)
SMOKING_LEVELS <- c("Never" = "Never", "Former smoker" = "Former", "Current smoker" = "Current")

load_smoking <- function(gp2ids) {
  if (!file.exists(CLINICAL_DATA)) stop(CLINICAL_DATA, " not found (CLINICAL_DATA in config.R)")
  readr::read_csv(CLINICAL_DATA,
                  col_types  = readr::cols(.default = "c"),
                  col_select = c("GP2ID", "visit_month", "smoking_status"),
                  progress   = FALSE) %>%
    filter(GP2ID %in% gp2ids, smoking_status %in% names(SMOKING_LEVELS)) %>%
    mutate(visit_rank = match(as.numeric(visit_month), EWAS_CLINICAL_VISIT_MONTHS)) %>%
    filter(!is.na(visit_rank)) %>%
    arrange(GP2ID, visit_rank) %>%
    distinct(GP2ID, .keep_all = TRUE) %>%
    transmute(GP2ID,
              smoking_status = unname(SMOKING_LEVELS[smoking_status]),
              smoking_ever   = ifelse(smoking_status == "Never", "Never", "Ever"))
}

# --- Model -------------------------------------------------------------------

# Covariate table for the model, one row per sample in `targets` order:
# EWAS_COVARIATES (minus `exclude`, e.g. age in the age model), plus phenotype
# if `phenotype` is TRUE. Categorical covariates are factors
build_ewas_covariates <- function(targets, exclude = character(0), phenotype = FALSE) {
  covariates <- setdiff(c(EWAS_COVARIATES, if (phenotype) "phenotype"), exclude)
  missing <- setdiff(covariates, names(targets))
  if (length(missing)) {
    stop("Covariates not in SAMPLE_SHEET_FINAL: ", paste(missing, collapse = ", "),
         " (EWAS_COVARIATES in ewas/config.R)")
  }
  categorical <- c("sex", "chip_row", "phenotype", "combat_batch", "race")
  as_tibble(lapply(targets[covariates], function(x) {
    if (is.character(x) || is.logical(x)) factor(x) else x
  })) %>%
    mutate(across(any_of(categorical), factor))
}

# Fit M value ~ all columns of `model_df` for every CpG (limma lmFit + eBayes).
# `variable` is the column of interest; the other columns are covariates.
# Samples with a missing value in model_df are dropped; then categorical
# covariates with a single value among the remaining samples are left out and
# unused factor levels dropped, so every design column can be estimated
fit_ewas <- function(mVals, model_df, variable) {
  keep <- complete.cases(model_df)
  model_df <- droplevels(model_df[keep, , drop = FALSE])
  if (nrow(model_df) < EWAS_MIN_SAMPLES) {
    stop("Too few samples with complete data for the EWAS: ", nrow(model_df),
         " (EWAS_MIN_SAMPLES = ", EWAS_MIN_SAMPLES, ")")
  }

  single <- names(model_df)[sapply(model_df, function(x)
    !is.numeric(x) && length(unique(x)) < 2)]
  if (variable %in% single) stop(variable, " has a single value among the samples")
  if (length(single)) {
    cat("Covariates left out (one value only):", single, "\n")
    model_df <- model_df[, setdiff(names(model_df), single), drop = FALSE]
  }

  design <- model.matrix(~ ., data = model_df)
  if (qr(design)$rank < ncol(design)) {
    stop("Design matrix is not full rank; check for collinear covariates: ",
         paste(colnames(design), collapse = ", "))
  }
  cat("Model: M value ~", paste(names(model_df), collapse = " + "), "\n")
  cat("Samples in model:", nrow(model_df), "| CpGs:", nrow(mVals), "\n")

  fit <- eBayes(lmFit(mVals[, keep, drop = FALSE], design))
  list(fit = fit, design = design, model_df = model_df,
       n_samples = nrow(model_df),
       covariates = setdiff(names(model_df), variable))
}

# One coefficient's per-CpG results, sorted by p: logFC (effect on the M value
# scale per unit or vs the reference level), AveExpr (mean M value), t, P.Value,
# FDR (BH across CpGs), B
ewas_table <- function(fit, coef) {
  stopifnot(coef %in% colnames(fit$coefficients))
  topTable(fit, coef = coef, number = Inf, sort.by = "P") %>%
    rownames_to_column("cpg") %>%
    rename(FDR = adj.P.Val)
}

# Genomic inflation factor from p-values (median chi-square / its expected
# value under the null, 0.456); ~1 when test statistics aren't inflated
genomic_lambda <- function(p) {
  median(qchisq(p, df = 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, df = 1)
}

# One-row summary of a model's results, collected by 4_plots.R into
# EWAS_MODEL_SUMMARY
summarize_ewas <- function(model, label, coefficient, fit, results, groups = NA_character_) {
  tibble(model       = model,
         label       = label,
         coefficient = coefficient,
         n_samples   = fit$n_samples,
         groups      = groups,
         covariates  = paste(fit$covariates, collapse = ", "),
         n_cpg       = nrow(results),
         n_fdr_sig   = sum(results$FDR < EWAS_FDR, na.rm = TRUE),
         fdr         = EWAS_FDR,
         min_p       = min(results$P.Value, na.rm = TRUE),
         lambda      = genomic_lambda(results$P.Value))
}

# --- Outputs -----------------------------------------------------------------

# A model's output files in DIR_RESULTS_EWAS (see Outputs in ewas/config.R)
ewas_output_files <- function(prefix) {
  list(results  = file.path(DIR_RESULTS_EWAS, paste0(prefix, "_results.csv.gz")),
       top_hits = file.path(DIR_RESULTS_EWAS, paste0(prefix, "_top_hits.csv")),
       summary  = file.path(DIR_RESULTS_EWAS, paste0(prefix, "_summary.csv")))
}

# Remove a model's outputs (and figures) from an earlier run, so a skipped
# model doesn't leave stale results for 4_plots.R. `prefix` is a regex matched
# at the start of the file name
remove_ewas_outputs <- function(prefix) {
  old <- c(list.files(DIR_RESULTS_EWAS, pattern = paste0("^", prefix, "_(results\\.csv\\.gz|top_hits\\.csv|summary\\.csv)$"),
                      full.names = TRUE),
           list.files(DIR_FIGURES_EWAS, pattern = paste0("^", prefix, "_.*\\.png$"), full.names = TRUE))
  if (length(old)) {
    cat("Removing outputs from an earlier run:", basename(old), "\n")
    file.remove(old)
  }
  invisible(old)
}

# Write a model's per-CpG results merged with the CpG annotation, its top hits
# and its summary row; returns the file paths
write_ewas_result_set <- function(results, prefix, cpg_annotation, summary_row) {
  files  <- ewas_output_files(prefix)
  merged <- left_join(results, cpg_annotation, by = "cpg")
  stopifnot(nrow(merged) == nrow(results))   # one annotation row per CpG

  fwrite(merged, files$results)
  write.csv(head(merged, EWAS_TOP_HITS), files$top_hits, row.names = FALSE)   # already sorted by p
  write.csv(summary_row, files$summary, row.names = FALSE)

  cat("CpGs with FDR <", EWAS_FDR, ":", summary_row$n_fdr_sig, "of", summary_row$n_cpg,
      "| min p:", signif(summary_row$min_p, 3), "| lambda:", round(summary_row$lambda, 3), "\n")
  cat("Results saved to:", files$results, "\n")
  cat("Top hits saved to:", files$top_hits, "\n")
  cat("Summary saved to:", files$summary, "\n")
  invisible(files)
}
