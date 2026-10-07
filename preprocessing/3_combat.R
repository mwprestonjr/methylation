# =============================================================================
# Preprocessing 3: ComBat batch correction
# Author: GP2 Subtypes and Mechanisms - M.E., M.P.
# Date: May 6, 2026
# Updated: Sept 30, 2026
# Description: Removes batch effects (COMBAT_BATCH_VAR) from the QC-passed M
#              values with ComBat, keeping the variation explained by
#              COMBAT_PROTECT. Saves the corrected M values and the final
#              sample sheet (QC-passed sheet + cell proportions) used by later
#              analyses (e.g. mQTL). Check the result with
#              2_variation_sources.R <data source> combat
# Usage:       Rscript preprocessing/3_combat.R <data source>   (see config.R)
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

library(tidyverse)
library(sva)                      # ComBat

# Load shared configuration and M value helpers
source("config.R")
source("preprocessing/config.R")
source("R/mvalues.R")

# --- 1. Load data ------------------------------------------------------------

cat("Loading QC-passed M values...\n")
# ComBat can't handle infinite values; files from before Script 02 capped
# them may still contain some (R/mvalues.R)
mVals <- cap_infinite_m(readRDS(file.path(DIR_RESULTS, "mVals.rds")))
cat("M values dimensions:", dim(mVals), "\n")

targets <- read.csv(SAMPLE_SHEET_QC,
                    colClasses = c(GP2ID            = "character",
                                   clinical_id      = "character",
                                   Sentrix_ID       = "character",
                                   Sentrix_Position = "character",
                                   Basename         = "character")) %>%
  mutate(Sample = basename(Basename)) %>%
  left_join(read.csv(CELL_PROPORTIONS), by = "Sample")

# Sample sheet in M value column order
targets <- targets[match(colnames(mVals), targets$Sample), ]
stopifnot(all(targets$Sample == colnames(mVals)))
cat("Samples:", nrow(targets), "\n")

# --- 2. Batch and protected variables ----------------------------------------

stopifnot(COMBAT_BATCH_VAR %in% names(targets), all(COMBAT_PROTECT %in% names(targets)))
batch <- as.character(targets[[COMBAT_BATCH_VAR]])
cat("\nBatch variable:", COMBAT_BATCH_VAR, "-", length(unique(batch)), "batches\n")
batch_sizes <- table(batch)
cat("Batch sizes: min", min(batch_sizes), "/ median", median(batch_sizes),
    "/ max", max(batch_sizes), "\n")

# A batch of one sample can't be corrected: its batch mean is the sample's
# own profile, so correcting it would remove the sample's own signal (e.g.
# genetic effects). Use a coarser batch variable (e.g. plate instead of chip)
if (any(batch_sizes < 2)) {
  stop(sum(batch_sizes < 2), " of ", length(batch_sizes), " batches of ",
       COMBAT_BATCH_VAR, " have a single sample - choose a batch variable with ",
       "larger batches (COMBAT_BATCH_VAR in preprocessing/config.R)")
}
if (median(batch_sizes) < 4) {
  cat("WARNING: median batch size is", median(batch_sizes),
      "- ComBat estimates are unstable for very small batches\n")
}

# ComBat keeps the variation explained by the protected variables. Samples
# with missing values can't be in the model
protect_cols <- COMBAT_PROTECT[sapply(COMBAT_PROTECT, function(v)
  length(unique(na.omit(targets[[v]]))) > 1)]
cat("Protected variables:", if (length(protect_cols)) paste(protect_cols, collapse = ", ") else "none", "\n")
missing <- !complete.cases(targets[, protect_cols, drop = FALSE])
if (any(missing)) {
  stop(sum(missing), " samples have missing values for COMBAT_PROTECT (",
       paste(protect_cols, collapse = ", "), "): ",
       paste(targets$Sample[missing], collapse = ", "))
}
mod <- if (length(protect_cols)) {
  model.matrix(reformulate(protect_cols), data = targets)
} else {
  NULL
}

# If batch and phenotype overlap strongly, protecting phenotype can inflate
# phenotype differences (Nygaard et al. 2016); see confounder_summary.csv
# from 2_variation_sources.R
if ("phenotype" %in% protect_cols) {
  t <- suppressWarnings(chisq.test(table(targets$phenotype, batch)))
  cat("Batch vs phenotype (chi-square): p =", format(t$p.value, digits = 3), "\n")
}

# --- 3. ComBat ---------------------------------------------------------------

cat("\nRunning ComBat...\n")
combat_mVals <- ComBat(dat         = as.matrix(mVals),
                       batch       = batch,
                       mod         = mod,
                       par.prior   = TRUE,
                       prior.plots = FALSE)
cat("Corrected M values dimensions:", dim(combat_mVals), "\n")
rm(mVals); invisible(gc())

# --- 4. Save outputs ---------------------------------------------------------

saveRDS(combat_mVals, COMBAT_MVALS)
cat("\nBatch-corrected M values saved to:", COMBAT_MVALS, "\n")

# Final sample sheet: QC-passed sheet + cell proportions, plus the batch used
targets <- targets %>% mutate(combat_batch = batch)
write.csv(targets, SAMPLE_SHEET_FINAL, row.names = FALSE)
cat("Final sample sheet saved to:", SAMPLE_SHEET_FINAL, "\n")

combat_summary <- data.frame(
  metric = c("Samples", "CpGs", "Batch variable", "Batches",
             "Smallest batch", "Protected variables"),
  value  = c(ncol(combat_mVals), nrow(combat_mVals), COMBAT_BATCH_VAR,
             length(batch_sizes), min(batch_sizes),
             paste(protect_cols, collapse = ", "))
)
write.csv(combat_summary, file.path(DIR_RESULTS, "combat_summary.csv"), row.names = FALSE)

cat("\nNext: Rscript preprocessing/2_variation_sources.R", DATA_SOURCE, "combat\n")
