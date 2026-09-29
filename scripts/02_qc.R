# =============================================================================
# Quality Control
# Author: GP2 Subtypes and Mechanisms - M.E., M.P.
# Date: April 29, 2026
# Updated: Sept 29, 2026
# Description: Reads the idat files in the sample sheet from Script 01 into
#              minfi, runs QC following the minfi user guide, optionally
#              computes SeSAMe per-sample QC stats (USE_SESAME_QC), filters
#              failed samples and probes, normalizes using preprocessFunnorm,
#              and saves QC-passed data plus the data for the QC figures
#              (drawn by 03_qc_plots.R). Works for any array in
#              R/array_profiles.R; one array type per run
# Usage:       Rscript scripts/02_qc.R <data source>   (see config.R)
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------
library(minfi)
library(tidyverse)
library(maxprobes)

# Load shared configuration and array profiles
source("config.R")
source("R/array_profiles.R")

# --- 0. Load sample sheet ----------------------------------------------------

cat("Loading sample sheet...\n")
targets <- read.csv(SAMPLE_SHEET,
                    colClasses = c(Sentrix_ID = "character",
                                   Sentrix_Position = "character",
                                   clinical_id = "character",
                                   Basename = "character"))
# Quick checks
cat("Samples loaded:", nrow(targets), "\n")
cat("Samples per dataset and array:\n")
print(table(targets$Dataset, targets$Array))
cat("Phenotype breakdown:\n")
print(table(targets$phenotype, useNA = "ifany"))
cat("Female:", sum(targets$sex == "Female", na.rm = TRUE), "\n")
cat("Male:", sum(targets$sex == "Male", na.rm = TRUE), "\n")
cat("\nSex by phenotype breakdown:\n")
print(table(targets$phenotype, targets$sex))

# minfi can't read different array types into one RGChannelSet, and the
# annotation and cross-reactive probes depend on the array
array_type <- unique(targets$Array)
if (length(array_type) != 1 || !array_type %in% names(ARRAY_PROFILES)) {
  stop("Sample sheet must contain exactly one supported array type (",
       paste(names(ARRAY_PROFILES), collapse = ", "), "); found: ",
       paste(array_type, collapse = ", "), " - see the array table from Script 01")
}
profile <- ARRAY_PROFILES[[array_type]]
library(profile$anno_pkg, character.only = TRUE)
cat("Array:", array_type, "- annotation:", profile$anno_pkg, "\n")

# --- 1. SeSAMe QC stats ---------------------------------------------------

# Per-sample detection (pOOBAH), intensity, dye bias and beta distribution.
# Only runs when USE_SESAME_QC is TRUE (see config.R), and then its
# detection rate is used for sample removal in 3b.
# NOTE: sesame's bisulfite conversion score (bisConversionControl) fails on
# EPICv2 in sesame 1.24, so it is not included
if (USE_SESAME_QC) {
  cat("\nComputing SeSAMe QC stats \n")
  # The idats are on a gcsfuse mount, where reads occasionally fail with
  # "error reading from connection" under parallel load; retry those reads
  read_idat_pair <- function(b, attempts = 3) {
    for (i in seq_len(attempts)) {
      sdf <- try(sesame::readIDATpair(b), silent = TRUE)
      if (!inherits(sdf, "try-error")) return(sdf)
      Sys.sleep(5 * i)
    }
    stop(sprintf("reading %s failed after %d attempts: %s", b, attempts, sdf))
  }

  # mc.preschedule = FALSE runs each sample as its own job, so one error
  # only affects that sample instead of every sample on the same core
  sesame_qc <- parallel::mclapply(targets$Basename, function(b) {
    sdf <- read_idat_pair(b)
    as.data.frame(sesame::sesameQC_getStats(sesame::sesameQC_calcStats(sdf)))
  }, mc.cores = max(1, parallel::detectCores() - 1), mc.preschedule = FALSE)

  # mclapply returns errors instead of stopping; report which samples failed
  sesame_failed <- map_lgl(sesame_qc, inherits, "try-error")
  if (any(sesame_failed)) {
    stop("SeSAMe QC failed for: ", paste(targets$Basename[sesame_failed], collapse = ", "),
         "\n", paste(unique(unlist(sesame_qc[sesame_failed])), collapse = "\n"))
  }

  sesame_qc <- bind_cols(
    targets %>% select(GP2ID, GP2sampleID, clinical_id, Dataset, Sentrix_ID, Sentrix_Position),
    bind_rows(sesame_qc)
  )
  write.csv(sesame_qc, file.path(DIR_RESULTS, "qc_sesame_stats.csv"), row.names = FALSE)
  cat("SeSAMe QC stats saved\n")
} else {
  cat("\nSkipping SeSAMe QC (USE_SESAME_QC = FALSE); using minfi detection p-values\n")
}

# --- 2. Read IDAT files into minfi ------------------------------------------------------

cat("\nReading IDAT files into minfi...\n")
rgSet <- read.metharray.exp(targets = targets,
                        verbose   = TRUE,
                        force     = TRUE,
                        extended = TRUE)

# annotate the RGChannelSet with the array and annotation version from the
# profile (minfi's own guess is wrong for EPICv2)
annotation(rgSet) <- profile$annotation

# Verify metadata attached correctly
cat("RGChannelSet dimensions:", dim(rgSet), "\n")
cat("Array type:", annotation(rgSet)["array"], "\n")
cat("Sample metadata columns:", ncol(pData(rgSet)), "\n")
cat("First few GP2IDs:", head(pData(rgSet)$GP2ID), "\n")

# --- 3. Initial QC -----------------------------------------------------------

cat("\nRunning initial QC...\n")

# 3a. Median methylated vs unmethylated intensity per sample (plotted in Script 03)
mSet <- preprocessRaw(rgSet)
qc   <- getQC(mSet)

# 3b. Sex prediction
cat("\nPredicting sex from methylation data...\n")
mSet_mapped   <- mapToGenome(mSet)
sex_predicted <- getSex(mSet_mapped, cutoff = -2)

# Compare predicted vs reported sex (R12: "Female"/"Male")
sex_check <- data.frame(
  GP2ID         = targets$GP2ID,
  clinical_id   = targets$clinical_id,
  Dataset       = targets$Dataset,
  reported_sex  = targets$sex,
  predicted_sex = sex_predicted$predictedSex
) %>%
  mutate(
    reported_sex_label = case_when(
      reported_sex == "Female" ~ "F",
      reported_sex == "Male"   ~ "M",
      TRUE                     ~ "Unknown"
    ),
    sex_discordant = reported_sex_label != predicted_sex &
                     reported_sex_label != "Unknown"
  )

cat("Sex discordant samples:", sum(sex_check$sex_discordant, na.rm = TRUE), "\n")
if (sum(sex_check$sex_discordant, na.rm = TRUE) > 0) {
  cat("Discordant samples:\n")
  print(sex_check %>% filter(sex_discordant))
}

# 3c. Beta density curves BEFORE normalization, computed as in
# minfi::densityPlot. Saved so the density figures can be redrawn without the
# full beta matrices. Curves are per sample, so compute them for all samples
# now and keep the QC-passed ones later; mSet can then be freed before the
# detection p-values
density_curves <- function(b) {
  d <- apply(b, 2, function(x) density(as.vector(x), na.rm = TRUE))
  list(x = sapply(d, `[[`, "x"), y = sapply(d, `[[`, "y"))
}
density_before <- density_curves(getBeta(mSet))
sample_names   <- colnames(mSet)   # idat basename, matches the beta matrix columns

rm(mSet, mSet_mapped)
invisible(gc())

# 3d. Detection p-values
# Computed in sample chunks to limit peak memory (detectionP is per-sample)
cat("\nComputing detection p-values...\n")
chunks <- split(seq_len(ncol(rgSet)), ceiling(seq_len(ncol(rgSet)) / DETP_CHUNK_SIZE))
detP <- do.call(cbind, lapply(chunks, function(i) {
  cat("  samples", min(i), "-", max(i), "\n")
  detectionP(rgSet[, i])
}))

# compute fraction of failed probes per sample
failed <- detP > DETECTION_P_THRESHOLD
cat("Fraction of failed probes per sample:\n")
print(round(colMeans(failed), 4))
cat("Probes failed in >50% of samples:",
    sum(rowMeans(failed) > 0.5), "\n")

# Mean detection p-value per sample for plotting and sample-level QC
mean_detP <- colMeans(detP)
cat("Mean detection p-value range:",
    round(min(mean_detP), 6), "to", round(max(mean_detP), 6), "\n")

# Identify failed samples, using the method chosen by USE_SESAME_QC in config.R
# (sesame_qc rows are in the same order as targets and the rgSet columns)
if (USE_SESAME_QC) {
  failed_samples <- sesame_qc$frac_dt_cg < SESAME_MIN_FRAC_DETECTED
  cat("Samples failing SeSAMe detection (fraction of cg probes detected <",
      SESAME_MIN_FRAC_DETECTED, "):", sum(failed_samples), "\n")
} else {
  failed_samples <- mean_detP > DETECTION_P_THRESHOLD
  cat("Samples failing minfi detection p-value threshold (mean p >",
      DETECTION_P_THRESHOLD, "):", sum(failed_samples), "\n")
}
if (sum(failed_samples) > 0) {
  cat("Failed samples:\n")
  print(targets[failed_samples, c("GP2ID", "clinical_id", "Dataset")])
}

# --- 4. Remove failed samples ------------------------------------------------

cat("\nRemoving failed samples...\n")
sex_discordant    <- sex_check$sex_discordant %in% TRUE
samples_to_remove <- failed_samples |
                     (REMOVE_SEX_DISCORDANT & sex_discordant)
cat("Removing sex discordant samples:", REMOVE_SEX_DISCORDANT, "\n")
cat("Total samples removed:", sum(samples_to_remove), "\n")
cat("Samples remaining:", sum(!samples_to_remove), "\n")

keep          <- !samples_to_remove
targets_clean <- targets[keep, ]

# Per-sample QC metrics for all input samples: everything plotted in the
# qc_01 to qc_03 figures, plus which samples failed and were removed
qc_metrics <- data.frame(
  Sample = sample_names,
  targets %>% select(GP2ID, GP2sampleID, clinical_id, phenotype, sex,
                     Dataset, Batch, Sentrix_ID, Sentrix_Position),
  mMed               = qc$mMed,                 # qc_01
  uMed               = qc$uMed,                 # qc_01
  mean_detP          = mean_detP,               # qc_02
  frac_dt_cg         = if (USE_SESAME_QC) sesame_qc$frac_dt_cg else NA,   # qc_02b
  xMed               = sex_predicted$xMed,      # qc_03
  yMed               = sex_predicted$yMed,      # qc_03
  predicted_sex      = sex_predicted$predictedSex,
  reported_sex_label = sex_check$reported_sex_label,
  sex_discordant     = sex_check$sex_discordant,
  failed_detection   = failed_samples,
  removed            = samples_to_remove,
  row.names = NULL
)
write.csv(qc_metrics, file.path(DIR_RESULTS, "qc_sample_metrics.csv"), row.names = FALSE)
cat("Per-sample QC metrics saved\n")

# Density curves of the QC-passed samples only
density_before <- list(x = density_before$x[, keep, drop = FALSE],
                       y = density_before$y[, keep, drop = FALSE])

# --- 5. Normalization --------------------------------------------------------

# Funnorm needs ~100 MB/sample of working memory on top of its input, so work
# out the probe filters that need the big objects (bead counts from the extended
# rgSet, detection p-values) now, then free everything before normalizing.
# Both filters are per probe, so computing them here gives the same result
# as after normalization.

# Probe annotation (probe name, bead addresses, chr, ...) for the array
ann <- getAnnotation(rgSet)

# Low bead count probes (used in 6a).
# getNBeads() rows are bead addresses, not probe names; map them to probes.
# Type I probes use two addresses (A and B), so take the lower count of the two
nbeads   <- getNBeads(rgSet)[, keep, drop = FALSE]
nbeads_A <- nbeads[match(as.character(ann$AddressA), rownames(nbeads)), , drop = FALSE]
nbeads_B <- nbeads[match(as.character(ann$AddressB), rownames(nbeads)), , drop = FALSE]
nbeads_probe <- ifelse(is.na(nbeads_B), nbeads_A, pmin(nbeads_A, nbeads_B))
rownames(nbeads_probe) <- ann$Name
# Probes with < MIN_BEADS beads in > 5% of samples
low_bead_probes <- rowSums(nbeads_probe < MIN_BEADS, na.rm = TRUE) > (0.05 * ncol(nbeads_probe))
low_bead_probes <- names(low_bead_probes)[low_bead_probes]

# Failed detection probes (used in 6c): detection p > DETECTION_P_THRESHOLD
# in more than FAILED_SAMPLE_CUTOFF of the kept samples
detP_clean <- detP[, keep, drop = FALSE]
failed_probe_names <- rownames(detP_clean)[
  rowSums(detP_clean > DETECTION_P_THRESHOLD) > (FAILED_SAMPLE_CUTOFF * ncol(detP_clean))]

# Slim RGChannelSet for Funnorm: only the Red/Green signal it uses (the
# extended rgSet also carries NBeads and SD matrices, ~2.5x the size)
rgSet_clean <- RGChannelSet(Green      = getGreen(rgSet)[, keep, drop = FALSE],
                            Red        = getRed(rgSet)[, keep, drop = FALSE],
                            colData    = colData(rgSet)[keep, ],
                            annotation = annotation(rgSet))

rm(rgSet, qc, detP, detP_clean, failed,
   nbeads, nbeads_A, nbeads_B, nbeads_probe)
invisible(gc())

# Apply functional normalization
cat("\nNormalizing with preprocessFunnorm...\n")
mSetSq <- preprocessFunnorm(rgSet_clean)
cat("Normalized object dimensions:", dim(mSetSq), "\n")
rm(rgSet_clean); invisible(gc())

# Save unfiltered betas for methylation clocks (clock CpGs may be removed by probe filters)
beta_norm <- getBeta(mSetSq)
saveRDS(beta_norm, BVALS_UNFILTERED)
cat("Unfiltered beta values saved\n")

# Density curves AFTER normalization
density_after <- density_curves(beta_norm)
rm(beta_norm); invisible(gc())

# Save everything needed to redraw the QC figures: per-sample metrics, density
# curves (columns named by Sample, as in qc_metrics) and the thresholds used
saveRDS(list(samples        = qc_metrics,
             density_before = density_before,
             density_after  = density_after,
             params         = list(DATA_SOURCE              = DATA_SOURCE,
                                   ARRAY                    = array_type,
                                   DETECTION_P_THRESHOLD    = DETECTION_P_THRESHOLD,
                                   USE_SESAME_QC            = USE_SESAME_QC,
                                   SESAME_MIN_FRAC_DETECTED = SESAME_MIN_FRAC_DETECTED)),
        QC_FIGURE_DATA)
cat("QC figure data saved to:", QC_FIGURE_DATA, "\n")

# --- 6. Probe Filtering -------------------------------------------------------

cat("\nFiltering probes...\n")
n_probes_start <- nrow(mSetSq)

# 6a. Remove probes with low bead count (list computed in section 5)
cat("Removing low bead count probes...\n")
mSetSq <- mSetSq[!rownames(mSetSq) %in% low_bead_probes, ]
n_after_beads <- nrow(mSetSq)
cat("Probes removed (low bead count):",
    n_probes_start - n_after_beads, "\n")

# 6b. Remove cross-reactive probes (list from the array profile)
cat("Removing cross-reactive probes...\n")
n_before_xreact <- nrow(mSetSq)
cross_reactive <- profile$xreactive(ann)
mSetSq <- mSetSq[!rownames(mSetSq) %in% cross_reactive, ]
n_after_xreact <- nrow(mSetSq)
cat("Probes removed (cross-reactive):",
    n_before_xreact - n_after_xreact, "\n")

# 6c. Remove failed probes (list computed in section 5)
cat("Removing failed probes...\n")
mSetSq <- mSetSq[!rownames(mSetSq) %in% failed_probe_names, ]
n_after_failed <- nrow(mSetSq)
cat("Probes removed (failed detection):", n_after_xreact - n_after_failed, "\n")

# 6d. Remove SNP-overlapping probes
cat("Removing SNP-overlapping probes...\n")
mSetSq <- dropLociWithSnps(mSetSq)
n_after_snp <- nrow(mSetSq)
cat("Probes removed (SNP-overlapping):", n_after_failed - n_after_snp, "\n")

# 6e. Remove sex chromosome probes
cat("Removing sex chromosome probes...\n")
sex_probes <- ann$Name[ann$chr %in% c("chrX", "chrY")]
mSetSq     <- mSetSq[!rownames(mSetSq) %in% sex_probes, ]
n_after_sex <- nrow(mSetSq)
cat("Probes removed (sex chromosomes):", n_after_snp - n_after_sex, "\n")

cat("Total probes removed:", n_probes_start - n_after_sex, "\n")
cat("Final probe count:", n_after_sex, "\n")

# --- 7. Calculate M/Beta values and Save QC-passed data -------------------------

cat("\nSaving QC-passed data...\n")

cat("\nCalculating M and Beta values...\n")
mVals <- getM(mSetSq)
bVals <- getBeta(mSetSq)

cat("M values dimensions:", dim(mVals), "\n")
cat("Beta values dimensions:", dim(bVals), "\n")

# Save as RDS
saveRDS(mVals, file.path(DIR_RESULTS, "mVals.rds"))
saveRDS(bVals, file.path(DIR_RESULTS, "bVals.rds"))
cat("M and Beta values saved\n")

# Save normalized filtered GenomicRatioSet
saveRDS(mSetSq, MSET_QC)

# Save cleaned sample sheet
write.csv(targets_clean, SAMPLE_SHEET_QC, row.names = FALSE)

# Save QC summary (progressive counts)
qc_summary <- data.frame(
  step = c("Input samples",
           if (USE_SESAME_QC) "Failed detection (SeSAMe pOOBAH)" else "Failed detection (minfi mean p-value)",
           if (REMOVE_SEX_DISCORDANT) "Sex discordant (removed)" else "Sex discordant (kept)",
           "Final samples",
           "Input probes",
           "Low bead count probes removed",
           "Cross-reactive probes removed",
           "Failed probes removed",
           "SNP probes removed",
           "Sex chromosome probes removed",
           "Final probes"),
  n    = c(nrow(targets),
           sum(failed_samples),
           sum(sex_discordant),
           nrow(targets_clean),
           n_probes_start,
           n_probes_start - n_after_beads,
           n_before_xreact - n_after_xreact,
           n_after_xreact - n_after_failed,
           n_after_failed - n_after_snp,
           n_after_snp - n_after_sex,
           n_after_sex)
)

write.csv(qc_summary,
          file.path(DIR_RESULTS, "qc_summary.csv"),
          row.names = FALSE)

cat("QC summary saved\n")
cat("\nQC pipeline complete!\n")
cat("Results saved to:", DIR_RESULTS, "\n")
