# =============================================================================
# Quality Control
# Author: GP2 Subtypes and Mechanisms - M.E., M.P.
# Date: Sept 22, 2026
# Updated: Sept 24, 2026
# Description: Reads idat files of all datasets (merged sample sheet from
#              Script 01) into minfi, runs QC following the minfi
#              user guide, filters failed samples and probes, normalizes
#              using preprocessFunnorm, and saves QC-passed data
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------
library(minfi)
library(tidyverse)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(maxprobes)

# Load shared configuration
source("~/methylation/scripts/00_config.R")

# --- 1. Load sample sheet ----------------------------------------------------

cat("Loading sample sheet...\n")
targets <- read.csv(SAMPLE_SHEET,
                    colClasses = c(Sentrix_ID = "character",
                                   Sentrix_Position = "character",
                                   Basename = "character")) 
# Quick checks
cat("Samples loaded:", nrow(targets), "\n")
cat("Samples per dataset and array:\n")
print(table(targets$Dataset, targets$Array))
cat("PD:", sum(targets$GP2_phenotype == "PD"), "\n")
cat("Control:", sum(targets$GP2_phenotype == "Control"), "\n")
cat("Other phenotype:", sum(targets$GP2_phenotype != "PD" & targets$GP2_phenotype != "Control"), "\n")
cat("Female:", sum(targets$sex == "Female"), "\n")
cat("Male:", sum(targets$sex == "Male"), "\n")
cat("\nSex by Diagnosis breakdown:\n")
print(table(targets$GP2_phenotype, targets$sex))

# The annotation and probe filters below are EPICv2-specific, and minfi can't
# read different array types into one RGChannelSet
if (!all(targets$Array == "EPICv2")) {
  stop("Sample sheet contains non-EPICv2 arrays (",
       paste(setdiff(unique(targets$Array), "EPICv2"), collapse = ", "),
       "); only EPICv2 is supported for now - see the array table from Script 01")
}

# TEMP reduce size for testing (need ot increase machine memory)
# take the first samples of each dataset so the merge is exercised
# targets <- targets %>%
#   group_by(Dataset) %>%
#   slice_head(n = N_SAMPLES_TESTING) %>%
#   ungroup() %>%
#   as.data.frame()

# --- 2. Read IDAT files ------------------------------------------------------

cat("\nReading IDAT files into minfi...\n")
rgSet <- read.metharray.exp(targets = targets,
                        verbose   = TRUE,
                        force     = TRUE, 
                        extended = TRUE)

# Verify metadata attached correctly 
cat("RGChannelSet dimensions:", dim(rgSet), "\n")
# cat("Array type:", annotation(rgSet)["array"], "\n")
cat("Sample metadata columns:", ncol(pData(rgSet)), "\n")
cat("First few GP2IDs:", head(pData(rgSet)$GP2ID), "\n")
cat("Phenotype breakdown:\n")
print(table(pData(rgSet)$GP2_phenotype))

# annotate the RGChannelSet with the appropriate array and annotation version
annotation(rgSet) <- c(array = "IlluminaHumanMethylationEPICv2",
                       annotation = "20a1.hg38")

# --- 3. Initial QC -----------------------------------------------------------

cat("\nRunning initial QC...\n")

# 3a. QC plot - median methylated vs unmethylated intensity per sample
cat("Generating QC plot (median intensities)...\n")
mSet <- preprocessRaw(rgSet)
qc   <- getQC(mSet)

png(file.path(DIR_RESULTS, "qc_01_median_intensities.png"), width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
plotQC(qc)
dev.off()
cat("QC plot saved\n")

# 3b. Detection p-values
cat("Computing detection p-values...\n")
detP <- detectionP(rgSet)

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

# Plot mean detection p-value per sample, sorted, on a log scale so the
# threshold (usually far above typical values) stays visible.
# Only failing samples are labelled (by GP2ID)
ord      <- order(mean_detP, decreasing = TRUE)
detP_ord <- mean_detP[ord]
failing  <- detP_ord > DETECTION_P_THRESHOLD
y_lim    <- c(10^floor(log10(min(detP_ord))),
              10 * max(detP_ord, DETECTION_P_THRESHOLD))
y_ticks  <- 10^seq(log10(y_lim[1]), log10(y_lim[2]))

png(file.path(DIR_RESULTS, "qc_02_detection_pvalues.png"), width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
par(mar = c(3, 5.5, 3, 1), las = 1)
plot(detP_ord,
     log  = "y",
     ylim = y_lim,
     pch  = 16,
     cex  = 0.8,
     col  = ifelse(failing, "red", "grey30"),
     xaxt = "n",
     yaxt = "n",
     xlab = "",
     ylab = "",
     main = "Mean Detection P-value per Sample")
axis(2, at = y_ticks, labels = format(y_ticks, scientific = TRUE))
title(ylab = "Mean detection p-value (log scale)", line = 4)
mtext(sprintf("Samples, sorted (n = %d)", length(detP_ord)), side = 1, line = 1)
abline(h   = DETECTION_P_THRESHOLD,
       col = "red",
       lty = 2)
if (any(failing)) {
  text(which(failing), detP_ord[failing],
       labels = targets$GP2ID[ord][failing],
       pos    = 4,
       cex    = 0.6,
       col    = "red")
}
legend("topright",
       legend = sprintf("Threshold (%g)", DETECTION_P_THRESHOLD),
       lty    = 2,
       col    = "red",
       bty    = "n")
dev.off()
cat("Detection p-value plot saved\n")

# Identify failed samples
failed_samples <- mean_detP > DETECTION_P_THRESHOLD
cat("Samples failing detection p-value threshold:", sum(failed_samples), "\n")
if (sum(failed_samples) > 0) {
  cat("Failed samples:\n")
  print(targets[failed_samples, c("GP2ID", "Dataset")])
}

# 3c. Sex prediction
cat("\nPredicting sex from methylation data...\n")
mSet_mapped   <- mapToGenome(mSet)
sex_predicted <- getSex(mSet_mapped, cutoff = -2)

# Compare predicted vs reported sex (R12: "Female"/"Male")
sex_check <- data.frame(
  GP2ID         = targets$GP2ID,
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

# Plot sex prediction
png(file.path(DIR_RESULTS, "qc_03_sex_prediction.png"), width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
plot(sex_predicted$xMed, sex_predicted$yMed,
     col  = ifelse(sex_predicted$predictedSex == "F", "hotpink", "steelblue"),
     pch  = 16,
     xlab = "X chromosome median intensity",
     ylab = "Y chromosome median intensity",
     main = "Sex Prediction")
# Add text labels for discordant samples
discordant_idx <- !is.na(sex_check$sex_discordant) & sex_check$sex_discordant
if (any(discordant_idx)) {
  text(sex_predicted$xMed[discordant_idx], 
       sex_predicted$yMed[discordant_idx],
       labels = sex_check$GP2ID[discordant_idx],
       pos    = 3,
       cex    = 0.7,
       col    = "red")
}
legend("topright", 
       legend = c("Female", "Male", "Discordant"),
       col    = c("hotpink", "steelblue", "red"), 
       pch    = c(16, 16, 16))
dev.off()
cat("Sex prediction plot saved\n")

# --- 4. Remove failed samples ------------------------------------------------

cat("\nRemoving failed samples...\n")
samples_to_remove <- failed_samples
cat("Total samples removed:", sum(samples_to_remove), "\n")
cat("Samples remaining:", sum(!samples_to_remove), "\n")

rgSet_clean   <- rgSet[, !samples_to_remove]
detP_clean    <- detP[,  !samples_to_remove]
targets_clean <- targets[!samples_to_remove, ]

# --- 5. Normalization --------------------------------------------------------

cat("\nNormalizing with preprocessFunnorm...\n")

# Density plot BEFORE normalization
png(file.path(DIR_RESULTS, "qc_04_density_before_normalization.png"), width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
densityPlot(getBeta(preprocessRaw(rgSet_clean)),
            sampGroups = targets_clean$GP2_phenotype,
            main       = "Beta Values - Before Normalization",
            legend     = TRUE)
dev.off()

# Same, grouped by dataset, to spot dataset-level shifts before normalization
png(file.path(DIR_RESULTS, "qc_04b_density_by_dataset.png"), width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
densityPlot(getBeta(preprocessRaw(rgSet_clean)),
            sampGroups = targets_clean$Dataset,
            main       = "Beta Values by Dataset - Before Normalization",
            legend     = TRUE)
dev.off()

# Apply functional normalization
mSetSq <- preprocessFunnorm(rgSet_clean)
cat("Normalized object dimensions:", dim(mSetSq), "\n")

# Save unfiltered betas for methylation clocks (clock CpGs may be removed by probe filters)
saveRDS(getBeta(mSetSq), BVALS_UNFILTERED)
cat("Unfiltered beta values saved\n")

# Density plot AFTER normalization
png(file.path(DIR_RESULTS, "qc_05_density_after_normalization.png"), width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
densityPlot(getBeta(mSetSq),
            sampGroups = targets_clean$GP2_phenotype,
            main       = "Beta Values - After Normalization",
            legend     = TRUE)
dev.off()
cat("Density plots saved\n")

# --- 6. Probe Filtering -------------------------------------------------------

cat("\nFiltering probes...\n")
n_probes_start <- nrow(mSetSq)

# EPICv2 probe annotation (probe name, bead addresses, chr, v1 name)
ann_v2 <- getAnnotation(rgSet_clean)

# 6a. Remove probes with low bead count
cat("Removing low bead count probes...\n")
# getNBeads() rows are bead addresses, not probe names; map them to probes.
# Type I probes use two addresses (A and B), so take the lower count of the two
nbeads   <- getNBeads(rgSet_clean)
nbeads_A <- nbeads[as.character(ann_v2$AddressA), , drop = FALSE]
nbeads_B <- nbeads[match(as.character(ann_v2$AddressB), rownames(nbeads)), , drop = FALSE]
nbeads_probe <- ifelse(is.na(nbeads_B), nbeads_A, pmin(nbeads_A, nbeads_B))
rownames(nbeads_probe) <- ann_v2$Name
# Remove probes with < MIN_BEADS beads in > 5% of samples
low_bead_probes <- rowSums(nbeads_probe < MIN_BEADS) > (0.05 * ncol(nbeads_probe))
# Only filter probes that exist in mSetSq
low_bead_probes <- names(low_bead_probes)[low_bead_probes]
mSetSq <- mSetSq[!rownames(mSetSq) %in% low_bead_probes, ]
n_after_beads <- nrow(mSetSq)
cat("Probes removed (low bead count):", 
    n_probes_start - n_after_beads, "\n")

# 6b. Remove cross-reactive probes
cat("Removing cross-reactive probes...\n")
n_before_xreact <- nrow(mSetSq)
# maxprobes only has an EPICv1 list; map it to EPICv2 names via EPICv1_Loci.
# NOTE: probes new on EPICv2 (no v1 equivalent) are not screened here
cross_reactive_v1 <- unlist(maxprobes::xreactive_probes(array_type = "EPIC"))
cross_reactive <- ann_v2$Name[ann_v2$EPICv1_Loci %in% cross_reactive_v1]
mSetSq <- mSetSq[!rownames(mSetSq) %in% cross_reactive, ]
n_after_xreact <- nrow(mSetSq)
cat("Probes removed (cross-reactive):", 
    n_before_xreact - n_after_xreact, "\n")

# 6c. Remove failed probes
cat("Removing failed probes...\n")
# Ensure detP matches current mSetSq probes after bead/xreact filtering
detP_clean <- detP_clean[rownames(detP_clean) %in% rownames(mSetSq), ]
detP_clean <- detP_clean[match(rownames(mSetSq), rownames(detP_clean)), ]
failed_probes <- rowSums(detP_clean > DETECTION_P_THRESHOLD) >
                 (FAILED_SAMPLE_CUTOFF * ncol(detP_clean))
mSetSq <- mSetSq[!failed_probes, ]
n_after_failed <- nrow(mSetSq)
cat("Probes removed (failed detection):", sum(failed_probes), "\n")

# 6d. Remove SNP-overlapping probes
cat("Removing SNP-overlapping probes...\n")
mSetSq <- dropLociWithSnps(mSetSq)
n_after_snp <- nrow(mSetSq)  
cat("Probes removed (SNP-overlapping):", n_after_failed - n_after_snp, "\n")

# 6e. Remove sex chromosome probes
cat("Removing sex chromosome probes...\n")
sex_probes <- ann_v2$Name[ann_v2$chr %in% c("chrX", "chrY")]
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

# Save QC summary (fix 3: use progressive counts)
qc_summary <- data.frame(
  step = c("Input samples",
           "Failed detection p-value",
           "Sex discordant",
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
           sum(sex_check$sex_discordant, na.rm = TRUE),
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