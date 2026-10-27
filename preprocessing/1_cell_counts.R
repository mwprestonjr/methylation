# =============================================================================
# Preprocessing 1: Cell type deconvolution
# Author: GP2 Subtypes and Mechanisms - M.E., M.P.
# Date: May 6, 2026
# Updated: Oct 16, 2026
# Description: Estimates blood cell type proportions (CELL_TYPES) for the
#              QC-passed samples with the IDOL reference (FlowSorted.Blood.EPIC,
#              Salas et al. 2018), from the noob-normalized betas saved by
#              qc/2_qc.R (BVALS_NOOB): the IDOL reference needs noob, not
#              Funnorm, values. Works for EPICv1 and EPICv2
# Usage:       Rscript preprocessing/1_cell_counts.R <data source>   (see config.R)
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

library(minfi)
library(tidyverse)
library(FlowSorted.Blood.EPIC)    # IDOL reference: IDOLOptimizedCpGs(.compTable)

# Load shared configuration and beta/M value helpers
source("config.R")
source("preprocessing/config.R")
source("R/mvalues.R")

# --- 1. Load QC-passed sample sheet ------------------------------------------

cat("Loading QC-passed sample sheet...\n")
targets <- read.csv(SAMPLE_SHEET_QC,
                    colClasses = c(GP2ID            = "character",
                                   clinical_id      = "character",
                                   Sentrix_ID       = "character",
                                   Sentrix_Position = "character",
                                   Basename         = "character"))
cat("Samples:", nrow(targets), "\n")

cat("Array:", unique(targets$Array), "\n")

# --- 2. Noob-normalized betas (qc/2_qc.R) ------------------------------------

# Saved before density outliers are removed, so subset to the QC-passed
# samples (columns are idat basenames)
cat("\nLoading noob-normalized betas...\n")
beta_noob <- readRDS(BVALS_NOOB)
samples   <- basename(targets$Basename)
stopifnot(all(samples %in% colnames(beta_noob)))
beta_noob <- beta_noob[, samples, drop = FALSE]
cat("Noob betas:", dim(beta_noob), "\n")

# --- 3. Estimate cell type proportions ---------------------------------------

# estimateCellCounts2 only supports 450k/EPICv1, so run its IDOL steps
# directly (the same for every array): noob betas, reduced to cg IDs,
# projected onto the IDOL reference
cat("\nEstimating cell type proportions with IDOL CpGs...\n")

# EPICv2 names carry a suffix (cg00000029_TC21); average replicates per cg ID
# (collapse_replicate_probes in R/mvalues.R) and name rows by the cg ID.
# EPICv1 names have no suffix, so this leaves them unchanged
cg_id     <- sub("_.*$", "", rownames(beta_noob))
beta_idol <- collapse_replicate_probes(beta_noob[cg_id %in% IDOLOptimizedCpGs, , drop = FALSE])$beta
rownames(beta_idol) <- sub("_.*$", "", rownames(beta_idol))
cat("IDOL CpGs found on array:", nrow(beta_idol), "/", length(IDOLOptimizedCpGs), "\n")
rm(beta_noob); invisible(gc())

cell_counts <- projectCellType_CP(
  Y           = beta_idol,
  coefWBC     = IDOLOptimizedCpGs.compTable[rownames(beta_idol), CELL_TYPES],
  lessThanOne = FALSE
)

# One row per sample, in sample sheet order; Sample = idat basename, which
# matches the columns of the M value matrices
cell_proportions <- data.frame(Sample = rownames(cell_counts), cell_counts,
                               row.names = NULL)
cell_proportions <- cell_proportions[match(basename(targets$Basename),
                                           cell_proportions$Sample), ]
stopifnot(all(cell_proportions$Sample == basename(targets$Basename)))

cat("Cell type proportions estimated for", nrow(cell_proportions), "samples\n")
cat("Mean proportion per cell type:\n")
print(round(colMeans(cell_proportions[, CELL_TYPES]), 4))

# --- 4. Plot and save --------------------------------------------------------

cell_long <- bind_cols(targets %>% select(phenotype), cell_proportions) %>%
  pivot_longer(cols = all_of(CELL_TYPES),
               names_to  = "cell_type",
               values_to = "proportion")

# plot boxplot of cell type proportions
png(file.path(DIR_FIGURES_PREPROCESSING, "cell_proportions.png"),
    width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
print(ggplot(cell_long, aes(x = cell_type, y = proportion)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(title = "Cell Type Proportions",
       x     = "Cell Type",
       y     = "Estimated Proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()

# plot boxplot of cell type proportions by phenotype
png(file.path(DIR_FIGURES_PREPROCESSING, "cell_proportions_by_phenotype.png"),
    width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
print(ggplot(cell_long, aes(x = cell_type, y = proportion, fill = phenotype)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  labs(title = "Cell Type Proportions by Phenotype",
       fill  = "Phenotype",
       x     = "Cell Type",
       y     = "Estimated Proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()

write.csv(cell_proportions, CELL_PROPORTIONS, row.names = FALSE)
cat("\nCell proportions saved to:", CELL_PROPORTIONS, "\n")
cat("Next: Rscript preprocessing/2_variation_sources.R", DATA_SOURCE, "\n")
