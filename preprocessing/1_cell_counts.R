# =============================================================================
# Preprocessing 1: Cell type deconvolution
# Author: GP2 Subtypes and Mechanisms - M.E., M.P.
# Date: May 6, 2026
# Updated: Sept 30, 2026
# Description: Estimates blood cell type proportions (CELL_TYPES) for the
#              QC-passed samples with the IDOL reference (FlowSorted.Blood.EPIC,
#              Salas et al. 2018). Re-reads the idat files, since this needs
#              noob-normalized raw intensities rather than the Funnorm output.
#              Works for any array in R/array_profiles.R
# Usage:       Rscript preprocessing/1_cell_counts.R <data source>   (see config.R)
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

library(minfi)
library(tidyverse)
library(FlowSorted.Blood.EPIC)    # IDOL reference: IDOLOptimizedCpGs(.compTable)

# Load shared configuration and array profiles
source("config.R")
source("R/array_profiles.R")

# --- 1. Load QC-passed sample sheet ------------------------------------------

cat("Loading QC-passed sample sheet...\n")
targets <- read.csv(SAMPLE_SHEET_QC,
                    colClasses = c(GP2ID            = "character",
                                   clinical_id      = "character",
                                   Sentrix_ID       = "character",
                                   Sentrix_Position = "character",
                                   Basename         = "character"))
cat("Samples:", nrow(targets), "\n")

array_type <- unique(targets$Array)
stopifnot(length(array_type) == 1, array_type %in% names(ARRAY_PROFILES))
profile <- ARRAY_PROFILES[[array_type]]
library(profile$anno_pkg, character.only = TRUE)
cat("Array:", array_type, "\n")

# --- 2. Re-read idat files ---------------------------------------------------

# Bead counts/SDs aren't needed here; extended = FALSE halves memory
cat("\nRe-reading idat files for cell type estimation...\n")
rgSet <- read.metharray.exp(targets  = targets,
                            verbose  = TRUE,
                            force    = TRUE,
                            extended = FALSE)
annotation(rgSet) <- profile$annotation

# --- 3. Estimate cell type proportions ---------------------------------------

# estimateCellCounts2 only supports 450k/EPICv1, so run its IDOL steps
# directly (the same for every array): noob-normalize, reduce probes to cg IDs,
# project onto the IDOL reference
cat("\nEstimating cell type proportions with IDOL CpGs...\n")
beta_noob <- getBeta(preprocessNoob(rgSet))
rm(rgSet); invisible(gc())

# EPICv2 names carry a suffix (cg00000029_TC21); average replicates per cg ID.
# EPICv1 names have no suffix, so this leaves them unchanged
cg_id     <- sub("_.*$", "", rownames(beta_noob))
keep      <- cg_id %in% IDOLOptimizedCpGs
beta_idol <- rowsum(beta_noob[keep, ], cg_id[keep]) /
             as.vector(rowsum(rep(1, sum(keep)), cg_id[keep]))
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

png(file.path(DIR_FIGURES, "cell_proportions.png"),
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
