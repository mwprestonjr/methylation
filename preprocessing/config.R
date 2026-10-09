# =============================================================================
# Configuration - preprocessing module (preprocessing/ scripts)
# Author: GP2 Subtypes and Mechanisms - M.E., M.P.
# Date: Oct 7, 2026
# Description: Settings used only by the preprocessing/ scripts, run per data
#              source after qc/: 1_cell_counts.R, 2_variation_sources.R (raw),
#              3_combat.R, then 2_variation_sources.R again on the ComBat
#              output. Source after config.R
# =============================================================================

# Cell type proportions (1_cell_counts.R -> 2, 3; cell types: CELL_TYPES in
# config.R)
CELL_PROPORTIONS <- file.path(DIR_RESULTS_PREPROCESSING, "cell_proportions.csv")

# Sources of variation: SVD of the VARIATION_TOP_CPGS most variable CpGs,
# testing the first VARIATION_N_PCS components against each sample variable
VARIATION_TOP_CPGS <- 50000
VARIATION_N_PCS    <- 10

# ComBat batch variable, per data source (a sample sheet column: "Sentrix_ID" =
# chip, "Batch" = plate for PPMI / delivery for Psomagen); choose it from
# 2_variation_sources.R. Every batch needs at least 2 samples
#   ppmi_p140: plate. Project 140 samples are spread over ~240 chips (median 2
#              samples per chip, many alone), too few per chip for ComBat; the
#              22 plates hold 4-40 each
#   psomagen:  chip. ~8 samples per chip, and chip explains more technical
#              variation than delivery
COMBAT_BATCH_VAR <- switch(DATA_SOURCE,
                           ppmi_p140 = "Batch",
                           psomagen  = "Sentrix_ID")

# Variables whose variation ComBat keeps (e.g. add "age", "sex", or the cell
# types)
COMBAT_PROTECT <- c("phenotype")
