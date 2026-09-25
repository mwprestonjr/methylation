# =============================================================================
# Configuration
# Author: GP2 Subtypes and Mechanisms - M.E., MP
# Date: April 23, 2026
# Updated: Sept 22, 2026
# Description: Defines shared paths and parameters used across all scripts
# =============================================================================

# --- Paths -------------------------------------------------------------------

# Set directories and paths
# Datasets to process together; each has <DIR_DELIVERY>/<ID>/<ID>_QC_table.csv
# and its idat files under <DIR_DELIVERY>/<ID>/<Sentrix_ID>/
DATASETS        <- c("AB00000952",
                     "AB00000963")
DIR_DELIVERY    <- "/mnt/psomagen_delivery/nba-samples"
DIR_OUTPUT      <- "/mnt/output/methylation" # path in which to save outputs and results
# FNAME_METADATA  <- "/mnt/output/metadata/R12_CURRENT_master_key_nba_wgs_20_06_2026.txt" # path to clinical and other metadata (R12)
FNAME_METADATA  <- "/mnt/output/metadata/INTERNAL_USE_ONLY_master_key_release12_final_vwb.csv" # path to clinical and other metadata (R12)

# Define/create derivative paths
DIR_RESULTS <- file.path(DIR_OUTPUT, "results")
dir.create(DIR_RESULTS, showWarnings = FALSE, recursive = TRUE)

# --- Analysis parameters -----------------------------------------------------

# QC thresholds
DETECTION_P_THRESHOLD <- 0.01   # Max detection p-value
MIN_BEADS             <- 3      # Minimum beads per probe
FAILED_SAMPLE_CUTOFF  <- 0.1    # Max fraction of failed probes per sample (ChAMP default)

# Sample-level detection QC method (Script 02)
#   TRUE  = SeSAMe pOOBAH: remove samples with < SESAME_MIN_FRAC_DETECTED of cg
#           probes detected; also saves qc_sesame_stats.csv.
#   FALSE = minfi detectionP: remove samples with mean p > DETECTION_P_THRESHOLD
# Probe-level filtering uses minfi detectionP either way
USE_SESAME_QC            <- TRUE
SESAME_MIN_FRAC_DETECTED <- 0.95   # Min fraction of cg probes detected (pOOBAH p < 0.05)

# Figure size for PNGs (16:9, fits a Google Slides slide)
FIG_WIDTH  <- 6    # inches
FIG_HEIGHT <- 4  # inches
FIG_RES    <- 300  # dpi

# --- Pipeline outputs --------------------------------------------------------
# These files are created by one script and read by the next
SAMPLE_SHEET        <- file.path(DIR_RESULTS, "sample_sheet.csv")
SAMPLE_SHEET_QC     <- file.path(DIR_RESULTS, "sample_sheet_qc_passed.csv")
MSET_QC             <- file.path(DIR_RESULTS, "mSetSq_qc_passed.rds")
BVALS_UNFILTERED    <- file.path(DIR_RESULTS, "bVals_unfiltered.rds")  # normalized, no probe filtering (for clocks)
