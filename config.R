# =============================================================================
# Configuration
# Author: GP2 Subtypes and Mechanisms - M.E., M.P.
# Date: April 23, 2026
# Updated: Sept 29, 2026
# Description: Defines shared paths and parameters used across all scripts.
#              Settings that differ by data source (input paths, output
#              directory) are chosen by DATA_SOURCE; everything else is shared.
#              Scripts are run from the repository root, e.g.
#                Rscript scripts/02_qc.R ppmi_p140
# =============================================================================

# --- Data source -------------------------------------------------------------

# Which data source this run processes. The 01 scripts set DATA_SOURCE
# themselves; Scripts 02+ take it as the first command-line argument. In an
# interactive session, set DATA_SOURCE before sourcing this file
DATA_SOURCES <- c("ppmi_p140", "psomagen")
if (!exists("DATA_SOURCE")) DATA_SOURCE <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(DATA_SOURCE) || !DATA_SOURCE %in% DATA_SOURCES) {
  stop("Set the data source, one of: ", paste(DATA_SOURCES, collapse = ", "),
       "\n  e.g. Rscript scripts/02_qc.R ", DATA_SOURCES[1])
}
cat("Data source:", DATA_SOURCE, "\n")

# --- Shared inputs -----------------------------------------------------------

# GP2 R12 master key: GP2ID, sex, race, phenotype, age at sample collection
FNAME_METADATA <- "/mnt/output/metadata/INTERNAL_USE_ONLY_master_key_release12_final_vwb.csv"

# --- Source-specific paths ---------------------------------------------------

if (DATA_SOURCE == "ppmi_p140") {
  # PPMI Project 140 (EPICv1). Sample sheet: 01_build_sample_sheet_ppmi.R
  PPMI_DIR           <- "/mnt/expansion_working/methylation/ppmi"
  P140_DIR           <- file.path(PPMI_DIR, "project_140")
  IDAT_DIR           <- file.path(P140_DIR, "idat")
  SUBJ_DIR           <- file.path(P140_DIR, "metadata", "subject_characteristics")
  LINK_LIST          <- file.path(PPMI_DIR, "documentation", "ppmi_140_link_list_20210607.csv")
  PARTICIPANT_STATUS <- file.path(SUBJ_DIR, "Participant_Status_23Apr2026.csv")
  AGE_AT_VISIT       <- file.path(SUBJ_DIR, "Age_at_visit_23Apr2026.csv")

  # Samples to include: PPMI enrollment cohort at the chosen visit
  #   COHORT 1 = Parkinson's Disease, 2 = Healthy Control,
  #   3 = SWEDD, 4 = Prodromal
  # (R12 GP2_phenotype can't be used here: it labels genetic-cohort unaffected
  # carriers "Control" and later converters "PD")
  PRIMARY_TIMEPOINT   <- "BL"     # Baseline only
  COHORTS_OF_INTEREST <- c(1, 2)  # PD vs Healthy Control

  # ComBat batch variable (see Preprocessing below): plate. Project 140
  # samples are spread over ~240 chips (median 2 samples per chip, many
  # alone), too few per chip for ComBat; the 22 plates hold 4-40 each
  COMBAT_BATCH_VAR <- "Batch"

  DIR_OUTPUT <- "/mnt/output/methylation/ppmi"
}

if (DATA_SOURCE == "psomagen") {
  # Psomagen deliveries (EPICv2). Sample sheet: 01_build_sample_sheet_psomagen.R
  # Datasets to process together; each has <DIR_DELIVERY>/<ID>/<ID>_QC_table.csv
  # and its idat files under <DIR_DELIVERY>/<ID>/<Sentrix_ID>/
  DATASETS     <- c("AB00000952",
                    "AB00000963")
  DIR_DELIVERY <- "/mnt/psomagen_delivery/nba-samples"

  # ComBat batch variable (see Preprocessing below): chip. ~8 samples per
  # chip, and chip explains more technical variation than delivery
  COMBAT_BATCH_VAR <- "Sentrix_ID"

  DIR_OUTPUT <- "/mnt/output/methylation/psomagen"
}

# Define/create derivative paths
DIR_RESULTS <- file.path(DIR_OUTPUT, "results")
DIR_FIGURES <- file.path(DIR_OUTPUT, "figures")
dir.create(DIR_RESULTS, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FIGURES, showWarnings = FALSE, recursive = TRUE)

# --- Analysis parameters -----------------------------------------------------

# QC thresholds
DETECTION_P_THRESHOLD <- 0.01   # Max detection p-value
MIN_BEADS             <- 3      # Minimum beads per probe
FAILED_SAMPLE_CUTOFF  <- 0.1    # Max fraction of failed probes per sample (ChAMP default)
REMOVE_SEX_DISCORDANT <- TRUE   # Remove samples whose predicted sex doesn't match reported sex

# Beta density outliers (Script 02), scored on the curves after normalization:
#   density_mid_peak = prominence of the largest peak between beta 0.15 and
#                      0.75 (a third bump, e.g. from a mixed sample); 0 = none
DENSITY_MID_PEAK_MAX    <- 0.05   # flag if density_mid_peak is above this
REMOVE_DENSITY_OUTLIERS <- FALSE  # TRUE = remove flagged samples; FALSE = flag only

# Sample-level detection QC method (Script 02)
#   TRUE  = SeSAMe pOOBAH: remove samples with < SESAME_MIN_FRAC_DETECTED of cg
#           probes detected; also saves qc_sesame_stats.csv.
#   FALSE = minfi detectionP: remove samples with mean p > DETECTION_P_THRESHOLD
# Probe-level filtering uses minfi detectionP either way
USE_SESAME_QC            <- TRUE
SESAME_MIN_FRAC_DETECTED <- 0.90   # Min fraction of cg probes detected (pOOBAH p < 0.05)

# Samples per detectionP() call in Script 02; computing all samples at once
# needs more memory than the VM has for a few hundred EPIC samples
DETP_CHUNK_SIZE <- 50

# Figure size for PNGs
FIG_WIDTH  <- 6    # inches
FIG_HEIGHT <- 4  # inches
FIG_RES    <- 300  # dpi

# --- Pipeline outputs --------------------------------------------------------
# These files are created by one script and read by the next
SAMPLE_SHEET        <- file.path(DIR_RESULTS, "sample_sheet.csv")
SAMPLE_SHEET_QC     <- file.path(DIR_RESULTS, "sample_sheet_qc_passed.csv")
MSET_QC             <- file.path(DIR_RESULTS, "mSetSq_qc_passed.rds")
BVALS_UNFILTERED    <- file.path(DIR_RESULTS, "bVals_unfiltered.rds")  # normalized, no probe filtering (for clocks)
QC_FIGURE_DATA      <- file.path(DIR_RESULTS, "qc_figure_data.rds")    # everything needed to redraw the Script 02 QC figures

# Columns every 01 script writes to SAMPLE_SHEET (source-specific extras may
# follow). Scripts 02+ only rely on these
SAMPLE_SHEET_COLUMNS <- c("GP2ID", "GP2sampleID", "clinical_id",
                          "phenotype", "sex", "race", "age",
                          "Dataset", "Batch", "Sentrix_ID", "Sentrix_Position",
                          "Array", "Basename")

# --- Preprocessing (preprocessing/ scripts) ----------------------------------
# Run per data source after qc/: 1_cell_counts.R, 2_variation_sources.R (raw),
# 3_combat.R, then 2_variation_sources.R again on the ComBat output

# Blood cell types estimated with the IDOL reference (FlowSorted.Blood.EPIC)
CELL_TYPES       <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")
CELL_PROPORTIONS <- file.path(DIR_RESULTS, "cell_proportions.csv")

# Sources of variation: SVD of the VARIATION_TOP_CPGS most variable CpGs,
# testing the first VARIATION_N_PCS components against each sample variable
VARIATION_TOP_CPGS <- 50000
VARIATION_N_PCS    <- 10

# ComBat: the batch variable is set per data source above (COMBAT_BATCH_VAR,
# a sample sheet column: "Sentrix_ID" = chip, "Batch" = plate for PPMI /
# delivery for Psomagen); choose it from 2_variation_sources.R. Every batch
# needs at least 2 samples. COMBAT_PROTECT: variables whose variation ComBat
# keeps (e.g. add "age", "sex", or the cell types)
COMBAT_PROTECT   <- c("phenotype")

COMBAT_MVALS       <- file.path(DIR_RESULTS, "combat_mVals.rds")
SAMPLE_SHEET_FINAL <- file.path(DIR_RESULTS, "sample_sheet_final.csv")  # QC-passed sheet + cell proportions
