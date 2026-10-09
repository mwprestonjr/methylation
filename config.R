# =============================================================================
# Configuration - shared
# Author: GP2 Subtypes and Mechanisms - M.E., M.P.
# Date: April 23, 2026
# Updated: Oct 7, 2026
# Description: Settings shared by all modules: the data sources (input paths,
#              output directory), figure settings, and the files one module
#              writes and another reads. Settings used by only one module are
#              in that module's config (qc/config.R, preprocessing/config.R,
#              mqtl/config.R). Scripts source this file, then their module's:
#                source("config.R"); source("qc/config.R")
#              and are run from the repository root, e.g.
#                Rscript qc/02_qc.R ppmi_p140
# =============================================================================

# --- Data source -------------------------------------------------------------

# Which data source this run processes. The 01 scripts set DATA_SOURCE
# themselves; the other scripts take it as the first command-line argument. In
# an interactive session, set DATA_SOURCE before sourcing this file
DATA_SOURCES <- c("ppmi_p140", "psomagen")
if (!exists("DATA_SOURCE")) DATA_SOURCE <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(DATA_SOURCE) || !DATA_SOURCE %in% DATA_SOURCES) {
  stop("Set the data source, one of: ", paste(DATA_SOURCES, collapse = ", "),
       "\n  e.g. Rscript qc/02_qc.R ", DATA_SOURCES[1])
}
cat("Data source:", DATA_SOURCE, "\n")

# --- Source-specific paths ---------------------------------------------------

if (DATA_SOURCE == "ppmi_p140") {
  # PPMI Project 140 (EPICv1). Sample sheet: qc/01_build_sample_sheet_ppmi.R
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

  DIR_OUTPUT <- "/mnt/output/methylation/ppmi"
}

if (DATA_SOURCE == "psomagen") {
  # Psomagen deliveries (EPICv2). Sample sheet: qc/01_build_sample_sheet_psomagen.R
  # Datasets to process together; each has <DIR_DELIVERY>/<ID>/<ID>_QC_table.csv
  # and its idat files under <DIR_DELIVERY>/<ID>/<Sentrix_ID>/
  DATASETS     <- c("AB00000952",
                    "AB00000963")
  DIR_DELIVERY <- "/mnt/psomagen_delivery/nba-samples"

  DIR_OUTPUT <- "/mnt/output/methylation/psomagen"
}

# Define/create derivative paths
DIR_RESULTS <- file.path(DIR_OUTPUT, "results")
DIR_FIGURES <- file.path(DIR_OUTPUT, "figures")
dir.create(DIR_RESULTS, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FIGURES, showWarnings = FALSE, recursive = TRUE)

# --- Figures -----------------------------------------------------------------

# Figure size for PNGs
FIG_WIDTH  <- 6    # inches
FIG_HEIGHT <- 4    # inches
FIG_RES    <- 300  # dpi

# --- Files and settings shared between modules -------------------------------

# GP2 R12 master key: GP2ID, sex, race, phenotype, ages (sample collection,
# onset, diagnosis), family history. Used for the sample sheets (qc/01
# scripts) and clinical traits (clocks/)
FNAME_METADATA <- "/mnt/output/metadata/INTERNAL_USE_ONLY_master_key_release12_final_vwb.csv"

# Columns every qc/01 script writes to the sample sheet (source-specific
# extras may follow). Later scripts only rely on these
SAMPLE_SHEET_COLUMNS <- c("GP2ID", "GP2sampleID", "clinical_id",
                          "phenotype", "sex", "race", "age",
                          "Dataset", "Batch", "Sentrix_ID", "Sentrix_Position",
                          "Array", "Basename")

# qc/ outputs (also <DIR_RESULTS>/mVals.rds and bVals.rds)
SAMPLE_SHEET_QC  <- file.path(DIR_RESULTS, "sample_sheet_qc_passed.csv")
MSET_QC          <- file.path(DIR_RESULTS, "mSetSq_qc_passed.rds")
BVALS_UNFILTERED <- file.path(DIR_RESULTS, "bVals_unfiltered.rds")  # normalized, no probe filtering (for clocks)

# preprocessing/ outputs
COMBAT_MVALS       <- file.path(DIR_RESULTS, "combat_mVals.rds")
SAMPLE_SHEET_FINAL <- file.path(DIR_RESULTS, "sample_sheet_final.csv")  # QC-passed sheet + cell proportions

# Blood cell types estimated with the IDOL reference (FlowSorted.Blood.EPIC);
# their proportions are covariates in later analyses
CELL_TYPES <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")
