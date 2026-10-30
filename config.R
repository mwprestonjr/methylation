# =============================================================================
# Configuration - shared
# Author: GP2 Subtypes and Mechanisms - M.E., M.P.
# Date: April 23, 2026
# Updated: Oct 30, 2026
# Description: Settings shared by all modules: the data sources (input paths,
#              output directory), figure settings, and the files one module
#              writes and another reads. Settings used by only one module are
#              in that module's config (qc/config.R, preprocessing/config.R,
#              clocks/config.R, ewas/config.R, mqtl/config.R). Scripts source
#              this file, then their module's:
#                source("config.R"); source("qc/config.R")
#              and are run from the repository root, e.g.
#                Rscript qc/2_qc.R ppmi_p140
# =============================================================================

# --- Data source -------------------------------------------------------------

# Which data source this run processes. The qc/1 scripts set DATA_SOURCE
# themselves; the other scripts take it as the first command-line argument. In
# an interactive session, set DATA_SOURCE before sourcing this file
DATA_SOURCES <- c("ppmi_p140", "psomagen")
if (!exists("DATA_SOURCE")) DATA_SOURCE <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(DATA_SOURCE) || !DATA_SOURCE %in% DATA_SOURCES) {
  stop("Set the data source, one of: ", paste(DATA_SOURCES, collapse = ", "),
       "\n  e.g. Rscript qc/2_qc.R ", DATA_SOURCES[1])
}
cat("Data source:", DATA_SOURCE, "\n")

# --- Source-specific paths ---------------------------------------------------

if (DATA_SOURCE == "ppmi_p140") {
  # PPMI Project 140 (EPICv1). Sample sheet: qc/1_build_sample_sheet_ppmi.R
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
}

if (DATA_SOURCE == "psomagen") {
  # Psomagen deliveries (EPICv2). Sample sheet: qc/1_build_sample_sheet_psomagen.R
  # Datasets to process together; each has <DIR_DELIVERY>/<ID>/<ID>_QC_table.csv
  # and its idat files under <DIR_DELIVERY>/<ID>/<Sentrix_ID>/
  DATASETS     <- c("AB00000952",
                    "AB00000963")
  DIR_DELIVERY <- "/mnt/psomagen_delivery/nba-samples"
}

# Output folder per data source. Listed for every source, not only
# DATA_SOURCE, so a script can read another source's outputs (e.g.
# EPICV2_SIM_COMPARISON in clocks/config.R)
DIR_OUTPUTS <- c(ppmi_p140 = "/mnt/output/methylation/ppmi",
                 psomagen  = "/mnt/output/methylation/psomagen")
stopifnot(all(DATA_SOURCES %in% names(DIR_OUTPUTS)))
DIR_OUTPUT <- DIR_OUTPUTS[[DATA_SOURCE]]

# Output folders: one per module under results/ and figures/
#   <DIR_OUTPUT>/results/<module>   (mqtl: one subfolder per ancestry,
#   <DIR_OUTPUT>/figures/<module>    set in R/mqtl_setup.R)
DIR_RESULTS <- file.path(DIR_OUTPUT, "results")
DIR_FIGURES <- file.path(DIR_OUTPUT, "figures")
DIR_RESULTS_QC            <- file.path(DIR_RESULTS, "qc")
DIR_RESULTS_PREPROCESSING <- file.path(DIR_RESULTS, "preprocessing")
DIR_RESULTS_CLOCKS        <- file.path(DIR_RESULTS, "clocks")
DIR_RESULTS_MQTL          <- file.path(DIR_RESULTS, "mqtl")
DIR_RESULTS_EWAS          <- file.path(DIR_RESULTS, "ewas")
DIR_FIGURES_QC            <- file.path(DIR_FIGURES, "qc")
DIR_FIGURES_PREPROCESSING <- file.path(DIR_FIGURES, "preprocessing")
DIR_FIGURES_CLOCKS        <- file.path(DIR_FIGURES, "clocks")
DIR_FIGURES_MQTL          <- file.path(DIR_FIGURES, "mqtl")
DIR_FIGURES_EWAS          <- file.path(DIR_FIGURES, "ewas")
for (d in c(DIR_RESULTS_QC, DIR_RESULTS_PREPROCESSING, DIR_RESULTS_CLOCKS, DIR_RESULTS_MQTL, DIR_RESULTS_EWAS,
            DIR_FIGURES_QC, DIR_FIGURES_PREPROCESSING, DIR_FIGURES_CLOCKS, DIR_FIGURES_MQTL, DIR_FIGURES_EWAS)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# --- Figures -----------------------------------------------------------------

# Figure size for PNGs
FIG_WIDTH  <- 6    # inches
FIG_HEIGHT <- 4    # inches
FIG_RES    <- 300  # dpi

# --- Files and settings shared between modules -------------------------------

# GP2 R12 master key: GP2ID, sex, race, phenotype, ages (onset, diagnosis),
# family history, genotype availability, GP2 QC and ancestry labels. Used for
# the sample sheets (qc/1 scripts), clinical traits (clocks/) and the mQTL
# sample map (mqtl/). From the gp2_release12 bucket, mounted at
# ~/gp2_release12 (see README.md).
# NOTE: this filename may change from time to time
GP2_MASTER_KEY <- "~/gp2_release12/nba/CURRENT_master_key_nba_wgs_28_06_2026.csv"

# GP2 R12 internal-use master key: as GP2_MASTER_KEY, plus age at sample
# collection. Used only for Psomagen age (qc/1_build_sample_sheet_psomagen.R).
GP2_INTERNAL_KEY <- "~/gp2_release12/INTERNAL_USE_ONLY_master_key_release12_final_vwb.csv"

# GP2 R12 extended clinical data: one row per participant (GP2ID) and visit
# (visit_month). Used for the clinical traits (clocks/) and smoking status
# (ewas/). From the gp2_release12 bucket, mounted at ~/gp2_release12 (see
# README.md)
CLINICAL_DATA <- "~/gp2_release12/clinical/r12_extended_clinical_data.csv"

# Columns every qc/1 script writes to the sample sheet (source-specific
# extras may follow). Later scripts only rely on these
SAMPLE_SHEET_COLUMNS <- c("GP2ID", "GP2sampleID", "clinical_id",
                          "phenotype", "sex", "race", "age",
                          "Dataset", "Batch", "Sentrix_ID", "Sentrix_Position",
                          "Array", "Basename")

# qc/ outputs (also <DIR_RESULTS_QC>/bVals.rds)
SAMPLE_SHEET_QC  <- file.path(DIR_RESULTS_QC, "sample_sheet_qc_passed.csv")
MSET_QC          <- file.path(DIR_RESULTS_QC, "mSetSq_qc_passed.rds")
MVALS_QC         <- file.path(DIR_RESULTS_QC, "mVals.rds")             # normalized, probe-filtered, replicate probes averaged
REPLICATE_PROBES <- file.path(DIR_RESULTS_QC, "replicate_probes.csv")  # EPICv2 probes averaged per site in MVALS_QC (for mQTL)
BVALS_UNFILTERED <- file.path(DIR_RESULTS_QC, "bVals_unfiltered.rds")  # normalized, no probe filtering (for clocks)
BVALS_NOOB       <- file.path(DIR_RESULTS_QC, "bVals_noob.rds")        # noob-normalized, all probes (for cell deconvolution)

# preprocessing/ outputs (also CPG_ANNOTATION below); read by clocks/, ewas/
# and mqtl/
COMBAT_MVALS       <- file.path(DIR_RESULTS_PREPROCESSING, "combat_mVals.rds")
SAMPLE_SHEET_FINAL <- file.path(DIR_RESULTS_PREPROCESSING, "sample_sheet_final.csv")  # QC-passed sheet + cell proportions

# Blood cell types estimated with the IDOL reference (FlowSorted.Blood.EPIC);
# their proportions are covariates in later analyses
CELL_TYPES <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")

# --- CpG positions and gene annotation (R/cpg_annotation.R) ------------------

# hg19 -> hg38 liftover for arrays whose annotation is hg19 (EPICv1)
CHAIN_HG19_HG38 <- "/mnt/expansion_working/chain_files/hg19ToHg38.chain.gz"

# GENCODE genes (GRCh38) for annotating CpGs, built once per data source by
# preprocessing/4_cpg_annotation.R into CPG_ANNOTATION: the gene whose
# promoter (PROMOTER_UPSTREAM bp upstream to PROMOTER_DOWNSTREAM bp downstream
# of the transcription start) or body contains the CpG, else the nearest gene
GENCODE_GTF         <- "/mnt/expansion_working/methylation/gencode.v50.annotation.gtf.gz"
GENE_TYPES          <- c("protein_coding", "lncRNA")
PROMOTER_UPSTREAM   <- 1500
PROMOTER_DOWNSTREAM <- 500
CPG_ANNOTATION      <- file.path(DIR_RESULTS_PREPROCESSING, "cpg_annotation.csv")
