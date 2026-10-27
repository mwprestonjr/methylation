# =============================================================================
# Configuration - QC module (qc/ scripts)
# Author: GP2 Subtypes and Mechanisms - M.E., M.P.
# Date: Oct 7, 2026
# Description: Settings used only by the qc/ scripts. Source after config.R
# =============================================================================

# --- Sample and probe QC thresholds (qc/2_qc.R) ------------------------------

DETECTION_P_THRESHOLD <- 0.01   # Max detection p-value
MIN_BEADS             <- 3      # Minimum beads per probe
FAILED_SAMPLE_CUTOFF  <- 0.1    # Max fraction of failed probes per sample (ChAMP default)
REMOVE_SEX_DISCORDANT <- TRUE   # Remove samples whose predicted sex doesn't match reported sex

# Beta density outliers, scored on the curves after normalization:
#   density_mid_peak = prominence of the largest peak between beta 0.15 and
#                      0.75 (a third bump, e.g. from a mixed sample); 0 = none
DENSITY_MID_PEAK_MAX    <- 0.05   # flag if density_mid_peak is above this
REMOVE_DENSITY_OUTLIERS <- FALSE  # TRUE = remove flagged samples; FALSE = flag only

# Sample-level detection QC method
#   TRUE  = SeSAMe pOOBAH: remove samples with < SESAME_MIN_FRAC_DETECTED of cg
#           probes detected; also saves qc_sesame_stats.csv.
#   FALSE = minfi detectionP: remove samples with mean p > DETECTION_P_THRESHOLD
# Probe-level filtering uses minfi detectionP either way
USE_SESAME_QC            <- TRUE
SESAME_MIN_FRAC_DETECTED <- 0.90   # Min fraction of cg probes detected (pOOBAH p < 0.05)

# Duplicate arrays: several QC-passed arrays of one person (same GP2ID, or
# clinical_id if there is no GP2ID). One is kept, the best by
# DUPLICATE_METRIC; ties go to the more recent chip (higher Sentrix ID)
#   "frac_dt_cg" = highest fraction of cg probes detected (SeSAMe pOOBAH;
#                  needs USE_SESAME_QC = TRUE)
#   "mean_detP"  = lowest mean minfi detection p-value
# Sample swap check: if any two arrays of the person have rs SNP probe betas
# correlating below DUPLICATE_MIN_SNP_COR (true replicates are ~0.95 or
# higher, different people much lower), or a different predicted sex, they
# can't all be that person, so all of them are removed
DUPLICATE_METRIC      <- "frac_dt_cg"
DUPLICATE_MIN_SNP_COR <- 0.9

# Samples per detectionP() call; computing all samples at once needs more
# memory than the VM has for a few hundred EPIC samples
DETP_CHUNK_SIZE <- 50

# --- Optional outputs (qc/2_qc.R) --------------------------------------------

# Outputs no later step reads; large (for ~400 EPICv1 samples about 2.3 GB and
# 4.6 GB). Keep them for ad-hoc analyses, or set FALSE to save space and time
SAVE_BVALS    <- TRUE   # bVals.rds: beta values, same samples/probes as mVals.rds
SAVE_MSET_QC  <- TRUE   # mSetSq_qc_passed.rds (MSET_QC): normalized, filtered GenomicRatioSet

# --- Files passed between qc/ scripts ----------------------------------------

SAMPLE_SHEET   <- file.path(DIR_RESULTS_QC, "sample_sheet.csv")      # 1 -> 2
QC_FIGURE_DATA <- file.path(DIR_RESULTS_QC, "qc_figure_data.rds")    # 2 -> 3: everything needed to redraw the QC figures
