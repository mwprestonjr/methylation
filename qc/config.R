# =============================================================================
# Configuration - QC module (qc/ scripts)
# Author: GP2 Subtypes and Mechanisms - M.E., M.P.
# Date: Oct 7, 2026
# Description: Settings used only by the qc/ scripts. Source after config.R
# =============================================================================

# --- Inputs ------------------------------------------------------------------

# GP2 R12 master key: GP2ID, sex, race, phenotype, age at sample collection
# (sample sheets, qc/01 scripts)
FNAME_METADATA <- "/mnt/output/metadata/INTERNAL_USE_ONLY_master_key_release12_final_vwb.csv"

# --- Sample and probe QC thresholds (qc/02_qc.R) -----------------------------

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

# Samples per detectionP() call; computing all samples at once needs more
# memory than the VM has for a few hundred EPIC samples
DETP_CHUNK_SIZE <- 50

# --- Files passed between qc/ scripts ----------------------------------------

SAMPLE_SHEET   <- file.path(DIR_RESULTS, "sample_sheet.csv")      # 01 -> 02
QC_FIGURE_DATA <- file.path(DIR_RESULTS, "qc_figure_data.rds")    # 02 -> 03: everything needed to redraw the QC figures
