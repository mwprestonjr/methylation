# =============================================================================
# Configuration - ewas module (ewas/ scripts)
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 16, 2026
# Description: Settings used only by the ewas/ scripts: 1_case_control.R,
#              2_age.R, 3_exposure.R, then 4_plots.R. Source after config.R.
#              Inputs are the preprocessing/ outputs (COMBAT_MVALS,
#              SAMPLE_SHEET_FINAL, CPG_ANNOTATION in config.R), plus ancestry
#              from GP2_MASTER_KEY, smoking from CLINICAL_DATA (3_exposure.R)
#              and PPMI_GENETICS (1_case_control.R); outputs go to
#              DIR_RESULTS_EWAS and DIR_FIGURES_EWAS (config.R)
# =============================================================================

# --- Model (all scripts) -----------------------------------------------------

# Every model is fit per CpG with limma on the ComBat-corrected M values:
#   M value ~ variable of interest + EWAS_COVARIATES (+ phenotype, see below)
# Covariates:
#   age, sex
#   cell proportions (CELL_TYPES in config.R, without Neu: they sum to ~1)
#   chip_row: position on the chip (R01-R08), not removed by batch correction
#   ancestry: GP2 genetic ancestry label (below)
# The ComBat batch (combat_batch) is not a covariate: the M values are already
# corrected for it, and adjusting again would spend one degree of freedom per
# batch (Psomagen: one per chip). Categorical covariates with a single value
# among the samples in a model are left out of it
EWAS_COVARIATES <- c("age", "sex", setdiff(CELL_TYPES, "Neu"), "chip_row", "ancestry")

# Ancestry: this GP2_MASTER_KEY (config.R) column, joined by GP2ID
# (load_ewas_inputs() in R/ewas_functions.R). The NBA label covers the most
# samples (the mQTL genotypes are NBA too, GENO_SOURCE in mqtl/config.R).
# Samples without a label are "Unknown" rather than dropped; labels with fewer
# than EWAS_MIN_ANCESTRY_N samples are pooled into "Other"
EWAS_ANCESTRY_COL   <- "nba_label"
EWAS_MIN_ANCESTRY_N <- 10

# The age and exposure models (2, 3) also adjust for phenotype, so a
# PD-control difference isn't picked up as an age or exposure effect (left out
# when all samples have the same phenotype)
EWAS_ADJUST_PHENOTYPE <- TRUE

# Fewest samples with complete data a model may be fit on
EWAS_MIN_SAMPLES <- 10

# Benjamini-Hochberg FDR threshold for calling a CpG significant (summaries)
EWAS_FDR <- 0.05

# Rows in each <model>_top_hits.csv (CpGs with the smallest p-values)
EWAS_TOP_HITS <- 100

# --- Case-control (1_case_control.R) -----------------------------------------

# Phenotypes compared (sample sheet phenotype column); the first is the
# reference, so effects are <second> - <first>. Other phenotypes are left
# out. Skipped if either group has fewer than EWAS_MIN_GROUP_N samples (e.g.
# Psomagen, which has no controls)
EWAS_GROUPS       <- c("Control", "PD")
EWAS_MIN_GROUP_N  <- 3

# Monogenic PD variants (pd_gene, load_pd_gene() in R/pd_genetics.R; PPMI
# only). Where they're known, 1_case_control.R also fits:
#   case_control_no_carriers: EWAS_GROUPS without carriers of any known
#     variant (pd_gene other than EWAS_NONCARRIER; samples without variant
#     data are kept)
#   carrier_<gene>: within EWAS_GROUPS[2] (PD), carriers of each gene in
#     EWAS_CARRIER_GENES vs non-carriers. Carriers of other genes, or of two
#     (LRRK2+GBA1), are left out
EWAS_NONCARRIER    <- "None"
EWAS_CARRIER_GENES <- c("LRRK2", "GBA1")
# Fewest carriers (and non-carriers) a carrier model is fit with. Below this
# the results are unstable: with 8 GBA1 carriers in PPMI, the top CpGs changed
# completely when ancestry was added as a covariate
EWAS_MIN_CARRIER_N <- 10

# --- Age (2_age.R) -----------------------------------------------------------

# Sample sheet column with age at the blood draw (SAMPLE_SHEET_COLUMNS); it's
# the variable of interest here, so it isn't also a covariate
EWAS_AGE_COL <- "age"

# --- Exposure (3_exposure.R) -------------------------------------------------

# Smoking, from each study's own questionnaire items in the GP2 extended
# clinical data (CLINICAL_DATA in config.R), joined by GP2ID (load_smoking()
# in R/ewas_functions.R). smoking_ever (Never, Ever) is regular smoking in
# both data sources; smoking_status (Never, Former, Current) only where
# current smoking is asked. Definitions per questionnaire and how they
# differ: ewas/README.md (Methods and decisions, Smoking)
# Visits used, in order of preference (the first with an answer). In R12 these
# questions are answered at baseline (0) and carried to later visits
EWAS_CLINICAL_VISIT_MONTHS <- c(0, -1.5)
# Per-sample smoking values used by 3_exposure.R (Sample, GP2ID, phenotype,
# smoking_source, smoking_status, smoking_ever)
EWAS_SMOKING_VALUES <- file.path(DIR_RESULTS_EWAS, "smoking_values.csv")

# Exposure tested, and how it's modeled:
#   continuous:  effect per unit
#   categorical: each level vs EWAS_EXPOSURE_REFERENCE, one result set per
#                level (levels with fewer than EWAS_MIN_GROUP_N samples are
#                left out)
# smoking_status or smoking_ever (above), or any other SAMPLE_SHEET_FINAL
# column; the step is skipped if the column isn't available. Per data source:
#   both: smoking_ever (see ewas/README.md for why not smoking_status)
EWAS_EXPOSURE_COL <- switch(DATA_SOURCE,
                            ppmi_p140 = "smoking_ever",
                            psomagen  = "smoking_ever")
EWAS_EXPOSURE_TYPE      <- "categorical"
EWAS_EXPOSURE_REFERENCE <- "Never"

# --- Plots (4_plots.R) -------------------------------------------------------

# Dashed line on the Manhattan plots (suggestive significance)
EWAS_PLOT_SUGGESTIVE_P <- 1e-5

# --- Outputs -----------------------------------------------------------------

# In DIR_RESULTS_EWAS, one set per model, named by the model's prefix
# (case_control, age, exposure_<column>[_<level>]); written by
# write_ewas_result_set() in R/ewas_functions.R:
#   <prefix>_results.csv.gz  every CpG: limma effect (logFC, M value scale),
#                            t, p, FDR, merged with CPG_ANNOTATION
#   <prefix>_top_hits.csv    the EWAS_TOP_HITS CpGs with the smallest p
#   <prefix>_summary.csv     one row: samples, CpGs, FDR hits, genomic inflation
# Figures in DIR_FIGURES_EWAS: <prefix>_{manhattan,qq,distributions}.png
# 4_plots.R
EWAS_MODEL_SUMMARY <- file.path(DIR_RESULTS_EWAS, "ewas_model_summary.csv")  # all <prefix>_summary.csv rows combined
