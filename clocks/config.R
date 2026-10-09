# =============================================================================
# Configuration - clocks module (clocks/ scripts)
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 8, 2026
# Description: Settings used only by the clocks/ scripts. Source after config.R.
#              Outputs go to <DIR_RESULTS>/clocks and <DIR_FIGURES>/clocks
# =============================================================================

# --- Estimation (1_estimate.R) -----------------------------------------------

# Minimum share of a clock's CpGs that must be on the array. methylclock
# leaves CpGs missing from the array out of the weighted sum (no imputation),
# which shifts the clock age (EPICv2 lacks some CpGs used by older clocks)
CLOCK_COV_THRESH <- 0.8    # methylclock (DNAmAge min.perc)
PACE_COV_THRESH  <- 0.8    # DunedinPACE (PACEProjector proportionOfProbesRequired)

# Adult clocks used for age acceleration, comparisons and plots
# (pediatric/telomere clocks are skipped), per array. All clocks are still
# estimated and saved in CLOCK_RESULTS, and their coverage in CLOCK_COVERAGE
#   ppmi_p140 (EPICv1): all six
#   psomagen  (EPICv2): without Hannum and skinHorvath. EPICv2 lacks CpGs they
#                       weight heavily (Hannum: 5 CpGs on EPICv1 but not
#                       EPICv2 add ~15 years to a typical sample), giving
#                       clock ages ~31 and ~10 years below chronological age
#                       and lower correlation with it
ADULT_CLOCKS <- switch(DATA_SOURCE,
                       ppmi_p140 = c("Horvath", "Hannum", "Levine", "skinHorvath", "EN", "BLUP"),
                       psomagen  = c("Horvath", "Levine", "EN", "BLUP"))

# --- Group comparison (2_compare.R; covariates also 3_associations.R) --------

# Phenotypes compared; the first is the reference group. The comparison is
# skipped if either group has fewer than CLOCK_MIN_GROUP_N samples (e.g. a
# data source without controls)
CLOCK_GROUPS      <- c("Control", "PD")
CLOCK_MIN_GROUP_N <- 2

# Covariates in every model (outcome ~ group or trait + covariates):
#   age, sex
#   cell proportions (CELL_TYPES in config.R, without Neu: they sum to ~1)
#   chip_row:     position on the chip (R01-R08), not removed by batch correction
#   combat_batch: the batch variable used by ComBat (plate for PPMI, chip for
#                 Psomagen); the clock betas aren't ComBat-corrected
# Categorical covariates with a single value among the samples in a model are
# left out of it
CLOCK_COVARIATES <- c("age", "sex", setdiff(CELL_TYPES, "Neu"), "chip_row", "combat_batch")

# --- Associations with clinical traits in cases (3_associations.R) ----------

# Cases analysed (phenotype in the sample sheet)
CLOCK_CASE_GROUP <- "PD"

# GP2 extended clinical data (R12): one row per participant and visit
CLINICAL_DATA <- "/mnt/output/metadata/clinical_r12_extended_clinical_data.csv"

# Visits used for the clinical traits, in order of preference: the first visit
# with a value is used. The methylation sample is from the baseline visit
# (visit_month 0); PPMI collects MoCA at screening (-1.5) rather than baseline
CLOCK_CLINICAL_VISIT_MONTHS <- c(0, -1.5)

# PPMI monogenic PD status (PATNO; VAR_GENE: LRRK2, GBA, PRKN, ...)
PPMI_GENETICS <- if (DATA_SOURCE == "ppmi_p140") {
  file.path(SUBJ_DIR, "iu_genetic_consensus_20251025_23Apr2026.csv")
}

# Disease duration = age at the blood draw - this R12 master key age. PPMI has
# age at diagnosis for all cases (age at onset for few); Psomagen only has age
# at onset
CLOCK_DURATION_FROM <- switch(DATA_SOURCE,
                              ppmi_p140 = "age_at_diagnosis",
                              psomagen  = "age_of_onset")

# Traits tested, one model each: <outcome> ~ trait + CLOCK_COVARIATES.
#   trait: a column of CLINICAL_DATA, or one built in R/clock_traits.R
#          (disease_duration, family_history_pd, pd_gene)
#   type:  continuous (effect per unit) or categorical (vs ref)
# Traits without enough data in a data source are skipped (see the minimums
# below), so the list can include traits only one data source has
CLOCK_TRAITS <- data.frame(
  trait = c("disease_duration", "family_history_pd", "pd_gene",
            "hoehn_and_yahr_stage",
            "mds_updrs_part_i_summary_score", "mds_updrs_part_ii_summary_score",
            "mds_updrs_part_iii_summary_score",
            "schwab_england_pct_adl_score", "moca_total_score",
            "dat_sbr_putamen_mean",
            "gds15_total_score", "ess_total_score", "rbd_summary_score"),
  type  = c("continuous", "categorical", "categorical",
            "continuous",
            "continuous", "continuous",
            "continuous",
            "continuous", "continuous",
            "continuous",
            "continuous", "continuous", "continuous"),
  ref   = c(NA, "No", "Idiopathic", rep(NA, 10)),
  label = c("Disease duration (years)", "Family history of PD", "Monogenic PD",
            "Hoehn & Yahr stage",
            "MDS-UPDRS I (non-motor)", "MDS-UPDRS II (motor daily living)",
            "MDS-UPDRS III (motor exam)",
            "Schwab & England ADL (%)", "MoCA",
            "DAT-SPECT putamen SBR",
            "GDS-15 (depression)", "Epworth sleepiness", "RBD questionnaire")
)

CLOCK_MIN_TRAIT_N <- 30   # cases with a value needed to test a trait
CLOCK_MIN_LEVEL_N <- 10   # categorical levels with fewer cases are left out (e.g. GBA1 carriers in PPMI)

# --- Outputs -----------------------------------------------------------------

# <DIR_RESULTS>/clocks and <DIR_FIGURES>/clocks
CLOCK_DIR_RESULTS <- file.path(DIR_RESULTS, "clocks")
CLOCK_DIR_FIGURES <- file.path(DIR_FIGURES, "clocks")
dir.create(CLOCK_DIR_RESULTS, showWarnings = FALSE, recursive = TRUE)
dir.create(CLOCK_DIR_FIGURES, showWarnings = FALSE, recursive = TRUE)

# 1_estimate.R
CLOCK_RESULTS  <- file.path(CLOCK_DIR_RESULTS, "clock_age_acceleration.csv")      # per sample: clock ages, pace, acceleration, covariates
CLOCK_ACCURACY <- file.path(CLOCK_DIR_RESULTS, "clock_accuracy.csv")              # per clock: r, MAE, mean offset vs chronological age
CLOCK_COVERAGE <- file.path(CLOCK_DIR_RESULTS, "clock_cpg_coverage.csv")          # per clock: CpGs on the array, whether estimated
# 2_compare.R
CLOCK_GROUP_EFFECTS <- file.path(CLOCK_DIR_RESULTS, "clock_acceleration_by_diagnosis.csv")  # per outcome: group effect from the models
# 3_associations.R
CLOCK_TRAIT_VALUES  <- file.path(CLOCK_DIR_RESULTS, "clock_case_traits.csv")         # per case: trait values
CLOCK_TRAIT_SUMMARY <- file.path(CLOCK_DIR_RESULTS, "clock_trait_summary.csv")       # per trait: n, distribution, tested or why skipped
CLOCK_TRAIT_EFFECTS <- file.path(CLOCK_DIR_RESULTS, "clock_trait_associations.csv")  # per outcome x trait: effect, CI, p, FDR
