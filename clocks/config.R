# =============================================================================
# Configuration - clocks module (clocks/ scripts)
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 8, 2026
# Description: Settings used only by the clocks/ scripts. Source after config.R
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

# --- Group comparison (2_compare.R) ------------------------------------------

# Phenotypes compared; the first is the reference group. The comparison is
# skipped if either group has fewer than CLOCK_MIN_GROUP_N samples (e.g. a
# data source without controls)
CLOCK_GROUPS      <- c("Control", "PD")
CLOCK_MIN_GROUP_N <- 2

# Covariates in every comparison model (outcome ~ group + covariates):
#   age, sex
#   cell proportions (CELL_TYPES in config.R, without Neu: they sum to ~1)
#   chip_row:     position on the chip (R01-R08), not removed by batch correction
#   combat_batch: the batch variable used by ComBat (plate for PPMI, chip for
#                 Psomagen); the clock betas aren't ComBat-corrected
# Categorical covariates with a single value among the compared samples are
# left out of the model
CLOCK_COVARIATES <- c("age", "sex", setdiff(CELL_TYPES, "Neu"), "chip_row", "combat_batch")

# --- Outputs -----------------------------------------------------------------

# 1_estimate.R
CLOCK_RESULTS  <- file.path(DIR_RESULTS, "clock_age_acceleration.csv")      # per sample: clock ages, pace, acceleration, covariates
CLOCK_ACCURACY <- file.path(DIR_RESULTS, "clock_accuracy.csv")              # per clock: r, MAE, mean offset vs chronological age
CLOCK_COVERAGE <- file.path(DIR_RESULTS, "clock_cpg_coverage.csv")          # per clock: CpGs on the array, whether estimated
# 2_compare.R
CLOCK_GROUP_EFFECTS <- file.path(DIR_RESULTS, "clock_acceleration_by_diagnosis.csv")  # per outcome: group effect from the models
