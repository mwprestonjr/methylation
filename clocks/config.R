# =============================================================================
# Configuration - clocks module (clocks/ scripts)
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 8, 2026
# Description: Settings used only by the clocks/ scripts. Source after config.R.
#              Outputs go to DIR_RESULTS_CLOCKS and DIR_FIGURES_CLOCKS (config.R)
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
#   psomagen  (EPICv2): all six. EPICv2 lacks CpGs that Hannum and skinHorvath
#                       weight heavily (Hannum: 5 CpGs on EPICv1 but not
#                       EPICv2 add ~15 years to a typical sample), so their
#                       clock ages are ~31 and ~10 years below chronological
#                       age. 5_epicv2_simulation.R shows this is a near-constant
#                       offset: with EPICv2 CpGs only, PPMI's correlation with
#                       age is unchanged and age acceleration correlates
#                       r = 0.99 with the full-array values
ADULT_CLOCKS <- switch(DATA_SOURCE,
                       ppmi_p140 = c("Horvath", "Hannum", "Levine", "skinHorvath", "EN", "BLUP"),
                       psomagen  = c("Horvath", "Hannum", "Levine", "skinHorvath", "EN", "BLUP"))

# --- Group comparison (2_compare.R; covariates also 3_associations.R) --------

# FDR (Benjamini-Hochberg) threshold for calling results significant, in
# 2_compare.R (across outcomes), 3_associations.R (across outcome x trait
# tests) and the 4_plots.R figures
CLOCK_FDR <- 0.05

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

# Sensitivity analysis (2_compare.R, 3_associations.R): the same models without
# the cell proportions. Part of what the extrinsic clocks (Hannum, Levine)
# measure is immune cell composition, which adjusting for cells removes; this
# shows whether a difference is driven by cell composition
CLOCK_COVARIATES_NO_CELLS <- setdiff(CLOCK_COVARIATES, CELL_TYPES)

# --- Associations with clinical traits in cases (3_associations.R) ----------

# Cases analyzed (phenotype in the sample sheet)
CLOCK_CASE_GROUP <- "PD"

# Clinical traits come from the GP2 extended clinical data (CLINICAL_DATA in
# config.R; also used by ewas/)

# Visits used for the clinical traits, in order of preference: the first visit
# with a value is used. The methylation sample is from the baseline visit
# (visit_month 0); PPMI collects MoCA at screening (-1.5) rather than baseline
CLOCK_CLINICAL_VISIT_MONTHS <- c(0, -1.5)

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

# --- EPICv2 simulation (5_epicv2_simulation.R) -------------------------------

# An EPICv1 data source is re-estimated with only the CpGs on EPICv2 (from the
# EPICv2 annotation package, ARRAY_PROFILES in R/array_profiles.R) and compared
# with its full-array estimates and with the EPICv2 data source below, to
# separate the effect of the array's missing CpGs from other differences
# between the cohorts (population, age precision, disease stage)
EPICV2_SIM_SOURCE     <- "psomagen"
EPICV2_SIM_COMPARISON <- file.path(DIR_OUTPUTS[[EPICV2_SIM_SOURCE]],
                                   "results", "clocks")   # its DIR_RESULTS_CLOCKS (config.R)

# --- Outputs -----------------------------------------------------------------

# In DIR_RESULTS_CLOCKS; figures in DIR_FIGURES_CLOCKS (config.R)
# 1_estimate.R
CLOCK_RESULTS  <- file.path(DIR_RESULTS_CLOCKS, "clock_age_acceleration.csv")      # per sample: clock ages, pace, acceleration, covariates
CLOCK_ACCURACY <- file.path(DIR_RESULTS_CLOCKS, "clock_accuracy.csv")              # per clock: r, MAE, mean offset vs chronological age
CLOCK_COVERAGE <- file.path(DIR_RESULTS_CLOCKS, "clock_cpg_coverage.csv")          # per clock: CpGs on the array, whether estimated
# 2_compare.R
CLOCK_GROUP_EFFECTS <- file.path(DIR_RESULTS_CLOCKS, "clock_acceleration_by_diagnosis.csv")  # per outcome: group effect from the models
CLOCK_GROUP_EFFECTS_NO_CELLS <- file.path(DIR_RESULTS_CLOCKS, "clock_acceleration_by_diagnosis_no_cells.csv")  # same, sensitivity model without cell proportions
# 3_associations.R
CLOCK_TRAIT_VALUES  <- file.path(DIR_RESULTS_CLOCKS, "clock_case_traits.csv")         # per case: trait values
CLOCK_TRAIT_SUMMARY <- file.path(DIR_RESULTS_CLOCKS, "clock_trait_summary.csv")       # per trait: n, distribution, tested or why skipped
CLOCK_TRAIT_EFFECTS <- file.path(DIR_RESULTS_CLOCKS, "clock_trait_associations.csv")  # per outcome x trait: effect, CI, p, FDR
CLOCK_TRAIT_EFFECTS_NO_CELLS <- file.path(DIR_RESULTS_CLOCKS, "clock_trait_associations_no_cells.csv")  # same, sensitivity model without cell proportions
# 5_epicv2_simulation.R
EPICV2_SIM_RESULTS       <- file.path(DIR_RESULTS_CLOCKS, "clock_epicv2_sim_age_acceleration.csv")  # per sample: clock ages, pace, acceleration from EPICv2 CpGs only
EPICV2_SIM_COVERAGE      <- file.path(DIR_RESULTS_CLOCKS, "clock_epicv2_sim_cpg_coverage.csv")      # per clock: CpGs left with EPICv2 CpGs only
EPICV2_SIM_ACCURACY      <- file.path(DIR_RESULTS_CLOCKS, "clock_epicv2_sim_accuracy.csv")          # per clock x setting: r, MAE, offset vs chronological age
EPICV2_SIM_AGREEMENT     <- file.path(DIR_RESULTS_CLOCKS, "clock_epicv2_sim_agreement.csv")         # per outcome: EPICv2-CpG vs full-array estimates
EPICV2_SIM_GROUP_EFFECTS <- file.path(DIR_RESULTS_CLOCKS, "clock_epicv2_sim_by_diagnosis.csv")      # per outcome: group effect from the EPICv2-CpG estimates
