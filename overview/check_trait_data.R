# =============================================================================
# Overview: clinical trait and exposure data check
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Nov 4, 2026
# Description: Counts, for each data source, how many methylation samples have
#              a known value of each clinical trait and exposure in
#              EWAS_TRAIT_CHECK and EWAS_TRAIT_CHECK_MASTER (below),
#              to decide which traits are worth an EWAS. For each: samples with
#              a value at the first of EWAS_CLINICAL_VISIT_MONTHS (by
#              phenotype), samples with a value at any visit and after
#              baseline, the distribution, and whether there is enough data to
#              test. Fits no models and reads no M values. Not part of
#              `make`. Saves EWAS_TRAIT_CHECK_FILE (counts only)
# Usage:       Rscript overview/check_trait_data.R
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

library(tidyverse)

# config.R needs a data source; both are checked below, by DIR_OUTPUTS
DATA_SOURCE <- "ppmi_p140"
source("config.R")
source("ewas/config.R")

# --- Settings ----------------------------------------------------------------

# Clinical traits and exposures counted, by domain: CLINICAL_DATA columns
# (config.R). Wording: GP2 data dictionary (ewas/README.md, Smoking). Studies
# ask different questionnaires, so several items can measure one exposure
EWAS_TRAIT_CHECK <- list(
  Motor = c("hoehn_and_yahr_stage", "mds_updrs_part_ii_summary_score",
            "mds_updrs_part_iii_summary_score", "mds_updrs_part_iv_summary_score",
            "updrs_part_iii_summary_score", "schwab_england_pct_adl_score"),
  Cognition = c("moca_total_score", "moca_total_score_adjusted", "mmse_total_score",
                "scopa_cog_total_score", "pdq39_cognition_score",
                "clinical_dx_mild_cognitive_impairment", "clinical_dx_dementia",
                "hx_dementia_mci", "dementia_status_at_censoring",
                "age_at_dementia_status_censored"),
  `Mood and sleep` = c("mds_updrs_part_i_summary_score", "gds15_total_score",
                       "depress_test_score", "is_depressed", "hx_depression",
                       "clinical_dx_depression", "merq_depression", "hx_anxiety",
                       "merq_anxiety", "ess_total_score", "rbd_summary_score", "rbd1q",
                       "prodromal_rbd"),
  `Smell, autonomic` = c("smell_test_score", "smell_test_results", "hx_constipation"),
  `Imaging, biofluid` = c("dat_sbr_putamen_mean", "dat_sbr_striatum_mean", "availability_csf"),
  Medication = c("levodopa_use", "dopamine_agonist_use", "levodopa_mg_daily", "ledd_daily",
                 "anti_depressant_use", "medication_data_availability"),
  Smoking = c("smoking_status", "pdrfuq_smoke_100", "pdrfuq_smoke_regularly",
              "merq_cigarettes", "merq_cigarettes_current", "merq_cigarettes_years",
              "merq_cigarettes_daily", "m_merq_smoking_ever", "smoking_years",
              "smoking_pack_years", "merq_smoker_home"),
  Caffeine = c("pdrfuq_coffee", "pdrfuq_black_tea", "pdrfuq_green_tea",
               "pdrfuq_caffeinated_soda", "merq_caffeine", "merq_caffeine_current",
               "m_merq_caffeine_ever", "m_merq_caffeine_per_day"),
  Alcohol = c("pdrfuq_alcohol_100", "pdrfuq_alcohol_regularly", "m_merq_alcohol_ever"),
  `Head injury` = c("pdrfuq_concussion", "merq_head_injury", "m_merq_head_injury_ever",
                    "hx_head_trauma"),
  `Pesticides, chemicals` = c("merq_job_pesticides", "merq_home_pesticides",
                              "m_merq_job_pesticide", "merq_chemical_solvents",
                              "m_merq_job_solvent", "merq_heavy_metals",
                              "merq_chemicals_other", "m_merq_job_weld"),
  `Anti-inflammatory drugs` = c("pdrfuq_ibuprofen", "pdrfuq_aspirin",
                                "pdrfuq_other_anti_inflmmatory_med"),
  `Health, lifestyle` = c("BMI", "education_years", "education_level", "hx_diabetes",
                          "hx_hypertension", "hx_hyperlipidemia", "m_merq_exercise_per_week")
)
# GP2_MASTER_KEY columns counted
EWAS_TRAIT_CHECK_MASTER <- c("age_of_onset", "age_at_diagnosis", "family_history_pd")

# Answers that count as no value
EWAS_TRAIT_MISSING_CODES <- c("", "NA", "Unknown", "Don't Know", "Uncertain", "Refused",
                              "Prefer not to answer", "Not applicable", "-9", "-7")
# Enough data to test: at least EWAS_TRAIT_MIN_N samples and, for a categorical
# item, 2 levels with EWAS_TRAIT_MIN_LEVEL_N samples (as CLOCK_MIN_TRAIT_N and
# CLOCK_MIN_LEVEL_N in clocks/config.R). Numbers with up to
# EWAS_TRAIT_MAX_LEVELS distinct values count as categorical
EWAS_TRAIT_MIN_N       <- 30
EWAS_TRAIT_MIN_LEVEL_N <- 10
EWAS_TRAIT_MAX_LEVELS  <- 6
# Both data sources in one table, so next to their output folders
EWAS_TRAIT_CHECK_FILE <- file.path(dirname(DIR_OUTPUT), "reports", "ewas_trait_data_check.csv")

for (f in c(CLINICAL_DATA, GP2_MASTER_KEY)) {
  if (!file.exists(f)) stop(f, " not found (config.R) - is the gp2_release12 bucket mounted?")
}

items <- bind_rows(
  imap_dfr(EWAS_TRAIT_CHECK, ~ tibble(domain = .y, item = .x, file = "clinical")),
  tibble(domain = "Master key", item = EWAS_TRAIT_CHECK_MASTER, file = "master")
)

# --- 1. Samples --------------------------------------------------------------

# One row per person with a GP2ID in the final sample sheet of a data source
load_samples <- function(source) {
  f <- file.path(DIR_OUTPUTS[[source]], "results", "preprocessing", basename(SAMPLE_SHEET_FINAL))
  if (!file.exists(f)) stop(f, " not found - run make preprocessing SOURCE=", source)
  sheet <- read.csv(f, colClasses = c(GP2ID = "character"))
  with_id <- sheet %>% filter(!is.na(GP2ID), GP2ID != "") %>% distinct(GP2ID, phenotype)
  cat(source, ":", nrow(sheet), "samples,", nrow(with_id), "with a GP2ID\n")
  with_id
}

cat("Loading the final sample sheets...\n")
samples <- map(set_names(DATA_SOURCES), load_samples)
all_ids <- unlist(map(samples, "GP2ID"))

# --- 2. Clinical data and master key -----------------------------------------

cat("\nLoading", CLINICAL_DATA, "...\n")
clinical <- readr::read_csv(CLINICAL_DATA,
                            col_types  = readr::cols(.default = "c"),
                            col_select = any_of(c("GP2ID", "visit_month",
                                                  items$item[items$file == "clinical"])),
                            progress   = FALSE) %>%
  filter(GP2ID %in% all_ids) %>%
  mutate(visit_month = as.numeric(visit_month))

# The master key has one row per person, so its values count as baseline
cat("Loading", GP2_MASTER_KEY, "...\n")
master <- read.csv(GP2_MASTER_KEY, colClasses = "character", na.strings = c("", "NA")) %>%
  filter(GP2ID %in% all_ids) %>%
  distinct(GP2ID, .keep_all = TRUE) %>%
  select(GP2ID, any_of(items$item[items$file == "master"])) %>%
  mutate(visit_month = EWAS_CLINICAL_VISIT_MONTHS[1])

# --- 3. Count one item in one data source ------------------------------------

is_known <- function(x) {
  !is.na(x) & !tolower(trimws(x)) %in% tolower(EWAS_TRAIT_MISSING_CODES)
}

# Distribution of the known values, and whether they can be tested. Numbers
# with few distinct values (e.g. 0/1, H&Y stage) count as categories
describe_values <- function(values) {
  counts  <- sort(table(values), decreasing = TRUE)
  numbers <- suppressWarnings(as.numeric(values))
  if (length(values) == 0) {
    return(list(type = NA, distribution = NA, status = "no values"))
  }
  if (!anyNA(numbers) && length(counts) > EWAS_TRAIT_MAX_LEVELS) {
    return(list(type = "numeric",
                distribution = sprintf("mean %.3g, SD %.3g, range %.3g to %.3g", mean(numbers),
                                       sd(numbers), min(numbers), max(numbers)),
                status = if (length(values) < EWAS_TRAIT_MIN_N) {
                  paste("fewer than", EWAS_TRAIT_MIN_N, "samples")
                } else "testable"))
  }
  shown <- head(counts, EWAS_TRAIT_MAX_LEVELS)
  list(type = "categorical",
       distribution = paste0(paste(names(shown), shown, collapse = ", "),
                             if (length(counts) > length(shown)) ", ..."),
       status = if (length(values) < EWAS_TRAIT_MIN_N) {
         paste("fewer than", EWAS_TRAIT_MIN_N, "samples")
       } else if (sum(counts >= EWAS_TRAIT_MIN_LEVEL_N) < 2) {
         paste("fewer than 2 levels with", EWAS_TRAIT_MIN_LEVEL_N, "samples")
       } else "testable")
}

check_item <- function(source, domain, item, file) {
  row  <- tibble(data_source = source, domain = domain, item = item, file = file,
                 n_samples = nrow(samples[[source]]))
  data <- if (file == "master") master else clinical
  if (!item %in% names(data)) {
    return(row %>% mutate(n_known = 0L, status = "not in file"))
  }

  answers <- data %>%
    filter(GP2ID %in% samples[[source]]$GP2ID) %>%
    transmute(GP2ID, visit_month, value = .data[[item]]) %>%
    filter(is_known(value))

  # Baseline: the first visit of EWAS_CLINICAL_VISIT_MONTHS with a value
  baseline <- answers %>%
    mutate(visit_rank = match(visit_month, EWAS_CLINICAL_VISIT_MONTHS)) %>%
    filter(!is.na(visit_rank)) %>%
    arrange(GP2ID, visit_rank) %>%
    distinct(GP2ID, .keep_all = TRUE) %>%
    left_join(samples[[source]], by = "GP2ID")

  described <- describe_values(baseline$value)
  row %>%
    mutate(n_known      = nrow(baseline),
           pct_known    = round(100 * n_known / n_samples, 1),
           n_case       = sum(baseline$phenotype == EWAS_GROUPS[2]),
           n_control    = sum(baseline$phenotype == EWAS_GROUPS[1]),
           n_any_visit  = n_distinct(answers$GP2ID),
           n_after_baseline = n_distinct(answers$GP2ID[answers$visit_month > 0]),
           type         = described$type,
           distribution = described$distribution,
           status       = described$status)
}

# --- 4. Count every item -----------------------------------------------------

cat("\nCounting", nrow(items), "items...\n")
coverage <- map_dfr(DATA_SOURCES, function(source) {
  pmap_dfr(items, function(domain, item, file) check_item(source, domain, item, file))
})

# --- 5. Summary and save -----------------------------------------------------

# Known at baseline per data source, side by side
wide <- coverage %>%
  mutate(cell = ifelse(status == "testable", as.character(n_known),
                       ifelse(status == "not in file", "-", paste0(n_known, " (too few)")))) %>%
  select(domain, item, data_source, cell) %>%
  pivot_wider(names_from = data_source, values_from = cell)

cat("\nSamples with a known value at baseline (-: not in the file; too few:",
    "fewer than", EWAS_TRAIT_MIN_N, "samples, or fewer than 2 levels with",
    EWAS_TRAIT_MIN_LEVEL_N, "):\n")
print(as.data.frame(wide), row.names = FALSE, right = FALSE)

cat("\nTestable items per data source:\n")
print(coverage %>% filter(status == "testable") %>% count(data_source))

dir.create(dirname(EWAS_TRAIT_CHECK_FILE), showWarnings = FALSE, recursive = TRUE)
write.csv(coverage, EWAS_TRAIT_CHECK_FILE, row.names = FALSE)
cat("\nFull table (by phenotype, later visits, distributions) saved to:",
    EWAS_TRAIT_CHECK_FILE, "\n")
