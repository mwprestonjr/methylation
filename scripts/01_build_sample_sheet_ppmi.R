# =============================================================================
# Building the sample sheet - PPMI Project 140
# Author: GP2 Subtypes and Mechanisms - M.E.
# Date: April 23, 2026
# Updated: Sept 29, 2026
# Description: Filters the link list to one visit (PRIMARY_TIMEPOINT) and the
#              cohorts of interest, merges R12 (GP2ID, sex, race) and age at
#              that visit, constructs the Basename column for minfi, detects
#              array type per chip, and writes the standard sample sheet
#              (SAMPLE_SHEET_COLUMNS in config.R)
# Usage:       Rscript scripts/01_build_sample_sheet_ppmi.R
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

library(tidyverse)

# Load shared configuration, array and idat helpers
DATA_SOURCE <- "ppmi_p140"
source("config.R")
source("R/array_profiles.R")
source("R/idat_qc.R")

# --- 1. Load link list -------------------------------------------------------

# The link list is the only map from idat (Sentrix ID/position) to participant
# (PATNO) and visit (EVENT_ID)
cat("Loading link list...\n")
link_list <- read.csv(LINK_LIST, colClasses = c(SENTRIXID = "character", PATNO = "character"))
cat("Link list dimensions:", nrow(link_list), "rows x", ncol(link_list), "cols\n")
cat("Timepoints available:", paste(unique(link_list$EVENT_ID), collapse = ", "), "\n")

# --- 2. Filter to the primary timepoint --------------------------------------

# 2a. Filter to one visit (baseline by default)
cat("\nFiltering to", PRIMARY_TIMEPOINT, "samples...\n")
visit <- link_list %>%
  filter(EVENT_ID == PRIMARY_TIMEPOINT)
cat(PRIMARY_TIMEPOINT, "samples:", nrow(visit), "\n")

# 2b. Remove duplicates, keep most recent technical replicate
cat("\nChecking for duplicate PATNOs...\n")

# Identify duplicates
duplicate_counts <- visit %>%
  dplyr::group_by(PATNO) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(n > 1)

cat("PATNOs with multiple entries:", nrow(duplicate_counts), "\n")
print(duplicate_counts)

# Remove batch QC control samples (>2 entries)
# 40532, 40535, 40536 are likely batch QC controls or longitudinal samples with mislabeled timepoints
batch_controls <- duplicate_counts %>%
  dplyr::filter(n > 2) %>%
  dplyr::pull(PATNO)

cat("Removing batch QC controls:", paste(batch_controls, collapse = ", "), "\n")
visit <- visit %>%
  filter(!PATNO %in% batch_controls)

# For technical replicates (n=2), keep most recent SENTRIXID
# SENTRIXID assigned sequentially; higher SENTRIXID number = more recent chip
visit <- visit %>%
  group_by(PATNO) %>%
  arrange(desc(SENTRIXID)) %>%  # most recent first
  slice(1) %>%                  # keep first row per PATNO
  ungroup()

cat("Samples after deduplication:", nrow(visit), "\n")
cat("Unique PATNOs:", length(unique(visit$PATNO)), "\n")

# --- 3. Filter to cohorts of interest ----------------------------------------

# PPMI enrollment cohort defines who is included and their phenotype
cat("\nLoading participant status...\n")
participant_status <- read.csv(PARTICIPANT_STATUS,
                               colClasses = c(PATNO = "character")) %>%
  select(PATNO, COHORT, COHORT_DEFINITION, ENROLL_AGE)

sample_sheet <- visit %>%
  left_join(participant_status, by = "PATNO") %>%
  filter(COHORT %in% COHORTS_OF_INTEREST) %>%
  mutate(phenotype = case_when(COHORT == 1 ~ "PD",
                               COHORT == 2 ~ "Control",
                               COHORT == 3 ~ "SWEDD",
                               COHORT == 4 ~ "Prodromal"))

cat("Samples after filtering to cohorts", paste(COHORTS_OF_INTEREST, collapse = ", "), ":",
    nrow(sample_sheet), "\n")
print(table(sample_sheet$phenotype))

# --- 4. Merge R12 (GP2ID, sex, race) and age at visit ------------------------

cat("\nLoading R12...\n")
r12 <- read.csv(FNAME_METADATA, colClasses = c(clinical_id = "character")) %>%
  filter(str_starts(study, "PPMI"), clinical_id != "") %>%
  select(GP2ID, GP2sampleID, clinical_id,
         sex  = biological_sex_for_qc,
         race = race_for_qc)

# clinical_id (PATNO) should be unique within PPMI; keep the first if not
dup_r12 <- unique(r12$clinical_id[duplicated(r12$clinical_id)])
if (length(dup_r12) > 0) {
  cat("WARNING:", length(dup_r12), "PATNOs have more than one R12 row; keeping the first:\n")
  print(r12 %>% filter(clinical_id %in% dup_r12))
  r12 <- r12 %>% distinct(clinical_id, .keep_all = TRUE)
}

sample_sheet <- sample_sheet %>%
  left_join(r12, by = c("PATNO" = "clinical_id"), keep = TRUE)
cat("Samples not in R12:", sum(is.na(sample_sheet$GP2ID)), "\n")
cat("Sex breakdown:\n")
print(table(sample_sheet$sex, useNA = "always"))

# Age at the methylation visit, falling back to age at enrollment if the visit
# has no age (same thing at baseline). R12's age_at_sample_collection is age
# at the GP2 DNA draw, which is often a later visit
age_at_visit <- read.csv(AGE_AT_VISIT, colClasses = c(PATNO = "character")) %>%
  filter(EVENT_ID == PRIMARY_TIMEPOINT) %>%
  distinct(PATNO, .keep_all = TRUE) %>%
  select(PATNO, age = AGE_AT_VISIT)

sample_sheet <- sample_sheet %>%
  left_join(age_at_visit, by = "PATNO")
cat("Samples missing age at", PRIMARY_TIMEPOINT, ":", sum(is.na(sample_sheet$age)), "\n")
if (PRIMARY_TIMEPOINT == "BL") {
  sample_sheet <- sample_sheet %>%
    mutate(age = coalesce(age, as.numeric(ENROLL_AGE)))
  cat("Samples missing age after using ENROLL_AGE:", sum(is.na(sample_sheet$age)), "\n")
}

# --- 5. Map SENTRIXID to actual directory paths ------------------------------

cat("\nMapping SENTRIXID to actual directory paths...\n")

# Find all directories recursively and extract their names
all_dirs <- list.dirs(IDAT_DIR, recursive = TRUE, full.names = TRUE)

# Build a lookup table: SENTRIXID -> full path. The top-level folder under
# IDAT_DIR is the processing plate/delivery, used as the batch
sentrix_map <- tibble(
  full_path = all_dirs,
  SENTRIXID = basename(all_dirs)) %>%
  filter(SENTRIXID %in% sample_sheet$SENTRIXID) %>%
  distinct(SENTRIXID, .keep_all = TRUE) %>%
  mutate(Batch = str_extract(str_remove(full_path, paste0("^", IDAT_DIR, "/")), "^[^/]+"))

cat("Found", nrow(sentrix_map), "matching SENTRIXID directories\n")

# Check for any SENTRIXIDs in our sample sheet not found on disk
missing_sentrix <- sample_sheet %>%
  filter(!SENTRIXID %in% sentrix_map$SENTRIXID) %>%
  pull(SENTRIXID) %>%
  unique()

if (length(missing_sentrix) > 0) {
  cat("WARNING:", length(missing_sentrix),
      "SENTRIXIDs from link list not found on disk:\n")
  print(missing_sentrix)
} else {
  cat("All SENTRIXIDs found on disk!\n")
}

# Join sentrix_map to sample_sheet and construct Basename
sample_sheet <- sample_sheet %>%
  left_join(sentrix_map, by = "SENTRIXID") %>%
  mutate(
    Basename = file.path(full_path, paste0(SENTRIXID, "_", POSITION))
  ) %>%
  select(-full_path)

cat("Basename example:", sample_sheet$Basename[1], "\n")

# --- 6. Verify idat files exist ----------------------------------------------

sample_sheet <- check_idats(sample_sheet, id_cols = c("PATNO", "SENTRIXID", "POSITION"))

# --- 6b. Detect array type ---------------------------------------------------

# One idat per chip (detect_chip_arrays in R/array_profiles.R)
sample_sheet <- detect_chip_arrays(sample_sheet, chip_col = "SENTRIXID")

cat("Samples per array:\n")
print(table(sample_sheet$Array, useNA = "ifany"))

# --- 7. Final clean sample sheet ---------------------------------------------

# Keep only samples with both idat files; standard columns first, then the
# PPMI-specific ones
sample_sheet_final <- sample_sheet %>%
  filter(both_exist) %>%
  mutate(Dataset          = "PPMI_P140",
         Sentrix_ID       = SENTRIXID,
         Sentrix_Position = POSITION) %>%
  select(all_of(SAMPLE_SHEET_COLUMNS), PATNO, EVENT_ID, COHORT, COHORT_DEFINITION)

cat("\nFinal sample sheet:", nrow(sample_sheet_final), "samples ready for QC\n")

# --- 8. Save sample sheet ----------------------------------------------------

write.csv(sample_sheet_final, SAMPLE_SHEET, row.names = FALSE)
cat("Sample sheet saved to:", SAMPLE_SHEET, "\n")
