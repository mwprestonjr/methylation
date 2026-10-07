# =============================================================================
# Building the sample sheet - Psomagen deliveries
# Author: GP2 Subtypes and Mechanisms - M.E., M.P.
# Date: Sept 22, 2026
# Updated: Sept 29, 2026
# Description: merges the sample sheets of all datasets in DATASETS, constructs
#              Basename column for minfi, detects array type per chip, merges
#              clinical info from R12, and writes the standard sample sheet
#              (SAMPLE_SHEET_COLUMNS in config.R)
# Usage:       Rscript qc/01_build_sample_sheet_psomagen.R
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

library(tidyverse)

# Load shared configuration, array and idat helpers
DATA_SOURCE <- "psomagen"
source("config.R")
source("qc/config.R")
source("R/array_profiles.R")
source("R/idat_qc.R")

# --- 1. Load metadata --------------------------------------------------------

# read the sample sheet (QC table) of each dataset and set the Basename for minfi
read_dataset_sheet <- function(dataset) {
  dir_dataset <- file.path(DIR_DELIVERY, dataset)
  sheet <- read.csv(file.path(dir_dataset, paste0(dataset, "_QC_table.csv")),
                    colClasses = c(SentrixBarcode_A = "character"))

  sheet %>%
    rename(Sentrix_ID       = SentrixBarcode_A,
           Sentrix_Position = SentrixPosition_A) %>%
    mutate(Dataset  = dataset,
           Basename = file.path(dir_dataset,
                                Sentrix_ID,
                                paste0(Sentrix_ID, "_", Sentrix_Position),
                                paste0(Sentrix_ID, "_", Sentrix_Position)))
}

cat("Loading target sheets for", length(DATASETS), "datasets...\n")
sample_sheet <- map_dfr(DATASETS, read_dataset_sheet)
cat("Merged sample sheet dimensions:", nrow(sample_sheet), "rows x", ncol(sample_sheet), "cols\n")
cat("Samples per dataset:\n")
print(table(sample_sheet$Dataset))

# The same sample can appear in more than one dataset (e.g. re-runs)
dup_ids <- unique(sample_sheet$Sample_ID[duplicated(sample_sheet$Sample_ID)])
if (length(dup_ids) > 0) {
  cat("WARNING:", length(dup_ids), "Sample_IDs appear more than once:\n")
  print(sample_sheet %>% filter(Sample_ID %in% dup_ids) %>%
          select(Sample_ID, Dataset, Sentrix_ID, Sentrix_Position) %>%
          arrange(Sample_ID))
}

# Load R12 and combine (match 'GP2ID' in R12 to 'Sample_ID' in the sample sheet)
r12 <- read.csv(file.path(FNAME_METADATA), colClasses = c(clinical_id = "character"))
sample_sheet <- sample_sheet %>%
  left_join(r12, by = c("Sample_ID" = "GP2ID"), keep = TRUE)

# keep columns of interest, named as in SAMPLE_SHEET_COLUMNS. Methylation was
# measured on the GP2 DNA sample, so age at its collection is the age to use
sample_sheet <- sample_sheet %>%
  transmute(GP2ID, GP2sampleID, clinical_id,
            phenotype = GP2_phenotype,
            sex       = biological_sex_for_qc,
            race      = race_for_qc,
            age       = age_at_sample_collection,
            Dataset,
            Batch     = Dataset,
            Sentrix_ID, Sentrix_Position, Basename)

# --- 2. Verify idat files exist ----------------------------------------------

sample_sheet <- check_idats(sample_sheet)

# --- 2b. Detect array type ---------------------------------------------------

# Arrays can't be told apart from the QC table, so read one idat per chip
# (detect_chip_arrays in R/array_profiles.R)
sample_sheet <- detect_chip_arrays(sample_sheet, chip_col = "Sentrix_ID")

cat("Samples per dataset and array:\n")
print(table(sample_sheet$Dataset, sample_sheet$Array, useNA = "ifany"))

if (any(sample_sheet$Array %in% "Unknown")) {
  cat("WARNING: some chips are not a recognised methylation array",
      "(check DATASETS in config.R):\n")
  print(sample_sheet %>% filter(Array %in% "Unknown") %>% count(Dataset, Sentrix_ID))
}

# --- 3. Final clean sample sheet ---------------------------------------------

# Keep only samples with both idat files
sample_sheet_final <- sample_sheet %>%
  filter(both_exist) %>%
  select(all_of(SAMPLE_SHEET_COLUMNS))

cat("\nFinal sample sheet:", nrow(sample_sheet_final), "samples ready for QC\n")

# --- 4. Save sample sheet ----------------------------------------------------

write.csv(sample_sheet_final, SAMPLE_SHEET, row.names = FALSE)
cat("Sample sheet saved to:", SAMPLE_SHEET, "\n")
