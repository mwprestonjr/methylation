# =============================================================================
# Building the sample sheet
# Author: GP2 Subtypes and Mechanisms - M.E., M.P.
# Date: Sept 22, 2026
# Updated: Sept 24, 2026
# Description: merges the sample sheets of all datasets in DATASETS, constructs
#              Basename column for minfi, detects array type per chip, merges
#              clinical info from R12
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

library(tidyverse)

# Load shared configuration
source("~/methylation/scripts/00_config.R")

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
r12 <- read.csv(file.path(FNAME_METADATA))
sample_sheet <- sample_sheet %>%
  left_join(r12, by = c("Sample_ID" = "GP2ID"), keep = TRUE)

# keep columns of interest (GP2ID, GP2sampleID, GP2_phenotype, sex, race, age)
sample_sheet <- sample_sheet %>%
  select(GP2ID, GP2sampleID, GP2_phenotype,
         sex  = biological_sex_for_qc,
         race = race_for_qc,
         age  = age_at_sample_collection,
         Dataset, Sentrix_ID, Sentrix_Position, Basename)

# --- 2. Verify idat files exist ----------------------------------------------

cat("\nVerifying idat files exist on disk...\n")
sample_sheet <- sample_sheet %>%
  mutate(
    red_exists = file.exists(paste0(Basename, "_Red.idat")),
    grn_exists = file.exists(paste0(Basename, "_Grn.idat")),
    both_exist = red_exists & grn_exists
  )

cat("Samples with both Red and Green idat files:", 
    sum(sample_sheet$both_exist), "/", nrow(sample_sheet), "\n")

missing <- sample_sheet %>% filter(!both_exist)
if (nrow(missing) > 0) {
  cat("WARNING:", nrow(missing), "samples missing idat files:\n")
  print(missing %>% select(Basename))
} else {
  cat("All idat files found!\n")
}

# --- 2b. Detect array type ---------------------------------------------------

# Arrays can't be told apart from the QC table, so read one idat per chip and
# use the number of bead types (EPICv2: 1,105,209; EPICv1: 1,051,815-1,052,641;
# 450k: 622,399). Anything else (e.g. a genotyping chip) is "Unknown"
cat("\nDetecting array type per chip...\n")
array_from_idat <- function(basename) {
  n_beads <- nrow(illuminaio::readIDAT(paste0(basename, "_Grn.idat"))$Quants)
  case_when(between(n_beads, 1100000, 1110000) ~ "EPICv2",
            between(n_beads, 1045000, 1060000) ~ "EPICv1",
            between(n_beads,  615000,  625000) ~ "450k",
            TRUE                               ~ "Unknown")
}

chip_arrays <- sample_sheet %>%
  filter(both_exist) %>%
  distinct(Sentrix_ID, .keep_all = TRUE) %>%
  transmute(Sentrix_ID, Array = map_chr(Basename, array_from_idat))

sample_sheet <- sample_sheet %>%
  left_join(chip_arrays, by = "Sentrix_ID")

cat("Samples per dataset and array:\n")
print(table(sample_sheet$Dataset, sample_sheet$Array, useNA = "ifany"))

if (any(sample_sheet$Array %in% "Unknown")) {
  cat("WARNING: some chips are not a recognised methylation array",
      "(check DATASETS in 00_config.R):\n")
  print(sample_sheet %>% filter(Array %in% "Unknown") %>% count(Dataset, Sentrix_ID))
}

# --- 3. Final clean sample sheet ---------------------------------------------

# Keep only samples with both idat files
sample_sheet_final <- sample_sheet %>%
  filter(both_exist) %>%
  select(-red_exists, -grn_exists, -both_exist)

cat("\nFinal sample sheet:", nrow(sample_sheet_final), "samples ready for QC\n")

# --- 4. Save sample sheet ----------------------------------------------------

write.csv(sample_sheet_final, SAMPLE_SHEET, row.names = FALSE)
cat("Sample sheet saved to:", SAMPLE_SHEET, "\n")