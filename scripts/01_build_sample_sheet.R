# =============================================================================
# Building the sample sheet
# Author: GP2 Subtypes and Mechanisms - M.E., M.P.
# Date: Sept 22, 2026
# Description: constructs Basename column for minfi, merges clinical info from R12
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

library(tidyverse)

# Load shared configuration
source("~/methylation/scripts/00_config.R")

# --- 1. Load metadata --------------------------------------------------------

# read the sample sheet
cat("Loading target sheet...\n")
sample_sheet <- read.csv(file.path(DIR_DATASET, "AB00000952_QC_table.csv"))

# rename columns
names(sample_sheet)[names(sample_sheet) == "SentrixBarcode_A"]   <- "Sentrix_ID"
names(sample_sheet)[names(sample_sheet) == "SentrixPosition_A"] <- "Sentrix_Position"

# set Basename column for each sample
sample_sheet$Basename <- file.path(
  DIR_DATASET,
  sample_sheet$Sentrix_ID,
  paste0(sample_sheet$Sentrix_ID, "_", sample_sheet$Sentrix_Position),
  paste0(sample_sheet$Sentrix_ID, "_", sample_sheet$Sentrix_Position)
)
cat("Sample sheet dimensions:", nrow(sample_sheet), "rows x", ncol(sample_sheet), "cols\n")

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
         Sentrix_ID, Sentrix_Position, Basename)
  
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

# --- 3. Final clean sample sheet ---------------------------------------------

# Keep only samples with both idat files
sample_sheet_final <- sample_sheet %>%
  filter(both_exist) %>%
  select(-red_exists, -grn_exists, -both_exist)

cat("\nFinal sample sheet:", nrow(sample_sheet_final), "samples ready for QC\n")

# --- 4. Save sample sheet ----------------------------------------------------

write.csv(sample_sheet_final, SAMPLE_SHEET, row.names = FALSE)
cat("Sample sheet saved to:", SAMPLE_SHEET, "\n")