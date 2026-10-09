# =============================================================================
# mQTL Step 1: Link methylation samples to GP2 genotypes
# Author: GP2 Subtypes and Mechanisms - M.E.
# Date: Oct 7, 2026
# Description: Maps methylation samples to GP2 genotype sample IDs via the
#              GP2 master key (by GP2ID), keeps samples genotyped in
#              GENO_SOURCE that passed GP2 QC, restricts to MQTL_ANCESTRY, and
#              writes a plink2 --keep file for 2_mqtl_genotypes.sh
# Usage:       Rscript mqtl/1_mqtl_sample_map.R <data source> <ancestry>
#              e.g. psomagen AFR, or ppmi_p140 EUR,AJ (see R/mqtl_setup.R)
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

library(tidyverse)

# Load shared configuration and the mQTL ancestry/paths
source("config.R")
source("mqtl/config.R")
source("R/mqtl_setup.R")


# --- 1. Load methylation sample sheet ----------------------------------------

cat("Loading final sample sheet from preprocessing/3_combat.R...\n")
targets <- read.csv(SAMPLE_SHEET_FINAL,
                    colClasses = c(GP2ID       = "character",
                                   clinical_id = "character",
                                   Sentrix_ID  = "character",
                                   Basename    = "character"))
targets <- targets %>%
  mutate(meth_id = basename(Basename))   # matches colnames(combat_mVals)
cat("Methylation samples:", nrow(targets), "\n")

# --- 2. Load GP2 master key --------------------------------------------------

cat("\nLoading GP2 master key...\n")
master_key <- read.csv(GP2_MASTER_KEY, colClasses = "character")

ancestry_col <- paste0(GENO_SOURCE, "_label")
prune_col    <- paste0(GENO_SOURCE, "_prune_reason")

genotyped <- master_key %>%
  filter(.data[[GENO_SOURCE]] == "1.0",                             # genotyped
         is.na(.data[[prune_col]]) | .data[[prune_col]] == "") %>%  # passed GP2 QC
  transmute(GP2ID,
            ancestry = .data[[ancestry_col]]) %>%
  distinct(GP2ID, .keep_all = TRUE)   # one row per person
cat("GP2 participants with", GENO_SOURCE, "genotypes passing QC:", nrow(genotyped), "\n")

# --- 3. Join and restrict to ancestry ----------------------------------------

# Everything is matched on GP2ID: it links the sample sheet, the master key
# and the GP2 genotype files (whose sample IDs are GP2IDs)
sample_map <- targets %>%
  select(GP2ID, clinical_id, meth_id, phenotype) %>%
  inner_join(genotyped, by = "GP2ID")

cat("\nMethylation samples with genotypes:", nrow(sample_map),
    "of", nrow(targets), "\n")
cat("Ancestry breakdown:\n")
print(table(sample_map$ancestry, sample_map$phenotype))

sample_map <- sample_map %>% filter(ancestry %in% MQTL_ANCESTRY)
cat("\nSamples retained (", paste(MQTL_ANCESTRY, collapse = "+"), "):", nrow(sample_map), "\n")
if (nrow(sample_map) == 0) {
  stop("No ", DATA_SOURCE, " samples have ancestry ", paste(MQTL_ANCESTRY, collapse = "+"),
       " - see the ancestry breakdown above")
}

# --- 4. Save outputs ---------------------------------------------------------

write.csv(sample_map, MQTL_SAMPLE_MAP, row.names = FALSE)

# plink2 --keep file, with the sample IDs exactly as in the genotype file's
# .psam (the GP2 files have both FID and IID, each the GP2ID; plink2 matches
# on both when the .psam has FIDs)
psam <- read.delim(paste0(GENO_PFILE, ".psam"), check.names = FALSE,
                   colClasses = "character")
names(psam) <- sub("^#", "", names(psam))
id_cols <- intersect(c("FID", "IID"), names(psam))
keep <- psam[psam$IID %in% sample_map$GP2ID, id_cols, drop = FALSE]
cat("\nSamples in the genotype file (", basename(GENO_PFILE), "):",
    nrow(keep), "of", nrow(sample_map), "\n")
if (nrow(keep) == 0) {
  stop("None of the samples are in ", GENO_PFILE, ".psam - check GENO_PFILE_PATH/GENO_PFILE_NAME")
}
writeLines(c(paste0("#", paste(id_cols, collapse = "\t")),
             do.call(paste, c(keep, sep = "\t"))),
           MQTL_KEEP)

cat("\nSample map saved to:", MQTL_SAMPLE_MAP, "\n")
cat("plink2 keep file saved to:", MQTL_KEEP, "\n")
cat("Next: bash mqtl/2_mqtl_genotypes.sh", DATA_SOURCE, paste(MQTL_ANCESTRY, collapse = ","), "\n")
