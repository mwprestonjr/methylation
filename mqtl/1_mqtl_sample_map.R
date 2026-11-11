# =============================================================================
# mQTL Step 1: Link methylation samples to GP2 genotypes
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 7, 2026
# Updated: Nov 11, 2026
# Description: Maps methylation samples to GP2 genotype sample IDs via the
#              GP2 master key (by GP2ID), keeps samples genotyped in
#              GENO_SOURCE that passed GP2 QC, restricts to MQTL_ANCESTRY, and
#              writes a plink2 --keep file for 2_mqtl_genotypes.sh. With
#              imputed genotypes (MQTL_GENOTYPES), also finds each sample's
#              imputation batch and writes a --keep file per batch
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
# One array per person: two arrays with the same GP2ID would get the same
# genotypes and count as two independent samples. qc/2_qc.R keeps one, so
# this only fails on QC outputs from before that check
dup_gp2 <- unique(sample_map$GP2ID[duplicated(sample_map$GP2ID)])
if (length(dup_gp2) > 0) {
  print(sample_map %>% filter(GP2ID %in% dup_gp2) %>% arrange(GP2ID))
  stop(length(dup_gp2), " GP2IDs have more than one methylation sample (above) - ",
       "rerun qc/2_qc.R ", DATA_SOURCE, ", which keeps one array per person")
}

# --- 4. Imputed genotypes: samples per imputation batch ----------------------

# Each imputation batch has its own VCFs (MQTL_VCF_FILES, R/mqtl_setup.R). The
# VCF sample IDs are the GP2ID, or "<GP2ID>_<GP2ID>" (FID_IID joined when the
# VCFs were written); they're read from the chr22 header. Samples missing
# from the imputed data are left out, so every run step uses the same samples
vcf_sample_ids <- function(vcf) {
  if (!file.exists(vcf)) {
    stop(vcf, " not found - check MQTL_IMPUTED_VCFS in mqtl/config.R and that ",
         "the gp2_release12 bucket is mounted")
  }
  con <- gzfile(vcf, "r")
  on.exit(close(con))
  repeat {
    line <- readLines(con, n = 1)
    if (length(line) == 0) stop("No #CHROM header line in ", vcf)
    if (startsWith(line, "#CHROM")) return(strsplit(line, "\t")[[1]][-(1:9)])
  }
}

if (MQTL_GENOTYPES == "imputed") {
  # The array files have related individuals removed and give the genotype
  # PCs; the imputed VCFs keep relatives, so samples come from the array file
  array_psam <- read.delim(paste0(GENO_PFILE, ".psam"), check.names = FALSE,
                           colClasses = "character")
  array_iids <- array_psam[[which(sub("^#", "", names(array_psam)) == "IID")]]
  cat("\nSamples in the array genotype file (relatives removed):",
      sum(sample_map$GP2ID %in% array_iids), "of", nrow(sample_map), "\n")
  sample_map <- sample_map %>% filter(GP2ID %in% array_iids)

  cat("\nMatching samples to the imputed VCFs (chr22 headers)...\n")
  unlink(Sys.glob(gsub("{PART}", "*", MQTL_KEEP_IMPUTED, fixed = TRUE)))
  sample_map <- sample_map %>% mutate(imputed_part = NA_integer_, imputed_id = NA_character_)
  for (part in seq_along(MQTL_VCF_FILES)) {
    ids <- vcf_sample_ids(gsub("{CHR}", "22", MQTL_VCF_FILES[part], fixed = TRUE))
    hit <- coalesce(match(paste0(sample_map$GP2ID, "_", sample_map$GP2ID), ids),
                    match(sample_map$GP2ID, ids))
    new <- !is.na(hit) & is.na(sample_map$imputed_part)
    sample_map$imputed_part[new] <- part
    sample_map$imputed_id[new]   <- ids[hit[new]]
    cat("  batch", part, ":", length(ids), "samples,", sum(new), "of ours\n")
  }
  not_imputed <- is.na(sample_map$imputed_part)
  if (any(not_imputed)) {
    cat("WARNING:", sum(not_imputed), "samples aren't in the imputed VCFs and are left out\n")
    sample_map <- sample_map %>% filter(!not_imputed)
  }
  if (nrow(sample_map) == 0) {
    stop("None of the samples are in the imputed VCFs - check MQTL_IMPUTED_VCFS ",
         "in mqtl/config.R and the VCF sample ID format (#CHROM line)")
  }

  # plink2 --keep files (VCF IDs; the VCFs have no family IDs)
  for (part in sort(unique(sample_map$imputed_part))) {
    keep_file <- gsub("{PART}", part, MQTL_KEEP_IMPUTED, fixed = TRUE)
    writeLines(c("#IID", sample_map$imputed_id[sample_map$imputed_part == part]), keep_file)
    cat("  keep file for batch", part, ":", keep_file, "\n")
  }
}

# --- 5. Save outputs ---------------------------------------------------------

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
