# =============================================================================
# mQTL ancestry and output paths
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 7, 2026
# Updated: Nov 11, 2026
# Description: Sets MQTL_ANCESTRY, the genotype files (GENO_PFILE, and the
#              imputed MQTL_VCF_FILES) and the mQTL output paths for the mqtl/
#              scripts. The ancestry is the second command-line argument: a
#              GP2 master key label (wgs_label / nba_label), e.g.
#                Rscript mqtl/1_mqtl_sample_map.R psomagen AFR
#              Several labels separated by commas (e.g. EUR,AJ) only work
#              with a single genotype file covering all of them (see
#              GENO_PFILE_NAME in mqtl/config.R)
#              In an interactive session, set MQTL_ANCESTRY before sourcing.
#              Each ancestry (set) gets its own folders in DIR_RESULTS_MQTL
#              and DIR_FIGURES_MQTL, so runs for different ancestries don't
#              overwrite each other.
#              Source after config.R and mqtl/config.R
# =============================================================================

if (!exists("MQTL_ANCESTRY")) MQTL_ANCESTRY <- commandArgs(trailingOnly = TRUE)[2]
if (length(MQTL_ANCESTRY) == 0 || is.na(MQTL_ANCESTRY[1]) || MQTL_ANCESTRY[1] == "") {
  stop("Set the mQTL ancestry (GP2 master key label, e.g. EUR or AFR) as the ",
       "second argument, e.g. Rscript mqtl/1_mqtl_sample_map.R ", DATA_SOURCE, " EUR")
}
MQTL_ANCESTRY <- unlist(strsplit(MQTL_ANCESTRY, ","))
cat("mQTL ancestry:", paste(MQTL_ANCESTRY, collapse = "+"), "\n")

# Genotype file prefix for the ancestry (GENO_PFILE_PATH/GENO_PFILE_NAME in
# mqtl/config.R). Genotypes are released one file set per ancestry, so several
# ancestries would need their files merged first
if (grepl("{ANCESTRY}", GENO_PFILE_NAME, fixed = TRUE) && length(MQTL_ANCESTRY) > 1) {
  stop("Genotypes are one file per ancestry (GENO_PFILE_NAME); run one ancestry ",
       "at a time, or merge the files and point GENO_PFILE_NAME at the result")
}
GENO_PFILE <- file.path(GENO_PFILE_PATH,
                        gsub("{ANCESTRY}", MQTL_ANCESTRY, GENO_PFILE_NAME, fixed = TRUE))

# Imputed dosage VCFs for the ancestry (MQTL_IMPUTED_VCFS in mqtl/config.R):
# one path pattern per imputation batch, with {CHR} for the chromosome.
# MQTL_VCF_PATTERNS is the same, space-separated, for 2_mqtl_genotypes.sh
if (!MQTL_GENOTYPES %in% c("imputed", "array")) {
  stop("MQTL_GENOTYPES in mqtl/config.R must be \"imputed\" or \"array\", not \"",
       MQTL_GENOTYPES, "\"")
}
MQTL_VCF_FILES <- character(0)
if (MQTL_GENOTYPES == "imputed") {
  if (length(MQTL_ANCESTRY) > 1) {
    stop("Imputed genotypes are one set of files per ancestry; run one ancestry ",
         "at a time, or set MQTL_GENOTYPES <- \"array\" in mqtl/config.R")
  }
  MQTL_VCF_FILES <- if (MQTL_ANCESTRY %in% names(MQTL_IMPUTED_VCFS)) {
    MQTL_IMPUTED_VCFS[[MQTL_ANCESTRY]]
  } else {
    gsub("{ANCESTRY}", MQTL_ANCESTRY, MQTL_IMPUTED_VCFS$default, fixed = TRUE)
  }
  MQTL_VCF_FILES <- path.expand(MQTL_VCF_FILES)
  cat("Genotypes: imputed,", length(MQTL_VCF_FILES), "imputation batch(es)\n")
} else {
  cat("Genotypes: array\n")
}
MQTL_VCF_PATTERNS <- paste(MQTL_VCF_FILES, collapse = " ")

# Outputs per ancestry, and per genotype set (imputed runs get "_imputed", so
# the array results stay where they were):
#   <DIR_RESULTS_MQTL>/<ancestry>[_imputed]  sample map, genotypes/, cis-mQTL results
#   <DIR_FIGURES_MQTL>/<ancestry>[_imputed]  mQTL figures
MQTL_RUN     <- paste0(paste(MQTL_ANCESTRY, collapse = "_"),
                       if (MQTL_GENOTYPES == "imputed") "_imputed" else "")
MQTL_DIR     <- file.path(DIR_RESULTS_MQTL, MQTL_RUN)
MQTL_FIG_DIR <- file.path(DIR_FIGURES_MQTL, MQTL_RUN)
dir.create(MQTL_DIR,     showWarnings = FALSE, recursive = TRUE)
dir.create(MQTL_FIG_DIR, showWarnings = FALSE, recursive = TRUE)

MQTL_SAMPLE_MAP <- file.path(MQTL_DIR, "mqtl_sample_map.csv")
MQTL_KEEP       <- file.path(MQTL_DIR, "mqtl_keep.txt")
# plink2 --keep file per imputation batch (VCF sample IDs), written by step 1
# for the batches holding any of the samples; {PART} is the batch number
MQTL_KEEP_IMPUTED <- file.path(MQTL_DIR, "mqtl_keep_imputed_part{PART}.txt")
