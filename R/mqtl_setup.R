# =============================================================================
# mQTL ancestry and output paths
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 7, 2026
# Updated: Nov 17, 2026
# Description: Sets MQTL_ANCESTRY, the genotype files (GENO_PFILE, and for
#              imputed runs MQTL_VCF_FILES) and the mQTL output paths for the
#              mqtl/ scripts. The ancestry is the second command-line argument: a
#              GP2 master key label (wgs_label / nba_label), e.g.
#                Rscript mqtl/1_mqtl_sample_map.R psomagen AFR
#              One ancestry at a time: GP2 genotypes are one file set per
#              ancestry. In an interactive session, set MQTL_ANCESTRY before
#              sourcing. Each ancestry and genotype set (MQTL_GENOTYPES) gets
#              its own folders in DIR_RESULTS_MQTL and DIR_FIGURES_MQTL, so
#              runs don't overwrite each other.
#              Source after config.R and mqtl/config.R
# =============================================================================

if (!exists("MQTL_ANCESTRY")) MQTL_ANCESTRY <- commandArgs(trailingOnly = TRUE)[2]
if (length(MQTL_ANCESTRY) == 0 || is.na(MQTL_ANCESTRY[1]) || MQTL_ANCESTRY[1] == "") {
  stop("Set the mQTL ancestry (GP2 master key label, e.g. EUR or AFR) as the ",
       "second argument, e.g. Rscript mqtl/1_mqtl_sample_map.R ", DATA_SOURCE, " EUR")
}
MQTL_ANCESTRY <- unlist(strsplit(MQTL_ANCESTRY, ","))
cat("mQTL ancestry:", paste(MQTL_ANCESTRY, collapse = "+"), "\n")

# Genotype file prefix for the ancestry (MQTL_GENOTYPES, GENO_PFILE_WGS /
# GENO_PFILE_ARRAY in mqtl/config.R). Genotypes are released one file set per
# ancestry, so several ancestries would need their files merged first.
# Imputed runs use the array files for the samples, relatives and genotype PCs
if (!MQTL_GENOTYPES %in% c("wgs", "array", "imputed")) {
  stop("MQTL_GENOTYPES in mqtl/config.R must be \"wgs\", \"array\" or \"imputed\", not \"",
       MQTL_GENOTYPES, "\"")
}
if (length(MQTL_ANCESTRY) > 1) {
  stop("Genotypes are one file set per ancestry; run one ancestry at a time, ",
       "or merge the files and point GENO_PFILE_WGS / GENO_PFILE_ARRAY ",
       "(mqtl/config.R) at the result")
}
for_ancestry <- function(path) gsub("{ANCESTRY}", MQTL_ANCESTRY, path, fixed = TRUE)
GENO_PFILE <- for_ancestry(switch(MQTL_GENOTYPES, wgs = GENO_PFILE_WGS,
                                  array = GENO_PFILE_ARRAY, imputed = GENO_PFILE_ARRAY))
cat("Genotypes:", MQTL_GENOTYPES, "-", GENO_PFILE, "\n")

# Imputed dosage VCFs for the ancestry (MQTL_IMPUTED_VCFS in mqtl/config.R):
# one path pattern per imputation batch, with {CHR} for the chromosome.
# MQTL_VCF_PATTERNS is the same, space-separated, for 2_mqtl_genotypes.sh
MQTL_VCF_FILES <- character(0)
if (MQTL_GENOTYPES == "imputed") {
  MQTL_VCF_FILES <- if (MQTL_ANCESTRY %in% names(MQTL_IMPUTED_VCFS)) {
    MQTL_IMPUTED_VCFS[[MQTL_ANCESTRY]]
  } else {
    for_ancestry(MQTL_IMPUTED_VCFS$default)
  }
  MQTL_VCF_FILES <- path.expand(MQTL_VCF_FILES)
  cat("Imputed dosages:", length(MQTL_VCF_FILES), "imputation batch(es)\n")
}
MQTL_VCF_PATTERNS <- paste(MQTL_VCF_FILES, collapse = " ")

# Outputs per ancestry and genotype set (WGS and imputed runs get "_wgs" /
# "_imputed", so the array results stay where they were):
#   <DIR_RESULTS_MQTL>/<ancestry>[_wgs|_imputed]  sample map, genotypes/, cis-mQTL results
#   <DIR_FIGURES_MQTL>/<ancestry>[_wgs|_imputed]  mQTL figures
MQTL_RUN     <- paste0(MQTL_ANCESTRY, if (MQTL_GENOTYPES == "array") "" else paste0("_", MQTL_GENOTYPES))
MQTL_DIR     <- file.path(DIR_RESULTS_MQTL, MQTL_RUN)
MQTL_FIG_DIR <- file.path(DIR_FIGURES_MQTL, MQTL_RUN)
dir.create(MQTL_DIR,     showWarnings = FALSE, recursive = TRUE)
dir.create(MQTL_FIG_DIR, showWarnings = FALSE, recursive = TRUE)

MQTL_SAMPLE_MAP <- file.path(MQTL_DIR, "mqtl_sample_map.csv")
MQTL_KEEP       <- file.path(MQTL_DIR, "mqtl_keep.txt")
# plink2 --keep file per imputation batch (VCF sample IDs), written by step 1
# for the batches holding any of the samples; {PART} is the batch number
MQTL_KEEP_IMPUTED <- file.path(MQTL_DIR, "mqtl_keep_imputed_part{PART}.txt")
