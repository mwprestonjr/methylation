# =============================================================================
# mQTL ancestry and output paths
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 7, 2026
# Description: Sets MQTL_ANCESTRY, the genotype files (GENO_PFILE) and the
#              mQTL output paths for the mqtl/
#              scripts. The ancestry is the second command-line argument: a
#              GP2 master key label (wgs_label / nba_label), e.g.
#                Rscript mqtl/1_mqtl_sample_map.R psomagen AFR
#              Several labels separated by commas (e.g. EUR,AJ) only work
#              with a single genotype file covering all of them (see
#              GENO_PFILE_NAME in mqtl/config.R)
#              In an interactive session, set MQTL_ANCESTRY before sourcing.
#              Each ancestry (set) gets its own output folder, so runs for
#              different ancestries don't overwrite each other.
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

MQTL_DIR        <- file.path(DIR_RESULTS, "mqtl", paste(MQTL_ANCESTRY, collapse = "_"))
MQTL_SAMPLE_MAP <- file.path(MQTL_DIR, "mqtl_sample_map.csv")
MQTL_KEEP       <- file.path(MQTL_DIR, "mqtl_keep.txt")
