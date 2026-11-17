# =============================================================================
# mQTL ancestry and output paths
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 7, 2026
# Updated: Nov 12, 2026
# Description: Sets MQTL_ANCESTRY, the genotype files (GENO_PFILE) and the
#              mQTL output paths for the mqtl/
#              scripts. The ancestry is the second command-line argument: a
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
# ancestry, so several ancestries would need their files merged first
if (!MQTL_GENOTYPES %in% c("wgs", "array")) {
  stop("MQTL_GENOTYPES in mqtl/config.R must be \"wgs\" or \"array\", not \"",
       MQTL_GENOTYPES, "\"")
}
if (length(MQTL_ANCESTRY) > 1) {
  stop("Genotypes are one file set per ancestry; run one ancestry at a time, ",
       "or merge the files and point GENO_PFILE_", toupper(MQTL_GENOTYPES),
       " (mqtl/config.R) at the result")
}
for_ancestry <- function(path) gsub("{ANCESTRY}", MQTL_ANCESTRY, path, fixed = TRUE)
GENO_PFILE <- for_ancestry(switch(MQTL_GENOTYPES, wgs = GENO_PFILE_WGS, array = GENO_PFILE_ARRAY))
cat("Genotypes:", MQTL_GENOTYPES, "-", GENO_PFILE, "\n")

# Outputs per ancestry and genotype set (WGS runs get "_wgs", so the array
# results stay where they were):
#   <DIR_RESULTS_MQTL>/<ancestry>[_wgs]  sample map, genotypes/, cis-mQTL results
#   <DIR_FIGURES_MQTL>/<ancestry>[_wgs]  mQTL figures
MQTL_RUN     <- paste0(MQTL_ANCESTRY, if (MQTL_GENOTYPES == "wgs") "_wgs" else "")
MQTL_DIR     <- file.path(DIR_RESULTS_MQTL, MQTL_RUN)
MQTL_FIG_DIR <- file.path(DIR_FIGURES_MQTL, MQTL_RUN)
dir.create(MQTL_DIR,     showWarnings = FALSE, recursive = TRUE)
dir.create(MQTL_FIG_DIR, showWarnings = FALSE, recursive = TRUE)

MQTL_SAMPLE_MAP <- file.path(MQTL_DIR, "mqtl_sample_map.csv")
MQTL_KEEP       <- file.path(MQTL_DIR, "mqtl_keep.txt")
