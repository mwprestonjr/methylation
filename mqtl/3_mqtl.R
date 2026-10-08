# =============================================================================
# mQTL Step 3: cis-mQTL mapping
# Author: GP2 Subtypes and Mechanisms - M.E.
# Date: Oct 7, 2026
# Description: Tests SNP-CpG associations within +/- MQTL_CIS_WINDOW using
#              MatrixEQTL, one chromosome at a time. Methylation: ComBat
#              M values, inverse-normal transformed per CpG. Covariates: age,
#              sex, phenotype, cell proportions, chip row, genotype PCs and latent
#              methylation PCs. Probe coordinates come from the array's
#              annotation (R/array_profiles.R) and are lifted hg19 -> hg38
#              when needed (EPICv1) to match GP2 genotypes.
# Usage:       Rscript mqtl/3_mqtl.R <data source> <ancestry>   (see R/mqtl_setup.R)
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

library(minfi)                            # getAnnotation
library(tidyverse)
library(data.table)
library(matrixStats)                      # rowVars
library(MatrixEQTL)
library(GenomicRanges)
library(rtracklayer)                      # liftOver

# Load shared configuration, array profiles and mQTL functions
source("config.R")
source("mqtl/config.R")
source("R/mqtl_setup.R")
source("R/array_profiles.R")
source("R/mqtl_functions.R")

# --- 1. Methylation data and known covariates (R/mqtl_functions.R) -----------

prep     <- mqtl_prepare()
meth_int <- prep$meth_int
cpg_pos  <- prep$cpg_pos

# --- 2. Latent methylation PCs -----------------------------------------------

# MQTL_N_METH_PCS is tuned with mqtl/tune_meth_pcs.R
cat("Computing", MQTL_N_METH_PCS, "latent methylation PCs...\n")
meth_pcs <- latent_meth_pcs(meth_int, prep$known_covs, MQTL_N_METH_PCS)
all_covs <- cbind(prep$known_covs, meth_pcs)
cat("Covariates in model:", ncol(all_covs), "\n")
print(colnames(all_covs))
cvrt <- covariates_sliced(all_covs, colnames(meth_int))

# --- 3. cis-mQTL mapping, per chromosome -------------------------------------

cat("\n--- cis-mQTL mapping (window +/-", MQTL_CIS_WINDOW / 1e6, "Mb) ---\n")

cis_results <- list()
n_tests     <- 0

for (chr in 1:22) {
  cat("\nchr", chr, "\n", sep = "")
  geno <- load_chr_genotypes(prep$geno_dir, chr, prep$sample_map, colnames(meth_int))
  res  <- run_cis_chr(geno, meth_int, cpg_pos, chr, cvrt)
  n_tests <- n_tests + res$ntests
  cat("  cis tests:", res$ntests,
      "| saved (p <", MQTL_P_CIS_SAVE, "):", nrow(res$eqtls), "\n")
  cis_results[[chr]] <- res$eqtls
  rm(geno, res); invisible(gc())
}

# --- 4. Multiple testing and summary -----------------------------------------

cat("\n--- Combining results ---\n")
cat("Total cis tests:", n_tests, "\n")

summ <- summarise_cis(cis_results, n_tests, cpg_pos)
cis  <- summ$cis
lead <- summ$lead

# Significance thresholds (MQTL_FDR, MQTL_P_STRICT in mqtl/config.R)
fdr_label    <- paste0("FDR < ", MQTL_FDR)
strict_label <- paste0("p < ", MQTL_P_STRICT)
sig_fdr    <- lead %>% filter(fdr < MQTL_FDR)
sig_strict <- lead %>% filter(p < MQTL_P_STRICT)
cat("CpGs with a cis-mQTL at", fdr_label, ":", nrow(sig_fdr), "\n")
cat("CpGs with a cis-mQTL at", strict_label, ":", nrow(sig_strict), "\n")

# p-value that the FDR threshold corresponds to; warns if it lies beyond the
# saved range (MQTL_P_CIS_SAVE), i.e. the counts above are incomplete
fdr_p <- fdr_threshold(cis)
if (!is.na(fdr_p)) {
  cat(fdr_label, "threshold: p <=", signif(fdr_p, 3),
      "(within the saved range, p <", MQTL_P_CIS_SAVE, ")\n")
}

# Distance of lead SNPs from their CpG
png(file.path(MQTL_DIR, "mqtl_01_lead_snp_distance.png"),
    width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
print(ggplot(sig_fdr, aes(x = distance / 1e3)) +
  geom_histogram(bins = 100) +
  labs(title = paste0("Lead cis-mQTL SNP position relative to CpG (", fdr_label, ")"),
       x = "SNP - CpG distance (kb)", y = "CpGs") +
  theme_bw())
dev.off()

# --- 5. Save outputs ---------------------------------------------------------

cat("\n--- Saving outputs ---\n")
fwrite(cis,  file.path(MQTL_DIR, "cis_mqtl_all_saved.tsv.gz"), sep = "\t")
fwrite(lead, file.path(MQTL_DIR, "cis_mqtl_lead_per_cpg.tsv"),  sep = "\t")
write.csv(all_covs %>% mutate(meth_id = colnames(meth_int), .before = 1),
          file.path(MQTL_DIR, "mqtl_covariates.csv"), row.names = FALSE)

mqtl_summary <- data.frame(
  metric = c("Samples", "Ancestry", "CpGs tested", "Cis window (bp)",
             "Genotype PCs", "Methylation PCs", "Total cis tests",
             paste0("CpGs with cis-mQTL (", fdr_label, ")"),
             paste0("CpGs with cis-mQTL (", strict_label, ")"),
             "Saved pairs: p <", paste(fdr_label, "p-value threshold")),
  value  = c(ncol(meth_int), paste(MQTL_ANCESTRY, collapse = "+"), nrow(meth_int), MQTL_CIS_WINDOW,
             MQTL_N_GENO_PCS, MQTL_N_METH_PCS, n_tests,
             nrow(sig_fdr), nrow(sig_strict), MQTL_P_CIS_SAVE,
             if (is.na(fdr_p)) "beyond saved range (counts incomplete)" else signif(fdr_p, 3))
)
write.csv(mqtl_summary, file.path(MQTL_DIR, "mqtl_summary.csv"), row.names = FALSE)

cat("\nmQTL step 3 complete!\n")
cat("Outputs saved to:", MQTL_DIR, "\n")
