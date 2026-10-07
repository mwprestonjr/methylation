# =============================================================================
# mQTL: tune the number of latent methylation PCs
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 7, 2026
# Description: Runs the cis-mQTL scan of 3_mqtl.R once for each number of
#              latent methylation PCs in MQTL_N_METH_PCS_GRID and counts the
#              CpGs with a cis-mQTL for each. The usual choice for
#              MQTL_N_METH_PCS is where this curve levels off. The latent PCs
#              are nested, so they are computed once, and each chromosome's
#              genotypes are loaded once for all values. Needs steps 1 and 2
#              (sample map, genotypes); doesn't change the step 3 outputs.
#              Saves tune_meth_pcs.csv and mqtl_tune_meth_pcs.png in MQTL_DIR
# Usage:       Rscript mqtl/tune_meth_pcs.R <data source> <ancestry> [values]
#              e.g. psomagen AFR, or psomagen AFR 0,5,10 to override the grid
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

# Values to test: third command-line argument (comma-separated), or
# MQTL_N_METH_PCS_GRID in mqtl/config.R
grid_arg <- commandArgs(trailingOnly = TRUE)[3]
grid <- if (!is.na(grid_arg)) as.integer(strsplit(grid_arg, ",")[[1]]) else MQTL_N_METH_PCS_GRID
grid <- sort(unique(grid))
cat("Numbers of latent methylation PCs to test:", paste(grid, collapse = ", "), "\n")

# --- 1. Methylation data, covariates and latent PCs --------------------------

prep     <- mqtl_prepare()
meth_int <- prep$meth_int
cpg_pos  <- prep$cpg_pos

# The first k of max(grid) components are the k-component PCs, so compute once
cat("Computing", max(grid), "latent methylation PCs...\n")
pcs_max <- latent_meth_pcs(meth_int, prep$known_covs, max(grid))
cvrts <- lapply(grid, function(k) {
  covariates_sliced(cbind(prep$known_covs, pcs_max[, seq_len(k), drop = FALSE]),
                    colnames(meth_int))
})
names(cvrts) <- grid

# --- 2. cis-mQTL scan for every value, per chromosome ------------------------

cat("\n--- cis-mQTL mapping (window +/-", MQTL_CIS_WINDOW / 1e6, "Mb) ---\n")
cis_results <- setNames(rep(list(list()), length(grid)), grid)
n_tests     <- 0
start_time  <- Sys.time()

for (chr in 1:22) {
  cat("chr", chr, ":", sep = "")
  geno <- load_chr_genotypes(prep$geno_dir, chr, prep$sample_map, colnames(meth_int))
  for (k in as.character(grid)) {
    res <- run_cis_chr(geno, meth_int, cpg_pos, chr, cvrts[[k]])
    cis_results[[k]][[chr]] <- res$eqtls
    cat(" ", k, sep = "")
  }
  n_tests <- n_tests + res$ntests   # the same tests for every value
  cat("\n")
  rm(geno, res); invisible(gc())
}
cat("Scan time:", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1), "min\n")

# --- 3. Counts per value -----------------------------------------------------

tuning <- bind_rows(lapply(as.character(grid), function(k) {
  summ <- summarise_cis(cis_results[[k]], n_tests, cpg_pos)
  lead <- summ$lead
  # NA (with a warning) if the FDR threshold is beyond the saved range
  fdr_p <- fdr_threshold(summ$cis, label = paste0(k, " PCs: "))
  data.frame(n_meth_pcs       = as.integer(k),
             n_covariates     = ncol(prep$known_covs) + as.integer(k),
             cpgs_fdr05       = sum(lead$fdr < 0.05),
             cpgs_p1e8        = sum(lead$p < 1e-8),
             fdr05_p_threshold = fdr_p)
}))
tuning$n_tests <- n_tests

# Where the curve levels off: the smallest value within 1% of the maximum
# number of CpGs at FDR < 5%
best      <- tuning$n_meth_pcs[which.max(tuning$cpgs_fdr05)]
plateau   <- min(tuning$n_meth_pcs[tuning$cpgs_fdr05 >= 0.99 * max(tuning$cpgs_fdr05)])
tuning$most_cpgs   <- tuning$n_meth_pcs == best
tuning$plateau_99  <- tuning$n_meth_pcs == plateau
print(tuning)
cat("\nMost CpGs at FDR < 5%:", best, "PCs | smallest value within 1% of that:",
    plateau, "PCs | current MQTL_N_METH_PCS:", MQTL_N_METH_PCS, "\n")

# --- 4. Save and plot --------------------------------------------------------

write.csv(tuning, file.path(MQTL_DIR, "tune_meth_pcs.csv"), row.names = FALSE)

plot_df <- tuning %>%
  select(n_meth_pcs, `FDR < 5%` = cpgs_fdr05, `p < 1e-8` = cpgs_p1e8) %>%
  pivot_longer(-n_meth_pcs, names_to = "threshold", values_to = "cpgs")

png(file.path(MQTL_DIR, "mqtl_tune_meth_pcs.png"),
    width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
print(ggplot(plot_df, aes(n_meth_pcs, cpgs, colour = threshold)) +
  geom_vline(xintercept = MQTL_N_METH_PCS, linetype = "dashed", colour = "grey50") +
  geom_line() +
  geom_point(size = 2) +
  scale_x_continuous(breaks = grid) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = paste0("CpGs with a cis-mQTL by number of latent PCs - ",
                      DATA_SOURCE, " ", paste(MQTL_ANCESTRY, collapse = "+")),
       subtitle = "Dashed line: current MQTL_N_METH_PCS",
       x = "Latent methylation PCs", y = "CpGs with a cis-mQTL", colour = NULL) +
  theme_bw())
dev.off()

cat("\nResults saved to:", file.path(MQTL_DIR, "tune_meth_pcs.csv"), "\n")
cat("Plot saved to:", file.path(MQTL_DIR, "mqtl_tune_meth_pcs.png"), "\n")
cat("Set MQTL_N_METH_PCS in mqtl/config.R, then rerun 3_mqtl.R\n")
