# =============================================================================
# mQTL Step 4: permutations
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 27, 2026
# Description: CpG-level significance of the step 3 cis-mQTLs by
#              permutation, as FastQTL and tensorQTL (GTEx). Step 3's FDR over
#              SNP-CpG pairs treats the SNPs near a CpG as independent tests,
#              although LD correlates them. Here each CpG gets one p-value for
#              its lead SNP: the samples are shuffled MQTL_N_PERM times (the
#              same shuffles for every CpG), the best SNP of each shuffle
#              gives the CpG's null distribution, and a Beta distribution
#              fitted to it gives p_beta (cis_perm_chr, beta_approx_pval in
#              R/mqtl_functions.R). A CpG has a cis-mQTL if the q-value of
#              p_beta across CpGs is < MQTL_FDR (MQTL_PERM_FDR_METHOD:
#              Storey or Benjamini-Hochberg). The significant SNPs of each
#              such CpG are those below its nominal p-value threshold
#              (perm_significance).
#              Same model as step 3: the covariates are read from its
#              mqtl_covariates.csv, the same SNPs are tested (SNPs in the
#              probe excluded) and missing genotypes are set to the mean;
#              the lead SNPs are checked against step 3's
# Usage:       Rscript mqtl/4_mqtl_permutations.R <data source> <ancestry>
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

library(minfi)                            # getAnnotation
library(tidyverse)
library(data.table)
library(matrixStats)                      # colMaxs, colMins
library(MatrixEQTL)
library(GenomicRanges)
library(rtracklayer)                      # liftOver
library(parallel)                         # mclapply
library(RhpcBLASctl)                      # blas_set_num_threads

# Load shared configuration, array profiles and mQTL functions
source("config.R")
source("mqtl/config.R")
source("R/mqtl_setup.R")
source("R/array_profiles.R")
source("R/cpg_annotation.R")
source("R/mqtl_functions.R")

fdr_label <- paste0("q < ", MQTL_FDR, " (", MQTL_PERM_FDR_METHOD, ")")

# --- 1. Methylation data and step 3 covariates -------------------------------

prep     <- mqtl_prepare()
meth_int <- prep$meth_int
cpg_pos  <- prep$cpg_pos

# Covariates exactly as step 3 used them (incl. its latent methylation PCs),
# in the same sample order
covs_file <- file.path(MQTL_DIR, "mqtl_covariates.csv")
if (!file.exists(covs_file)) stop(covs_file, " not found - run mqtl/3_mqtl.R first")
covs <- read.csv(covs_file)
if (!identical(covs$meth_id, colnames(meth_int))) {
  stop("The samples in ", covs_file, " don't match the current inputs - rerun mqtl/3_mqtl.R")
}
C  <- as.matrix(covs[, setdiff(names(covs), "meth_id")])
qx <- qr(cbind(1, C))
stopifnot(qx$rank == ncol(C) + 1)          # MatrixEQTL stops on collinear covariates too
Q  <- qr.Q(qx)
df <- ncol(meth_int) - ncol(C) - 2         # residual df of methylation ~ SNP + covariates
cat("Samples:", ncol(meth_int), "| covariates:", ncol(C), "| residual df:", df, "\n")

# Methylation residuals after the covariates, scaled to unit length
cat("Residualizing methylation on the covariates...\n")
Y      <- residualize_rows(meth_int, Q)
y_norm <- sqrt(rowSums(Y^2))
Ys     <- Y / y_norm
# Only Ys is needed from here; prep's copy of the methylation is dropped too,
# so the forked workers start from a smaller parent
rm(Y, meth_int); prep$meth_int <- NULL; invisible(gc())

# --- 2. Permutations, per chromosome -----------------------------------------

# One set of shuffles (samples x permutations) for all CpGs, so CpGs keep
# their correlation with each other under the null, as in tensorQTL
set.seed(MQTL_PERM_SEED)
perms <- replicate(MQTL_N_PERM, sample.int(ncol(Ys)))

# Largest chromosomes first, so the parallel jobs finish at similar times
chrs <- as.integer(names(sort(table(cpg_pos$chr[cpg_pos$chr %in% as.character(1:22)]),
                              decreasing = TRUE)))
cat("\n--- Permutations:", MQTL_N_PERM, "per CpG, seed", MQTL_PERM_SEED, "|",
    MQTL_PERM_CORES, "cores ---\n")
start_time <- Sys.time()

# mc.preschedule = FALSE runs each chromosome as its own job. Each forked
# worker would otherwise run BLAS with one thread per CPU, so the workers
# would oversubscribe the CPUs many times over; one thread each instead
perm <- mclapply(chrs, function(chr) {
  blas_set_num_threads(1)
  t0     <- Sys.time()
  on_chr <- cpg_pos$chr == as.character(chr)
  geno   <- load_chr_genotypes(prep$geno_dir, chr, prep$sample_map, colnames(Ys), sliced = FALSE)
  res    <- cis_perm_chr(geno, Ys[on_chr, , drop = FALSE], y_norm[on_chr],
                         cpg_pos[on_chr, ], prep$probe_regions[on_chr, ], Q, perms, df)
  cat(sprintf("chr%d: %d CpGs, %d SNPs, %.1f min\n", chr, nrow(res), nrow(geno$geno),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  res
}, mc.cores = MQTL_PERM_CORES, mc.preschedule = FALSE)

# mclapply returns errors instead of stopping, and NULL for a worker that was
# killed (e.g. out of memory); report which chromosomes failed
failed <- vapply(perm, function(x) is.null(x) || inherits(x, "try-error"), logical(1))
if (any(failed)) {
  stop("Permutations failed for chr ", paste(chrs[failed], collapse = ", "),
       " (no result: worker killed, e.g. out of memory - lower MQTL_PERM_CORES or MQTL_PERM_BLOCK)",
       "\n", paste(unique(unlist(perm[failed])), collapse = "\n"))
}
perm <- bind_rows(perm) %>%
  left_join(cpg_pos %>% select(cpg = geneid, cpg_pos = left), by = "cpg") %>%
  mutate(distance = pos - cpg_pos)
cat("Permutation time:", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1), "min\n")
cat("CpGs tested:", nrow(perm), "\n")

# --- 3. Checks against step 3 ------------------------------------------------

# Same SNPs tested per CpG (mqtl_cpgs_tested.tsv.gz; missing for step 3
# outputs from before it was saved)
tested_file <- file.path(MQTL_DIR, "mqtl_cpgs_tested.tsv.gz")
if (file.exists(tested_file)) {
  s3 <- fread(tested_file)[n_snps > 0]
  m  <- match(s3$cpg, perm$cpg)
  if (anyNA(m) || nrow(s3) != nrow(perm) || any(s3$n_snps != perm$n_snps[m])) {
    stop("SNPs tested per CpG differ from step 3 (", tested_file, ") - rerun mqtl/3_mqtl.R")
  }
  cat("SNPs tested per CpG: same as step 3\n")
} else {
  warning(tested_file, " not found (step 3 run before it was saved): SNPs tested per CpG ",
          "not checked", call. = FALSE)
}

# Same lead p-values as step 3, for the CpGs step 3 saved (lead p <
# MQTL_P_CIS_SAVE); compared on the log scale, skipped where step 3's p is 0
# (below machine precision). The lead SNP itself can differ where SNPs tie
lead_s3 <- fread(file.path(MQTL_DIR, "cis_mqtl_lead_per_cpg.tsv"), select = c("cpg", "snp", "p", "fdr"))
cmp     <- merge(as.data.table(perm)[, .(cpg, snp, p)], lead_s3, by = "cpg", suffixes = c("", "_s3"))
ok      <- cmp$p_s3 > 0
off     <- !(abs(log10(cmp$p[ok]) / log10(cmp$p_s3[ok]) - 1) <= 1e-3)
cat("Lead p-value vs step 3:", nrow(cmp), "CpGs compared,", sum(off), "differ;",
    sum(cmp$snp != cmp$snp_s3), "with a different (tied) lead SNP\n")
if (nrow(cmp) < nrow(lead_s3) || any(off)) {
  warning(sum(off), " lead p-values differ from step 3, and ", nrow(lead_s3) - nrow(cmp),
          " step 3 CpGs are missing - check that steps 3 and 4 used the same inputs", call. = FALSE)
}

# Beta approximation against the direct permutation p-value, where the
# permutations resolve it (p_perm well above its floor 1 / (MQTL_N_PERM + 1))
resolved <- perm$p_perm > 10 / (MQTL_N_PERM + 1)
cat(sprintf("p_beta vs p_perm (p_perm > %.3g): correlation of log p %.4f, median ratio %.3f\n",
            10 / (MQTL_N_PERM + 1), cor(log(perm$p_beta[resolved]), log(perm$p_perm[resolved])),
            median(perm$p_beta[resolved] / perm$p_perm[resolved])))
cat("true_df: median", round(median(perm$true_df), 1), "(model df", df, ")\n")
# The permuted correlations are too small (a too-weak null, anti-conservative
# p-values) if true_df is well above the model's df
if (median(perm$true_df) > 1.05 * df) {
  warning("median true_df (", round(median(perm$true_df), 1), ") is well above the model df (",
          df, "): the permutation null may be too weak - check cis_perm_chr", call. = FALSE)
}

# --- 4. CpG-level significance -----------------------------------------------

sig  <- perm_significance(perm, df)
perm <- sig$perm
n_storey <- sum(perm$q_storey < MQTL_FDR)
n_bh     <- sum(perm$q_bh < MQTL_FDR)
n_sig    <- sum(perm$q < MQTL_FDR)
n_pair   <- sum(lead_s3$fdr < MQTL_FDR)
cat("\nStorey pi0:", round(sig$pi0, 3), "(estimated share of CpGs with a cis-mQTL:",
    round(1 - sig$pi0, 3), ")\n")
cat("CpGs with a cis-mQTL, q < ", MQTL_FDR, ": Storey ", n_storey, " | BH ", n_bh,
    " | used (", MQTL_PERM_FDR_METHOD, "): ", n_sig, "\n", sep = "")
cat("For comparison, step 3 (pair-level BH FDR <", MQTL_FDR, "):", n_pair, "\n")
cat("p_beta threshold:", signif(sig$p_beta_threshold, 3), "\n")

# --- 5. Significant SNP-CpG pairs --------------------------------------------

# Saved step 3 pairs below their CpG's nominal threshold. Step 3 saves only
# pairs with p < MQTL_P_CIS_SAVE, so a CpG whose threshold is above that may
# have significant pairs that weren't saved. These are CpGs with few SNPs in
# their window (a CpG with one SNP is significant at about the p_beta
# threshold); a handful is expected and only their pair lists are incomplete
sig_cpgs <- perm %>% filter(q < MQTL_FDR) %>% select(cpg, n_snps, p_nominal_threshold)
above    <- sig_cpgs$p_nominal_threshold >= MQTL_P_CIS_SAVE
n_above  <- sum(above)
if (n_above > 0) {
  msg <- sprintf(paste0("%d significant CpGs (%.2f%%; median %d SNPs in the window) have a ",
                        "nominal threshold >= MQTL_P_CIS_SAVE (%g), so their lists of ",
                        "significant pairs are incomplete"),
                 n_above, 100 * mean(above), as.integer(median(sig_cpgs$n_snps[above])),
                 MQTL_P_CIS_SAVE)
  if (mean(above) > 0.01) {
    warning(msg, ". Loosen MQTL_P_CIS_SAVE in mqtl/config.R and rerun steps 3 and 4",
            call. = FALSE)
  } else {
    cat(msg, "(expected for CpGs with few SNPs; per-CpG results are complete)\n")
  }
}
sig_cpgs <- select(sig_cpgs, -n_snps)
pairs <- fread(file.path(MQTL_DIR, "cis_mqtl_all_saved.tsv.gz"))
pairs <- merge(pairs, as.data.table(sig_cpgs), by = "cpg")[p < p_nominal_threshold]
setorder(pairs, p)
cat("Significant pairs:", nrow(pairs), "for", uniqueN(pairs$cpg), "CpGs\n")

# --- 6. Figures --------------------------------------------------------------

title_suffix <- paste0(" - ", DATA_SOURCE, " ", paste(MQTL_ANCESTRY, collapse = "+"))

# p_beta histogram: flat for CpGs without a cis-mQTL, at height pi0
png(file.path(MQTL_FIG_DIR, "mqtl_perm_01_pvalues.png"),
    width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
print(ggplot(perm, aes(p_beta)) +
  geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 1, 0.02), fill = "grey60") +
  geom_hline(yintercept = sig$pi0, linetype = "dashed", colour = "red") +
  labs(title = paste0("Permutation p-values (Beta approx.)", title_suffix),
       subtitle = sprintf("Dashed: Storey pi0 = %.3f; %s CpGs at %s", sig$pi0,
                          format(n_sig, big.mark = ","), fdr_label),
       x = "p_beta (lead SNP, per CpG)", y = "Density") +
  theme_bw())
dev.off()

# Beta approximation against the direct permutation p-value. -log10 p_beta is
# capped before binning, so the bins cover the plotted range rather than the
# strongest hits (p_beta down to ~1e-300)
y_cap <- 3 * log10(MQTL_N_PERM + 1)
png(file.path(MQTL_FIG_DIR, "mqtl_perm_02_beta_vs_empirical.png"),
    width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
print(ggplot(perm, aes(-log10(p_perm), pmin(-log10(p_beta), y_cap))) +
  geom_bin2d(bins = 80) +
  geom_abline(linetype = "dashed", colour = "red") +
  geom_vline(xintercept = log10(MQTL_N_PERM + 1), linetype = "dotted", colour = "grey40") +
  scale_fill_viridis_c(trans = "log10") +
  labs(title = paste0("Beta vs direct permutation p-value", title_suffix),
       subtitle = sprintf("Dotted: smallest p_perm (%d permutations); p_beta capped at 1e-%.0f",
                          MQTL_N_PERM, y_cap),
       x = "-log10 p_perm", y = "-log10 p_beta", fill = "CpGs") +
  theme_bw())
dev.off()

# --- 7. Save outputs ---------------------------------------------------------

cat("\n--- Saving outputs ---\n")
perm <- perm %>%
  select(cpg, chr, cpg_pos, n_snps, snp, pos, distance, counted, alt, beta, t_stat, p,
         p_perm, true_df, p_true_df, beta_shape1, beta_shape2, p_beta,
         q_storey, q_bh, q, p_nominal_threshold) %>%
  arrange(p_beta, p)
fwrite(perm,  file.path(MQTL_DIR, "cis_mqtl_perm_per_cpg.tsv.gz"), sep = "\t")
fwrite(pairs, file.path(MQTL_DIR, "cis_mqtl_perm_signif_pairs.tsv.gz"), sep = "\t")

perm_summary <- data.frame(
  metric = c("Samples", "Ancestry", "CpGs tested", "Permutations", "Seed",
             "FDR method", "Storey pi0",
             paste0("CpGs with cis-mQTL (Storey q < ", MQTL_FDR, ")"),
             paste0("CpGs with cis-mQTL (BH q < ", MQTL_FDR, ")"),
             paste0("CpGs with cis-mQTL (step 3, pair-level FDR < ", MQTL_FDR, ")"),
             "p_beta threshold", "Significant pairs",
             paste0("Significant CpGs with nominal threshold >= ", MQTL_P_CIS_SAVE)),
  value  = c(ncol(Ys), paste(MQTL_ANCESTRY, collapse = "+"), nrow(perm), MQTL_N_PERM,
             MQTL_PERM_SEED, MQTL_PERM_FDR_METHOD, signif(sig$pi0, 4),
             n_storey, n_bh, n_pair, signif(sig$p_beta_threshold, 4),
             nrow(pairs), n_above)
)
print(perm_summary)
write.csv(perm_summary, file.path(MQTL_DIR, "mqtl_perm_summary.csv"), row.names = FALSE)

cat("\nmQTL step 4 complete!\n")
cat("Results saved to:", MQTL_DIR, "| figures to:", MQTL_FIG_DIR, "\n")
