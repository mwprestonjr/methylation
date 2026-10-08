# =============================================================================
# mQTL Step 4: plots
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 8, 2026
# Description: Figures and a top-hits table from the step 3 results, saved in
#              MQTL_DIR:
#   mqtl_02_top_genotype_boxplots  methylation (beta) by genotype for the
#                                  MQTL_PLOT_N_TOP strongest mQTLs at distinct loci
#   mqtl_03_regional               -log10 p of every SNP within
#                                  +/- MQTL_PLOT_REGION of the MQTL_PLOT_N_REGIONAL
#                                  top CpGs, coloured by LD (r2) with the lead SNP
#   mqtl_04_genome                 lead p-value of each CpG across the genome
#   mqtl_05_distance_zoom          lead SNP - CpG distance within +/- 50 kb
#   mqtl_06_effect_vs_maf          effect size of the lead SNP by its allele frequency
#   mqtl_07_island_context         share of tested CpGs with a cis-mQTL by CpG
#                                  island context
#   mqtl_top_hits.csv              the top mQTLs with gene, island context and MAF
#              Gene and island context come from CPG_ANNOTATION
#              (preprocessing/4_cpg_annotation.R)
#              Step 3 saves only pairs with p < MQTL_P_CIS_SAVE, so the regional
#              plots recompute every SNP in the window with the step 3 model
#              (same transformed methylation and covariates)
# Usage:       Rscript mqtl/4_mqtl_plots.R <data source> <ancestry>
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

library(tidyverse)
library(data.table)

# Load shared configuration, array profiles and mQTL functions
source("config.R")
source("mqtl/config.R")
source("R/mqtl_setup.R")
source("R/array_profiles.R")
source("R/cpg_annotation.R")
source("R/mqtl_functions.R")

fig      <- function(name) file.path(MQTL_DIR, paste0("mqtl_", name, ".png"))
geno_dir <- file.path(MQTL_DIR, "genotypes")
fdr_label    <- paste0("FDR < ", MQTL_FDR)
strict_label <- paste0("p < ", MQTL_P_STRICT)
title_suffix <- paste0(" - ", DATA_SOURCE, " ", paste(MQTL_ANCESTRY, collapse = "+"))

# --- 1. Load results ---------------------------------------------------------

cat("Loading step 3 results...\n")
lead <- fread(file.path(MQTL_DIR, "cis_mqtl_lead_per_cpg.tsv"))
cis  <- fread(file.path(MQTL_DIR, "cis_mqtl_all_saved.tsv.gz"))
covs <- read.csv(file.path(MQTL_DIR, "mqtl_covariates.csv"))       # final samples, in model order
sample_map <- read.csv(MQTL_SAMPLE_MAP, colClasses = "character")
gp2id <- sample_map$GP2ID[match(covs$meth_id, sample_map$meth_id)]
stopifnot(!anyNA(gp2id))

# p-value the FDR threshold corresponds to (NA with a warning if beyond the
# saved range); used for threshold lines
fdr_p <- fdr_threshold(cis)
rm(cis); invisible(gc())

# p-values of 0 (below machine precision) are set to the smallest positive
# double so -log10 p is finite; ties are broken by |t|
lead[, p_plot := pmax(p, .Machine$double.xmin)]
lead <- lead[order(p, -abs(t_stat))]
sig <- lead[fdr < MQTL_FDR]
cat("CpGs with a cis-mQTL at", fdr_label, ":", nrow(sig), "\n")

# CpG annotation: island context and GENCODE gene, from
# preprocessing/4_cpg_annotation.R. Label: the gene whose promoter or body
# contains the CpG, or "near <gene>" for intergenic CpGs
if (!file.exists(CPG_ANNOTATION)) {
  stop(CPG_ANNOTATION, " not found - run preprocessing/4_cpg_annotation.R ", DATA_SOURCE)
}
ann <- fread(CPG_ANNOTATION, select = c("cpg", "island", "gene", "gene_type",
                                        "gene_region", "distance_to_gene"))
ann[, island := factor(island, levels = c("Island", "Shore", "Shelf", "OpenSea"))]
ann[, gene_label := ifelse(is.na(gene), NA_character_,
                           ifelse(gene_region == "intergenic", paste("near", gene), gene))]

# --- 2. Top mQTLs at distinct loci -------------------------------------------

# Strongest first; skip a CpG within MQTL_CIS_WINDOW of one already chosen, so
# each hit is a different locus
pick_loci <- function(d, n) {
  chosen <- integer(0)
  for (i in seq_len(nrow(d))) {
    near <- d$chr[chosen] == d$chr[i] & abs(d$cpg_pos[chosen] - d$cpg_pos[i]) < MQTL_CIS_WINDOW
    if (!any(near)) chosen <- c(chosen, i)
    if (length(chosen) == n) break
  }
  d[chosen]
}
top <- pick_loci(sig, max(MQTL_PLOT_N_TOP, MQTL_PLOT_N_REGIONAL))
top <- merge(top, ann[, .(cpg, island, gene, gene_type, gene_region, distance_to_gene, gene_label)],
             by = "cpg", all.x = TRUE, sort = FALSE)
top <- top[order(p, -abs(t_stat))]
cat("Top loci:", nrow(top), "\n")

# --- 3. Methylation and genotypes --------------------------------------------

cat("Loading methylation for the top CpGs...\n")
combat     <- readRDS(COMBAT_MVALS)
tested     <- rownames(combat)                 # CpGs tested (EPICv1: a few hundred fewer after liftover)
m_top      <- combat[top$cpg, covs$meth_id, drop = FALSE]
rm(combat); invisible(gc())
beta_top   <- 2^m_top / (1 + 2^m_top)          # M -> beta, for plotting
# The step 3 model uses the inverse-normal transformed M values
int_transform <- function(x) qnorm((rank(x, ties.method = "average") - 0.5) / length(x))
int_top    <- t(apply(m_top, 1, int_transform))

# Genotypes: dosages for the lead SNPs of all significant CpGs (for MAF), and
# every SNP in the regional windows. .traw sample columns are "<FID>_<IID>"
cat("Reading genotypes...\n")
regional  <- head(top, MQTL_PLOT_N_REGIONAL)
lead_maf  <- list()
lead_geno <- list()
win_geno  <- list()
for (chr in 1:22) {
  traw <- fread(file.path(geno_dir, paste0("geno_chr", chr, ".traw")))
  ids  <- names(traw)[-(1:6)]
  col  <- sapply(gp2id, function(id) which(endsWith(ids, paste0("_", id))))
  stopifnot(is.numeric(col), length(col) == length(gp2id))
  g    <- as.matrix(traw[, ids[col], with = FALSE])
  rownames(g) <- traw$SNP

  # MAF of lead SNPs in these samples
  rows <- which(traw$SNP %in% sig$snp)
  if (length(rows)) {
    f <- rowMeans(g[rows, , drop = FALSE], na.rm = TRUE) / 2
    lead_maf[[chr]] <- data.table(snp = traw$SNP[rows], maf = pmin(f, 1 - f))
  }
  # Lead SNP genotypes for the top hits
  for (s in intersect(top$snp, traw$SNP)) lead_geno[[s]] <- g[s, ]
  # All SNPs within the window of each regional CpG (widened if needed so it
  # includes the lead SNP, which can be up to MQTL_CIS_WINDOW away)
  for (i in which(regional$chr == chr)) {
    half <- max(MQTL_PLOT_REGION, abs(regional$distance[i]) + 10e3)
    w <- which(abs(traw$POS - regional$cpg_pos[i]) <= half)
    win_geno[[regional$cpg[i]]] <- list(g = g[w, , drop = FALSE], pos = traw$POS[w])
  }
  rm(traw, g); invisible(gc())
}
lead_maf <- rbindlist(lead_maf)
top <- merge(top, lead_maf, by = "snp", all.x = TRUE, sort = FALSE)
top <- top[order(p, -abs(t_stat))]

# Genotype labels: copies of the counted allele -> letters, e.g. counted C,
# other allele G: 0 = GG, 1 = CG, 2 = CC
geno_label <- function(copies, counted, alt) {
  c(paste0(alt, alt), paste0(counted, alt), paste0(counted, counted))[copies + 1]
}

# --- 4. Genotype boxplots of the top hits ------------------------------------

box_df <- rbindlist(lapply(seq_len(min(MQTL_PLOT_N_TOP, nrow(top))), function(i) {
  h <- top[i]
  data.table(panel    = sprintf("%s%s\n%s  p = %.1e",
                                h$cpg, ifelse(is.na(h$gene_label), "", paste0(" (", h$gene_label, ")")),
                                h$snp, h$p_plot),
             order    = i,
             copies   = round(lead_geno[[h$snp]]),
             genotype = geno_label(round(lead_geno[[h$snp]]), h$counted, h$alt),
             beta     = beta_top[h$cpg, ])
}))
box_df <- box_df[!is.na(copies)]
box_df[, panel := factor(panel, levels = unique(panel[order(order)]))]
# x positions are panel-specific ("<panel>|<genotype>", ordered by copies of
# the counted allele) so each panel keeps its own genotype order; the axis
# shows only the genotype
box_df[, x := paste(order, genotype, sep = "|")]
box_df[, x := factor(x, levels = unique(x[order(order, copies)]))]

png(fig("02_top_genotype_boxplots"), width = FIG_WIDTH * 2, height = FIG_HEIGHT * 2.2,
    units = "in", res = FIG_RES)
print(ggplot(box_df, aes(x, beta)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(width = 0.15, size = 0.6, alpha = 0.5) +
  facet_wrap(~ panel, scales = "free", ncol = 4) +
  scale_x_discrete(labels = function(x) sub("^[0-9]+\\|", "", x)) +
  labs(title = paste0("Top cis-mQTLs (distinct loci)", title_suffix),
       x = "Genotype of lead SNP", y = "Methylation (beta)") +
  theme_bw(base_size = 8) +
  theme(strip.text = element_text(size = 6)))
dev.off()

# --- 5. Regional plots -------------------------------------------------------

# Every SNP in the window, tested with the step 3 model: transformed
# methylation ~ SNP + all step 3 covariates (incl. latent PCs)
C <- as.matrix(covs[, setdiff(names(covs), "meth_id")])
# Returns NA for a SNP the model can't be fitted with (no variation in these
# samples, or identical to a covariate); MatrixEQTL skips these too
snp_test <- function(y, g) {
  ok  <- !is.na(g)
  if (length(unique(g[ok])) < 2) return(NA_real_)
  X   <- cbind(1, g[ok], C[ok, , drop = FALSE])
  fit <- lm.fit(X, y[ok])
  if (fit$rank < ncol(X)) return(NA_real_)
  df  <- sum(ok) - fit$rank
  se  <- sqrt(sum(fit$residuals^2) / df * chol2inv(qr.R(fit$qr))[2, 2])
  2 * pt(-abs(fit$coefficients[2] / se), df)
}

reg_df <- rbindlist(lapply(seq_len(nrow(regional)), function(i) {
  h  <- regional[i]
  wg <- win_geno[[h$cpg]]
  y  <- int_top[h$cpg, ]
  p  <- apply(wg$g, 1, function(g) snp_test(y, g))
  r2 <- apply(wg$g, 1, function(g) suppressWarnings(cor(g, wg$g[h$snp, ], use = "complete.obs"))^2)
  # The recomputed lead p-value should match step 3
  cat(sprintf("  %s lead SNP: step 3 p = %.3g, recomputed p = %.3g\n", h$cpg, h$p, p[h$snp]))
  if (any(is.na(p))) cat("   ", sum(is.na(p)), "SNP(s) in the window skipped (no variation or collinear)\n")
  data.table(panel = sprintf("%s%s (lead %s)", h$cpg,
                             ifelse(is.na(h$gene_label), "", paste0(" ", h$gene_label)), h$snp),
             order = i, pos = wg$pos, cpg_pos = h$cpg_pos,
             logp  = -log10(pmax(p, .Machine$double.xmin)),
             r2    = cut(r2, c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
                         labels = c("< 0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "> 0.8")),
             is_lead = rownames(wg$g) == h$snp)
}))
reg_df[, panel := factor(panel, levels = unique(panel[order(order)]))]

p_reg <- ggplot(reg_df, aes(pos / 1e6, logp)) +
  geom_vline(aes(xintercept = cpg_pos / 1e6), linetype = "dashed", colour = "grey50") +
  geom_point(aes(colour = r2), size = 1) +
  geom_point(data = reg_df[is_lead == TRUE], shape = 23, size = 2.5, fill = "purple") +
  scale_colour_manual(values = c("navy", "lightblue", "green3", "orange", "red"), drop = FALSE) +
  facet_wrap(~ panel, scales = "free", ncol = 2) +
  labs(title = paste0("Regional plots of top cis-mQTLs", title_suffix),
       subtitle = "Dashed line: CpG; diamond: lead SNP",
       x = "Position (Mb, hg38)", y = "-log10 p", colour = "LD r2 with lead") +
  theme_bw(base_size = 8)
if (!is.na(fdr_p)) p_reg <- p_reg + geom_hline(yintercept = -log10(fdr_p), linetype = "dotted", colour = "red")
png(fig("03_regional"), width = FIG_WIDTH * 2, height = FIG_HEIGHT * 1.8, units = "in", res = FIG_RES)
print(p_reg)
dev.off()

# --- 6. Genome-wide lead p-values --------------------------------------------

# One point per CpG (lead SNP p-value; CpGs whose best p is above
# MQTL_P_CIS_SAVE aren't saved and so not shown), chromosomes end to end
chr_len <- lead[, .(len = max(cpg_pos)), by = chr][order(chr)]
chr_len[, offset := cumsum(as.numeric(len)) - len]
gw <- merge(lead[, .(chr, cpg_pos, p_plot)], chr_len[, .(chr, offset)], by = "chr")
gw[, x := cpg_pos + offset]
gw[, logp := pmin(-log10(p_plot), 50)]             # cap so extreme hits don't flatten the rest
p_gw <- ggplot(gw, aes(x, logp, colour = factor(chr %% 2))) +
  geom_point(size = 0.2, alpha = 0.5) +
  geom_hline(yintercept = -log10(MQTL_P_STRICT), linetype = "dashed", colour = "grey30") +
  scale_colour_manual(values = c("grey40", "steelblue"), guide = "none") +
  scale_x_continuous(breaks = chr_len$offset + chr_len$len / 2, labels = chr_len$chr) +
  labs(title = paste0("Lead cis-mQTL p-value per CpG", title_suffix),
       subtitle = paste0("Dashed: ", strict_label, "; dotted: ", fdr_label, "; -log10 p capped at 50"),
       x = "Chromosome", y = "-log10 p (lead SNP)") +
  theme_bw(base_size = 8) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
if (!is.na(fdr_p)) p_gw <- p_gw + geom_hline(yintercept = -log10(fdr_p), linetype = "dotted", colour = "red")
png(fig("04_genome"), width = FIG_WIDTH * 2, height = FIG_HEIGHT, units = "in", res = FIG_RES)
print(p_gw)
dev.off()
rm(gw); invisible(gc())

# --- 7. Lead SNP distance, zoomed --------------------------------------------

png(fig("05_distance_zoom"), width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
print(ggplot(sig[abs(distance) <= 50e3], aes(distance / 1e3)) +
  geom_histogram(binwidth = 1, boundary = 0) +
  labs(title = "Lead SNP - CpG distance within 50 kb",
       subtitle = sprintf("%s %s: %s CpGs at %s; %.0f%% of them within 50 kb", DATA_SOURCE,
                          paste(MQTL_ANCESTRY, collapse = "+"), format(nrow(sig), big.mark = ","),
                          fdr_label, 100 * mean(abs(sig$distance) <= 50e3)),
       x = "SNP - CpG distance (kb)", y = "CpGs") +
  theme_bw())
dev.off()

# --- 8. Effect size by allele frequency --------------------------------------

maf_df <- merge(sig[, .(snp, beta)], lead_maf, by = "snp")
png(fig("06_effect_vs_maf"), width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
print(ggplot(maf_df, aes(maf, abs(beta))) +
  geom_bin2d(bins = 60) +
  scale_fill_viridis_c(trans = "log10") +
  labs(title = "Effect size by allele frequency of the lead SNP",
       subtitle = paste0(DATA_SOURCE, " ", paste(MQTL_ANCESTRY, collapse = "+"), ": CpGs at ", fdr_label),
       x = "Minor allele frequency", y = "|effect| (SD per allele)",
       fill = "CpGs") +
  theme_bw())
dev.off()

# --- 9. mQTLs by CpG island context ------------------------------------------

isl <- data.table(cpg = tested)[, has_mqtl := cpg %in% sig$cpg]
isl <- merge(isl, ann[, .(cpg, island)], by = "cpg")
isl_sum <- isl[!is.na(island), .(tested = .N, with_mqtl = sum(has_mqtl)), by = island]
isl_sum[, pct := 100 * with_mqtl / tested]
print(isl_sum[order(island)])
png(fig("07_island_context"), width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
print(ggplot(isl_sum, aes(island, pct)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = sprintf("%.1f%%\n(%s)", pct, format(with_mqtl, big.mark = ",", trim = TRUE))),
            vjust = -0.2, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(title = "CpGs with a cis-mQTL by island context",
       subtitle = paste0(DATA_SOURCE, " ", paste(MQTL_ANCESTRY, collapse = "+"),
                         ": share of tested CpGs with a cis-mQTL at ", fdr_label),
       x = NULL, y = "% of tested CpGs") +
  theme_bw())
dev.off()

# --- 10. Top hits table ------------------------------------------------------

top_out <- top[, .(cpg, gene, gene_type, gene_region, distance_to_gene, island,
                   chr, cpg_pos, snp, counted, alt, maf,
                   distance, beta, t_stat, p, fdr)]
fwrite(top_out, file.path(MQTL_DIR, "mqtl_top_hits.csv"))
cat("\nFigures and mqtl_top_hits.csv saved to:", MQTL_DIR, "\n")
