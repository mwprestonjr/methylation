# =============================================================================
# EWAS plot functions
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 16, 2026
# Description: Plotting helpers for ewas/4_plots.R: Manhattan, QQ, and effect /
#              p-value distribution plots of one model's per-CpG results
#              (<prefix>_results.csv.gz), drawn in base R (fast for ~900k
#              points) to DIR_FIGURES_EWAS
# =============================================================================

# Largest p-value with FDR < EWAS_FDR (NA if none), for the significance line
fdr_p_cutoff <- function(df) {
  sig <- df$P.Value[df$FDR < EWAS_FDR]
  if (length(sig)) max(sig) else NA_real_
}

# -log10(p) along the genome (hg38, CPG_ANNOTATION chr/pos_hg38; X and Y after
# the autosomes). Dashed line: EWAS_PLOT_SUGGESTIVE_P; solid line: FDR <
# EWAS_FDR. CpGs without a position (not annotated, or failed liftover) are
# left out
plot_ewas_manhattan <- function(df, main, file) {
  chr <- toupper(sub("^chr", "", as.character(df$chr), ignore.case = TRUE))
  chr_num <- suppressWarnings(as.numeric(chr))
  chr_num[chr == "X"] <- 23
  chr_num[chr == "Y"] <- 24
  keep <- !is.na(chr_num) & !is.na(df$pos_hg38)
  if (!any(keep)) {
    cat("No CpG positions; skipping the Manhattan plot\n")
    return(invisible(FALSE))
  }

  d <- data.frame(chr = chr_num[keep], pos = as.numeric(df$pos_hg38[keep]),
                  logp = -log10(df$P.Value[keep]))
  d <- d[order(d$chr, d$pos), ]

  # Cumulative position: each chromosome starts where the previous one ends
  chrs   <- sort(unique(d$chr))
  ends   <- tapply(d$pos, d$chr, max)[as.character(chrs)]
  offset <- setNames(c(0, cumsum(as.numeric(ends))[-length(chrs)]), chrs)
  d$pos_cum <- d$pos + offset[as.character(d$chr)]
  ticks  <- tapply(d$pos_cum, d$chr, median)[as.character(chrs)]
  labels <- ifelse(chrs == 23, "X", ifelse(chrs == 24, "Y", chrs))
  cols   <- c("gray30", "steelblue4")[(match(d$chr, chrs) %% 2) + 1]

  png(file, width = 10, height = 5, units = "in", res = FIG_RES)
  par(mar = c(5, 5, 4, 2))
  plot(d$pos_cum, d$logp, pch = 16, cex = 0.4, col = cols, xaxt = "n",
       xlab = "Chromosome (hg38)", ylab = expression(-log[10](italic(p))),
       main = main)
  axis(1, at = ticks, labels = labels, cex.axis = 0.7)
  abline(h = -log10(EWAS_PLOT_SUGGESTIVE_P), col = "firebrick", lty = 2)
  p_fdr <- fdr_p_cutoff(df)
  if (!is.na(p_fdr)) abline(h = -log10(p_fdr), col = "firebrick", lty = 1)
  legend("topright", bty = "n", cex = 0.8, col = "firebrick", lty = c(2, 1),
         legend = c(paste("p =", format(EWAS_PLOT_SUGGESTIVE_P)),
                    paste("FDR =", EWAS_FDR)))
  dev.off()
  invisible(TRUE)
}

# Observed vs expected -log10(p) under the null, with the genomic inflation
# factor (lambda) in the title
plot_ewas_qq <- function(df, main, file) {
  p <- df$P.Value[is.finite(df$P.Value) & df$P.Value > 0 & df$P.Value <= 1]
  if (length(p) == 0) return(invisible(FALSE))

  obs <- -log10(sort(p))
  exp <- -log10(ppoints(length(p)))

  png(file, width = 5, height = 5, units = "in", res = FIG_RES)
  par(mar = c(5, 5, 4, 2))
  plot(exp, obs, pch = 16, cex = 0.5, col = "navy",
       xlab = expression(Expected ~ -log[10](italic(p))),
       ylab = expression(Observed ~ -log[10](italic(p))),
       main = paste0(main, ", lambda = ", round(genomic_lambda(p), 3)), cex.main = 0.9)
  abline(0, 1, col = "firebrick", lty = 2)
  dev.off()
  invisible(TRUE)
}

# Histograms of the effects (logFC, M value scale) and of -log10(p)
plot_ewas_distributions <- function(df, main, file) {
  png(file, width = 10, height = 5, units = "in", res = FIG_RES)
  # Outer top margin holds `main`, shared by both panels
  par(mfrow = c(1, 2), mar = c(5, 5, 3, 2), oma = c(0, 0, 2.5, 0))
  hist(df$logFC, breaks = 100, col = "gray70", border = "white",
       main = "Effects", xlab = "Effect (M value scale)")
  abline(v = 0, lty = 2)
  # Solid line: FDR < EWAS_FDR (as on the Manhattan plot), with the CpG count
  hist(-log10(df$P.Value), breaks = 100, col = "steelblue", border = "white",
       main = "P-values", xlab = expression(-log[10](italic(p))))
  cutoff <- fdr_p_cutoff(df)
  if (!is.na(cutoff)) abline(v = -log10(cutoff), col = "firebrick")
  legend("topright", bty = "n",
         legend = sprintf("%s of %s CpGs\nFDR < %s",
                          format(sum(df$FDR < EWAS_FDR, na.rm = TRUE), big.mark = ","),
                          format(nrow(df), big.mark = ","), EWAS_FDR))
  mtext(main, side = 3, outer = TRUE, line = 0.5, cex = 1.3, font = 2)
  dev.off()
  invisible(TRUE)
}
