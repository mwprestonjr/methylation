# =============================================================================
# QC Plots
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Sept 25, 2026
# Description: Plots the QC results of Script 02 from its saved figure data
#              (QC_FIGURE_DATA) and saves the figures to DIR_FIGURES:
#   qc_01  Median methylated vs unmethylated intensity per sample (minfi
#          plotQC); low-intensity samples are flagged in red
#   qc_02  Mean detection p-value per sample, sorted, log scale, with the
#          DETECTION_P_THRESHOLD line
#   qc_02b Fraction of cg probes detected (SeSAMe pOOBAH) per sample, sorted,
#          coloured by dataset, with the SESAME_MIN_FRAC_DETECTED line
#          (only if Script 02 ran with USE_SESAME_QC)
#   qc_03  Sex prediction: X vs Y chromosome median intensity, fill =
#          predicted sex, border = reported sex
#   qc_04  Beta value densities before normalization, by phenotype
#   qc_04b Beta value densities before normalization, by dataset
#   qc_05  Beta value densities after normalization (Funnorm), by phenotype
#          (qc_04 to qc_05 include only samples that passed QC)
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------
# Load shared configuration
source("~/methylation/scripts/00_config.R")

cat("Loading QC figure data...\n")
fig_data <- readRDS(QC_FIGURE_DATA)
qc_metrics <- fig_data$samples
params     <- fig_data$params
cat("Samples loaded:", nrow(qc_metrics), "\n")

# --- 1. Median intensities (qc_01) -------------------------------------------

png(file.path(DIR_FIGURES, "qc_01_median_intensities.png"), width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
minfi::plotQC(qc_metrics[, c("mMed", "uMed")])
dev.off()
cat("QC plot saved\n")

# --- 2. Detection p-values (qc_02) -------------------------------------------

# Mean detection p-value per sample, sorted, on a log scale so the threshold
# (usually far above typical values) stays visible
DETECTION_P_THRESHOLD <- params$DETECTION_P_THRESHOLD
mean_detP <- qc_metrics$mean_detP
ord      <- order(mean_detP, decreasing = TRUE)
detP_ord <- mean_detP[ord]
failing  <- detP_ord > DETECTION_P_THRESHOLD
y_lim    <- c(10^floor(log10(min(detP_ord))),
              10 * max(detP_ord, DETECTION_P_THRESHOLD))
y_ticks  <- 10^seq(log10(y_lim[1]), log10(y_lim[2]))

png(file.path(DIR_FIGURES, "qc_02_detection_pvalues.png"), width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
par(mar = c(3, 5.5, 3, 1), las = 1)
plot(detP_ord,
     log  = "y",
     ylim = y_lim,
     pch  = 16,
     cex  = 0.8,
     col  = ifelse(failing, "red", "grey30"),
     xaxt = "n",
     yaxt = "n",
     xlab = "",
     ylab = "",
     main = "Mean Detection P-value per Sample")
axis(2, at = y_ticks, labels = format(y_ticks, scientific = TRUE))
title(ylab = "Mean detection p-value (log scale)", line = 4)
mtext(sprintf("Samples, sorted (n = %d)", length(detP_ord)), side = 1, line = 1)
abline(h   = DETECTION_P_THRESHOLD,
       col = "red",
       lty = 2)
legend("topright",
       legend = sprintf("Threshold (%g)", DETECTION_P_THRESHOLD),
       lty    = 2,
       col    = "red",
       bty    = "n")
dev.off()
cat("Detection p-value plot saved\n")

# --- 3. SeSAMe detection (qc_02b) --------------------------------------------

# Fraction of cg probes detected per sample, sorted, coloured by dataset.
# Only available if Script 02 ran with USE_SESAME_QC
if (params$USE_SESAME_QC) {
  SESAME_MIN_FRAC_DETECTED <- params$SESAME_MIN_FRAC_DETECTED
  ord       <- order(qc_metrics$frac_dt_cg)
  frac_ord  <- qc_metrics$frac_dt_cg[ord]
  datasets  <- sort(unique(qc_metrics$Dataset))
  ds_colors <- setNames(palette.colors(length(datasets) + 1, "Okabe-Ito")[-1], datasets)

  png(file.path(DIR_FIGURES, "qc_02b_sesame_detection.png"), width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
  par(mar = c(3, 4.5, 3, 1), las = 1)
  plot(frac_ord,
       ylim = range(c(frac_ord, SESAME_MIN_FRAC_DETECTED, 1)),
       pch  = 16,
       cex  = 0.8,
       col  = ds_colors[qc_metrics$Dataset[ord]],
       xaxt = "n",
       xlab = "",
       ylab = "Fraction of cg probes detected",
       main = "SeSAMe Detection (pOOBAH) per Sample")
  mtext(sprintf("Samples, sorted (n = %d)", length(frac_ord)), side = 1, line = 1)
  abline(h   = SESAME_MIN_FRAC_DETECTED,
         col = "red",
         lty = 2)
  legend("bottomright",
         legend = c(datasets, sprintf("Threshold (%g)", SESAME_MIN_FRAC_DETECTED)),
         col    = c(ds_colors, "red"),
         pch    = c(rep(16, length(datasets)), NA),
         lty    = c(rep(NA, length(datasets)), 2),
         bty    = "n")
  dev.off()
  cat("SeSAMe detection plot saved\n")
} else {
  cat("Skipping SeSAMe detection plot (Script 02 ran with USE_SESAME_QC = FALSE)\n")
}

# --- 4. Sex prediction (qc_03) -----------------------------------------------

# Fill = predicted sex, border = reported sex; discordant samples drawn on top
sex_colors <- c(F = "hotpink", M = "steelblue", Unknown = "grey60")
fill_col   <- sex_colors[qc_metrics$predicted_sex]
border_col <- sex_colors[qc_metrics$reported_sex_label]
discordant_on_top <- order(!is.na(qc_metrics$sex_discordant) & qc_metrics$sex_discordant)

png(file.path(DIR_FIGURES, "qc_03_sex_prediction.png"), width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
plot(qc_metrics$xMed[discordant_on_top], qc_metrics$yMed[discordant_on_top],
     pch  = 21,
     cex  = 1.4,
     bg   = fill_col[discordant_on_top],
     col  = border_col[discordant_on_top],
     lwd  = 1.5,
     xlab = "X chromosome median intensity",
     ylab = "Y chromosome median intensity",
     main = "Sex Prediction")
mtext("Fill: predicted sex   Border: reported sex", side = 3, line = 0.3, cex = 0.7)
legend("right",
       legend = c("Female", "Male", "Reported unknown", "Discordant"),
       pch    = 21,
       pt.bg  = c("hotpink", "steelblue", "hotpink", "hotpink"),
       col    = c("hotpink", "steelblue", "grey60", "steelblue"),
       pt.cex = 1.4,
       pt.lwd = 1.5,
       cex    = 0.8)
dev.off()
cat("Sex prediction plot saved\n")

# --- 5. Beta density plots (qc_04, qc_04b, qc_05) ----------------------------

# Same drawing as minfi::densityPlot, but from the saved density curves
# (x and y matrices, one column per sample) instead of the beta matrix
density_plot <- function(curves, sampGroups, main,
                         pal = RColorBrewer::brewer.pal(8, "Dark2")) {
  sampGroups <- as.factor(sampGroups)
  plot(x    = 0,
       type = "n",
       xlim = range(curves$x),
       ylim = range(curves$y),
       ylab = "Density",
       xlab = "Beta",
       main = main)
  abline(h = 0, col = "grey80")
  for (i in seq_len(ncol(curves$x))) {
    lines(curves$x[, i], curves$y[, i], col = pal[sampGroups[i]])
  }
  if (length(levels(sampGroups)) > 1) {
    legend("topright", legend = levels(sampGroups), text.col = pal)
  }
}

# Sample metadata in the column order of each set of curves (kept samples only)
groups_for <- function(curves, var) {
  qc_metrics[[var]][match(colnames(curves$x), qc_metrics$Sample)]
}

png(file.path(DIR_FIGURES, "qc_04_density_before_normalization.png"), width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
density_plot(fig_data$density_before,
             sampGroups = groups_for(fig_data$density_before, "GP2_phenotype"),
             main       = "Beta Values - Before Normalization")
dev.off()

png(file.path(DIR_FIGURES, "qc_04b_density_by_dataset.png"), width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
density_plot(fig_data$density_before,
             sampGroups = groups_for(fig_data$density_before, "Dataset"),
             main       = "Beta Values by Dataset - Before Normalization")
dev.off()

png(file.path(DIR_FIGURES, "qc_05_density_after_normalization.png"), width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
density_plot(fig_data$density_after,
             sampGroups = groups_for(fig_data$density_after, "GP2_phenotype"),
             main       = "Beta Values - After Normalization")
dev.off()
cat("Density plots saved\n")

cat("\nQC plots complete!\n")
cat("Figures saved to:", DIR_FIGURES, "\n")
