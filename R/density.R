# =============================================================================
# Beta density curves
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Sept 29, 2026
# Description: Per-sample beta density curves (as saved by Script 02 and
#              drawn by Script 03) and the density outlier score computed
#              from them
# =============================================================================

# Per-sample beta density curves, computed as in minfi::densityPlot: a list
# of x and y matrices with one column per sample (named by the columns of b)
density_curves <- function(b) {
  d <- apply(b, 2, function(x) density(as.vector(x), na.rm = TRUE))
  list(x = sapply(d, `[[`, "x"), y = sapply(d, `[[`, "y"))
}

# Density outlier score per sample (see DENSITY_MID_PEAK_MAX in qc/config.R):
# prominence of the largest peak between beta 0.15 and 0.75, where a normal
# curve has none. Prominence = a peak's height above the higher of the lowest
# points on either side of it. Curves are put on a common beta grid first so
# every sample is scored at the same resolution
density_outlier_scores <- function(curves, grid = seq(0, 1, by = 0.005)) {
  mid <- grid[grid > 0.15 & grid < 0.75]
  peak <- sapply(seq_len(ncol(curves$x)), function(i) {
    y <- approx(curves$x[, i], curves$y[, i], xout = mid, rule = 2)$y
    peaks <- which(diff(sign(diff(y))) == -2) + 1
    if (length(peaks) == 0) return(0)
    max(sapply(peaks, function(j) y[j] - max(min(y[1:j]), min(y[j:length(y)]))))
  })
  data.frame(Sample           = colnames(curves$x),
             density_mid_peak = peak)
}
