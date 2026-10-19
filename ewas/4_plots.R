# =============================================================================
# EWAS 4: plots
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 16, 2026
# Description: For every model with results in DIR_RESULTS_EWAS
#              (<prefix>_results.csv.gz from steps 1-3), draws a Manhattan
#              plot, a QQ plot and effect / p-value distributions to
#              DIR_FIGURES_EWAS, and combines the models' <prefix>_summary.csv
#              rows into EWAS_MODEL_SUMMARY (ewas/config.R)
# Usage:       Rscript ewas/4_plots.R <data source>   (see config.R)
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

library(tidyverse)

# Load shared configuration, EWAS helpers and plot helpers
source("config.R")
source("ewas/config.R")
source("R/ewas_functions.R")
source("R/ewas_plot_functions.R")

# Models run: one results file each (skipped models leave none)
prefixes <- sub("_results\\.csv\\.gz$", "",
                list.files(DIR_RESULTS_EWAS, pattern = "_results\\.csv\\.gz$"))
if (length(prefixes) == 0) {
  stop("No EWAS results in ", DIR_RESULTS_EWAS, " - run steps 1-3 first")
}
cat("Models:", prefixes, "\n")

# --- 1. Figures per model ----------------------------------------------------

for (prefix in prefixes) {
  cat("\nPlotting", prefix, "...\n")
  files <- ewas_output_files(prefix)
  df    <- fread(files$results, select = c("cpg", "logFC", "P.Value", "FDR", "chr", "pos_hg38"),
                 colClasses = c(chr = "character"))
  label <- read.csv(files$summary)$label[1]
  title <- paste0("EWAS: ", label, " (", DATA_SOURCE, ")")

  plot_ewas_manhattan(df, title, file.path(DIR_FIGURES_EWAS, paste0(prefix, "_manhattan.png")))
  plot_ewas_qq(df, paste0(label, "\n", DATA_SOURCE), file.path(DIR_FIGURES_EWAS, paste0(prefix, "_qq.png")))
  plot_ewas_distributions(df, label, file.path(DIR_FIGURES_EWAS, paste0(prefix, "_distributions.png")))
  rm(df); invisible(gc())
}
cat("\nFigures saved to:", DIR_FIGURES_EWAS, "\n")

# --- 2. Combined summary -----------------------------------------------------

combined <- bind_rows(lapply(prefixes, function(p) read.csv(ewas_output_files(p)$summary)))
print(as.data.frame(combined %>% select(model, n_samples, n_cpg, n_fdr_sig, min_p, lambda)))
write.csv(combined, EWAS_MODEL_SUMMARY, row.names = FALSE)
cat("Combined EWAS summary saved to:", EWAS_MODEL_SUMMARY, "\n")
