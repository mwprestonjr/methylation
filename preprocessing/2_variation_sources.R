# =============================================================================
# Preprocessing 2: Sources of variation
# Author: GP2 Subtypes and Mechanisms - M.E., M.P.
# Date: May 6, 2026
# Updated: Sept 30, 2026
# Description: Tests which technical and biological sample variables drive
#              the main components of variation in the M values: SVD of the
#              VARIATION_TOP_CPGS most variable CpGs, each of the first
#              VARIATION_N_PCS components tested against every variable
#              (heatmap of p-values), plus PCA plots and tests of age, sex and
#              batch against phenotype. Run on the QC output ("raw") to choose
#              the ComBat batch variable, and again on the ComBat output
#              ("combat") to check batch was removed and biology kept
# Usage:       Rscript preprocessing/2_variation_sources.R <data source> [raw|combat]
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

library(tidyverse)
library(matrixStats)              # rowVars

# Load shared configuration and M value helpers
source("config.R")
source("preprocessing/config.R")
source("R/mvalues.R")

# Which M values to analyse: second command-line argument (default "raw").
# In an interactive session, set VARIATION_INPUT before sourcing this script
if (!exists("VARIATION_INPUT")) VARIATION_INPUT <- commandArgs(trailingOnly = TRUE)[2]
if (is.na(VARIATION_INPUT)) VARIATION_INPUT <- "raw"
stopifnot(VARIATION_INPUT %in% c("raw", "combat"))
cat("Input:", VARIATION_INPUT, "M values\n")

# --- 1. Load data ------------------------------------------------------------

mvals_file <- if (VARIATION_INPUT == "raw") file.path(DIR_RESULTS, "mVals.rds") else COMBAT_MVALS
cat("Loading", mvals_file, "...\n")
mVals <- cap_infinite_m(readRDS(mvals_file))   # see R/mvalues.R
cat("M values dimensions:", dim(mVals), "\n")

targets <- read.csv(SAMPLE_SHEET_QC,
                    colClasses = c(GP2ID            = "character",
                                   clinical_id      = "character",
                                   Sentrix_ID       = "character",
                                   Sentrix_Position = "character",
                                   Basename         = "character")) %>%
  mutate(Sample = basename(Basename)) %>%
  left_join(read.csv(CELL_PROPORTIONS), by = "Sample")

# Samples in M value column order (ComBat output may have fewer samples)
targets <- targets[match(colnames(mVals), targets$Sample), ]
stopifnot(all(targets$Sample == colnames(mVals)))
cat("Samples:", nrow(targets), "\n")

# Variables to test: technical (array, processing) and biological. Variables
# with a single value in this data set are skipped
targets <- targets %>%
  mutate(Chip     = Sentrix_ID,
         Chip_row = substr(Sentrix_Position, 1, 3))   # R01..R08 position on the chip
technical  <- c("Chip", "Chip_row", "Batch", "Dataset")
biological <- c("phenotype", "sex", "age", CELL_TYPES)
variables  <- c(technical, biological)
variables  <- variables[sapply(variables, function(v) length(unique(na.omit(targets[[v]]))) > 1)]
# Drop a variable that groups samples exactly like an earlier one (e.g. for
# Psomagen, Batch and Dataset are both the delivery)
same_grouping <- function(a, b) {
  length(unique(paste(a, b))) == length(unique(a)) &&
    length(unique(a)) == length(unique(b))
}
variables <- variables[!sapply(seq_along(variables), function(i) {
  any(sapply(seq_len(i - 1), function(j)
    same_grouping(targets[[variables[i]]], targets[[variables[j]]])))
})]
cat("Variables tested:", paste(variables, collapse = ", "), "\n")

# --- 2. SVD of the most variable CpGs ----------------------------------------

cat("\nSVD of the", VARIATION_TOP_CPGS, "most variable CpGs...\n")
top  <- order(rowVars(mVals), decreasing = TRUE)[seq_len(min(VARIATION_TOP_CPGS, nrow(mVals)))]
x    <- mVals[top, ]
x    <- x - rowMeans(x)                  # centre each CpG
sv   <- svd(x, nu = 0, nv = VARIATION_N_PCS)
var_explained <- sv$d^2 / sum(sv$d^2) * 100
pcs  <- sv$v
colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))
rownames(pcs) <- colnames(mVals)
rm(x); invisible(gc())
cat("Variance explained by the first", VARIATION_N_PCS, "components (%):\n")
print(round(var_explained[seq_len(VARIATION_N_PCS)], 2))

# Association of each component with each variable: Kruskal-Wallis test for
# categorical variables, linear regression for numeric ones (as in ChAMP's
# champ.SVD)
assoc_p <- sapply(colnames(pcs), function(pc) {
  sapply(variables, function(v) {
    val <- targets[[v]]
    ok  <- !is.na(val)
    if (is.numeric(val)) {
      summary(lm(pcs[ok, pc] ~ val[ok]))$coefficients[2, 4]
    } else {
      kruskal.test(pcs[ok, pc], factor(val[ok]))$p.value
    }
  })
})
cat("\nAssociation p-values (rows: variables, columns: components):\n")
print(signif(assoc_p, 2))

# --- 3. Figures --------------------------------------------------------------

fig <- function(name) file.path(DIR_FIGURES, paste0("variation_", VARIATION_INPUT, "_", name, ".png"))
input_label <- if (VARIATION_INPUT == "raw") "before ComBat" else "after ComBat"

# 3a. SVD association heatmap (-log10 p, binned as in champ.SVD)
p_bins <- c("p >= 0.05", "p < 0.05", "p < 0.01", "p < 1e-5", "p < 1e-10")
heat_df <- as.data.frame(as.table(assoc_p)) %>%
  setNames(c("variable", "component", "p")) %>%
  mutate(bin = cut(p, c(-Inf, 1e-10, 1e-5, 0.01, 0.05, Inf), right = FALSE,
                   labels = rev(p_bins)),
         bin = factor(bin, levels = p_bins),
         variable  = factor(variable, levels = rev(variables)),
         component = factor(component,
                            levels = colnames(pcs),
                            labels = sprintf("%s\n(%.1f%%)", colnames(pcs),
                                             var_explained[seq_len(ncol(pcs))])))

png(fig("svd_heatmap"), width = FIG_WIDTH * 1.4, height = FIG_HEIGHT * 1.2,
    units = "in", res = FIG_RES)
print(ggplot(heat_df, aes(component, variable, fill = bin)) +
  geom_tile(colour = "white") +
  scale_fill_manual(values = c("white", "#FDD0A2", "#FD8D3C", "#D94801", "#7F2704"),
                    drop = FALSE) +
  labs(title = paste("Sources of variation -", input_label),
       x = "Component (variance explained)", y = NULL, fill = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank()))
dev.off()

# 3b. Scree plot
png(fig("scree"), width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
print(ggplot(data.frame(pc = seq_len(20), var = var_explained[seq_len(20)]),
             aes(pc, var)) +
  geom_col(fill = "grey40") +
  labs(title = paste("Variance explained per component -", input_label),
       x = "Component", y = "Variance explained (%)") +
  theme_bw())
dev.off()

# 3c. PC1 vs PC2, coloured by the main technical and biological variables
pca_df <- bind_cols(as.data.frame(pcs[, 1:2]), targets)
axis_lab <- function(i) sprintf("PC%d (%.1f%%)", i, var_explained[i])
for (v in intersect(c("Chip", "Batch", "phenotype", "sex", "age"), variables)) {
  p <- ggplot(pca_df, aes(PC1, PC2, colour = .data[[v]])) +
    geom_point(size = 1.5, alpha = 0.7) +
    labs(title = paste0("PCA ", input_label, " - coloured by ", v),
         x = axis_lab(1), y = axis_lab(2), colour = v) +
    theme_bw()
  # Many-level variables (e.g. ~50 chips) get no legend
  if (!is.numeric(pca_df[[v]]) && length(unique(pca_df[[v]])) > 12) {
    p <- p + theme(legend.position = "none")
  }
  png(fig(paste0("pca_by_", v)), width = FIG_WIDTH, height = FIG_HEIGHT,
      units = "in", res = FIG_RES)
  print(p)
  dev.off()
}
cat("Figures saved to:", DIR_FIGURES, "\n")

# --- 4. Tests of age, sex and batch against phenotype ------------------------

# A batch variable that differs by phenotype can't be corrected without
# affecting the phenotype signal; these don't depend on the M values
cat("\n--- Confounder tests (vs phenotype) ---\n")
confounder_test <- function(v) {
  d <- targets[!is.na(targets[[v]]), ]
  if (is.numeric(d[[v]])) {
    a <- anova(lm(d[[v]] ~ d$phenotype))
    data.frame(variable = v, test = "ANOVA", statistic = a$`F value`[1], p_value = a$`Pr(>F)`[1])
  } else {
    t <- suppressWarnings(chisq.test(table(d$phenotype, d[[v]])))
    data.frame(variable = v, test = "chi-square", statistic = unname(t$statistic), p_value = t$p.value)
  }
}
confounders <- bind_rows(lapply(intersect(c("age", "sex", "Batch", "Chip"), variables),
                                confounder_test)) %>%
  mutate(significant = p_value < 0.05)
print(confounders)

# --- 5. Save results ---------------------------------------------------------

# Small results behind the figures, so they can be redrawn without the M values
saveRDS(list(input         = VARIATION_INPUT,
             pcs           = pcs,
             var_explained = var_explained,
             assoc_p       = assoc_p,
             variables     = variables,
             confounders   = confounders),
        file.path(DIR_RESULTS, paste0("variation_sources_", VARIATION_INPUT, ".rds")))
write.csv(data.frame(variable = rownames(assoc_p), assoc_p, row.names = NULL),
          file.path(DIR_RESULTS, paste0("variation_svd_pvalues_", VARIATION_INPUT, ".csv")),
          row.names = FALSE)
write.csv(confounders, file.path(DIR_RESULTS, "confounder_summary.csv"), row.names = FALSE)
cat("\nResults saved to:", DIR_RESULTS, "\n")
if (VARIATION_INPUT == "raw") {
  cat("Next: set COMBAT_BATCH_VAR in preprocessing/config.R, then Rscript preprocessing/3_combat.R",
      DATA_SOURCE, "\n")
}
