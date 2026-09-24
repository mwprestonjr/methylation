# =============================================================================
# Batch Correction and Cell Type Deconvolution
# Author: GP2 Subtypes and Mechanisms - M.E., M.P.
# Date: Sept 23, 2026
# Description: Loads QC-passed methylation data, performs PCA exploration,
#              confounder analysis, cell type deconvolution, and batch
#              correction via ComBat. Saves corrected M values for use
#              in differential methylation analysis.
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

library(minfi)
library(tidyverse)
library(sva)                      # ComBat
library(FlowSorted.Blood.EPIC)    # cell type deconvolution
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

# Load shared configuration
source("~/methylation/scripts/00_config.R")

# --- 1. Load QC-passed data --------------------------------------------------

cat("Loading QC-passed data from Script 02...\n")

# Load the GenomicRatioSet (normalized, filtered)
mSetSq <- readRDS(MSET_QC)
cat("GenomicRatioSet dimensions:", dim(mSetSq), "\n")

# Load M values
mVals <- readRDS(file.path(DIR_RESULTS, "mVals.rds"))
cat("M values dimensions:", dim(mVals), "\n")

# Load QC-passed sample sheet
targets_clean <- read.csv(SAMPLE_SHEET_QC,
                          colClasses = c(GP2ID            = "character",
                                         Sentrix_ID       = "character",
                                         Sentrix_Position = "character",
                                         Basename         = "character"))
cat("Sample sheet rows:", nrow(targets_clean), "\n")

# Sample order consistency check 
# colnames of mVals should match the sample order in targets_clean
cat("\n--- Sample order consistency check ---\n")
cat("mVals columns match Basename order:", 
    all(colnames(mVals) == basename(targets_clean$Basename)), "\n")

# If not aligned, reorder targets_clean to match mVals columns
if (!all(colnames(mVals) == basename(targets_clean$Basename))) {
  cat("Reordering sample sheet to match M values column order...\n")
  targets_clean <- targets_clean[match(colnames(mVals), 
                                       basename(targets_clean$Basename)), ]
}

# Use the array chip (Sentrix_ID) as the batch variable
targets_clean <- targets_clean %>%
  mutate(Batch = Sentrix_ID)

cat("\nBatch breakdown:\n")
print(table(targets_clean$Batch))
cat("Number of unique batches:", length(unique(targets_clean$Batch)), "\n")

# Flag small batches for downstream awareness
small_batches <- table(targets_clean$Batch) %>%
  as.data.frame() %>%
  filter(Freq < 5)
if (nrow(small_batches) > 0) {
  cat("\nWARNING: Small batches (<5 samples):\n")
  print(small_batches)
}

# --- 2. PCA exploration BEFORE batch correction ------------------------------

cat("\n--- PCA exploration (before correction) ---\n")

# Use top 1000 most variable probes for PCA (standard EWAS practice)
n_top_probes <- 1000
probe_vars <- rowVars(mVals)
top_var_probes <- order(probe_vars, decreasing = TRUE)[1:n_top_probes]

# PCA on transposed matrix (samples as rows for prcomp)
pca <- prcomp(t(mVals[top_var_probes, ]), scale. = TRUE, center = TRUE)

# Calculate variance explained
var_explained <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
cat("Variance explained by PC1:", round(var_explained[1], 2), "%\n")
cat("Variance explained by PC2:", round(var_explained[2], 2), "%\n")

# Build a dataframe for plotting
pca_df <- data.frame(
  PC1       = pca$x[, 1],
  PC2       = pca$x[, 2],
  Batch         = targets_clean$Batch,
  GP2_phenotype = targets_clean$GP2_phenotype,
  sex           = targets_clean$sex,
  age           = targets_clean$age
)

# 2a. PCA colored by Batch
png(file.path(DIR_RESULTS, "qc_06_pca_before_by_batch.png"),
    width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    title = "PCA Before Batch Correction - Colored by Batch",
    x     = paste0("PC1 (", round(var_explained[1], 2), "%)"),
    y     = paste0("PC2 (", round(var_explained[2], 2), "%)")
  ) +
  theme_bw() +
  theme(legend.position = "right",
        legend.text = element_text(size = 7))
dev.off()

# 2b. PCA colored by Phenotype
png(file.path(DIR_RESULTS, "qc_07_pca_before_by_phenotype.png"),
    width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
ggplot(pca_df, aes(x = PC1, y = PC2, color = GP2_phenotype)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    title = "PCA Before Batch Correction - Colored by Phenotype",
    color = "Phenotype",
    x     = paste0("PC1 (", round(var_explained[1], 2), "%)"),
    y     = paste0("PC2 (", round(var_explained[2], 2), "%)")
  ) +
  theme_bw()
dev.off()

# 2c. PCA colored by Sex
png(file.path(DIR_RESULTS, "qc_08_pca_before_by_sex.png"),
    width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
ggplot(pca_df, aes(x = PC1, y = PC2, color = sex)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c("Female" = "hotpink", "Male" = "steelblue")) +
  labs(
    title = "PCA Before Batch Correction - Colored by Sex",
    color = "Sex",
    x     = paste0("PC1 (", round(var_explained[1], 2), "%)"),
    y     = paste0("PC2 (", round(var_explained[2], 2), "%)")
  ) +
  theme_bw()
dev.off()

# 2d. PCA colored by Age (continuous)
png(file.path(DIR_RESULTS, "qc_09_pca_before_by_age.png"),
    width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
ggplot(pca_df, aes(x = PC1, y = PC2, color = age)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_viridis_c() +
  labs(
    title = "PCA Before Batch Correction - Colored by Age",
    x     = paste0("PC1 (", round(var_explained[1], 2), "%)"),
    y     = paste0("PC2 (", round(var_explained[2], 2), "%)"),
    color = "Age at\nsample collection"
  ) +
  theme_bw()
dev.off()

cat("PCA plots saved\n")

# --- 3. Statistical confounder analysis --------------------------------------

cat("\n--- Statistical confounder analysis ---\n")

# 3a. Age by phenotype (one-way ANOVA; handles more than two phenotype groups)
age_anova <- anova(lm(age ~ GP2_phenotype, data = targets_clean))
age_test  <- list(statistic = age_anova$`F value`[1], p.value = age_anova$`Pr(>F)`[1])
cat("\nAge by Phenotype (ANOVA):\n")
cat("  Mean age per phenotype:\n")
print(round(tapply(targets_clean$age, targets_clean$GP2_phenotype, mean, na.rm = TRUE), 2))
cat("  F-statistic:", round(age_test$statistic, 3), "\n")
cat("  p-value:", format(age_test$p.value, scientific = TRUE, digits = 3), "\n")

# 3b. Sex by phenotype (chi-square)
sex_table <- table(targets_clean$GP2_phenotype, targets_clean$sex)
sex_test  <- chisq.test(sex_table)
cat("\nSex by Phenotype (chi-square):\n")
print(sex_table)
cat("  chi-square:", round(sex_test$statistic, 3), "\n")
cat("  p-value:", format(sex_test$p.value, scientific = TRUE, digits = 3), "\n")

# 3c. Batch by phenotype (chi-square)
batch_table <- table(targets_clean$GP2_phenotype, targets_clean$Batch)
batch_test  <- chisq.test(batch_table)
cat("\nBatch by Phenotype (chi-square):\n")
cat("  chi-square:", round(batch_test$statistic, 3), "\n")
cat("  p-value:", format(batch_test$p.value, scientific = TRUE, digits = 3), "\n")

# Save confounder analysis summary
confounder_summary <- data.frame(
  variable    = c("Age (continuous)", "Sex (categorical)", "Batch (categorical)"),
  test        = c("ANOVA", "chi-square", "chi-square"),
  statistic   = c(age_test$statistic, sex_test$statistic, batch_test$statistic),
  p_value     = c(age_test$p.value, sex_test$p.value, batch_test$p.value),
  significant = c(age_test$p.value, sex_test$p.value, batch_test$p.value) < 0.05
)
write.csv(confounder_summary,
          file.path(DIR_RESULTS, "confounder_summary.csv"),
          row.names = FALSE)
cat("\nConfounder analysis saved to confounder_summary.csv\n")

# --- 4. Cell type deconvolution ----------------------------------------------

cat("\n--- Cell type deconvolution ---\n")
cat("Using FlowSorted.Blood.EPIC (Salas et al. 2018 IDOL probes)\n")

# We need to load the original RGChannelSet for cell type estimation
# noob normalization requires the raw RGSet, not the normalized GenomicRatioSet
# We will need to reload the rgSet from idat files
# This is a known limitation - cell type deconvolution requires raw intensities

# NOTE: For this step we need to read idat files again
# Alternative: save rgSet from Script 02 (would require modification)
# For now, we'll re-read idat files

cat("Re-reading idat files for cell type estimation...\n")
rgSet <- read.metharray.exp(targets  = targets_clean,
                            verbose  = TRUE,
                            force    = TRUE,
                            extended = TRUE)

# annotate the RGChannelSet with the appropriate array and annotation version
annotation(rgSet) <- c(array = "IlluminaHumanMethylationEPICv2",
                       annotation = "20a1.hg38")

# estimateCellCounts2 only supports 450k/EPIC v1, so run its IDOL steps directly:
# noob-normalize, map EPICv2 probes to cg IDs, project onto the IDOL reference
cat("Estimating cell counts using IDOL probes...\n")
beta_noob <- getBeta(preprocessNoob(rgSet))

# EPICv2 names carry a suffix (cg00000029_TC21); average replicates per cg ID
cg_id     <- sub("_.*$", "", rownames(beta_noob))
keep      <- cg_id %in% IDOLOptimizedCpGs
beta_idol <- rowsum(beta_noob[keep, ], cg_id[keep]) /
             as.vector(rowsum(rep(1, sum(keep)), cg_id[keep]))
cat("IDOL CpGs found on array:", nrow(beta_idol), "/", length(IDOLOptimizedCpGs), "\n")

cell_counts <- projectCellType_CP(
  Y           = beta_idol,
  coefWBC     = IDOLOptimizedCpGs.compTable[rownames(beta_idol),
                                            c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")],
  lessThanOne = FALSE
)

cell_proportions <- as.data.frame(cell_counts)
cell_proportions$Basename_short <- rownames(cell_proportions)

cat("Cell type proportions estimated for", 
    nrow(cell_proportions), "samples\n")
cat("Cell type means across all samples:\n")
print(round(colMeans(cell_proportions[, 1:6]), 4))

# Add cell type proportions to targets_clean
# Sample order consistency check
cat("\nVerifying cell counts sample order matches targets_clean...\n")
cell_proportions <- cell_proportions[
  match(basename(targets_clean$Basename), cell_proportions$Basename_short), ]

stopifnot(all(cell_proportions$Basename_short == 
              basename(targets_clean$Basename)))
cat("Sample order verified\n")

# Bind cell counts to sample sheet
targets_clean <- bind_cols(targets_clean,
                            cell_proportions[, c("CD8T", "CD4T", "NK", 
                                                 "Bcell", "Mono", "Neu")])

# Plot cell type proportions by phenotype
cell_long <- targets_clean %>%
  select(GP2ID, GP2_phenotype, CD8T, CD4T, NK, Bcell, Mono, Neu) %>%
  pivot_longer(cols = c(CD8T, CD4T, NK, Bcell, Mono, Neu),
               names_to  = "cell_type",
               values_to = "proportion")

png(file.path(DIR_RESULTS, "qc_10_cell_proportions.png"),
    width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
ggplot(cell_long, aes(x = cell_type, y = proportion, fill = GP2_phenotype)) +
  geom_boxplot(alpha = 0.7) +
  labs(
    title = "Cell Type Proportions by Phenotype",
    fill  = "Phenotype",
    x     = "Cell Type",
    y     = "Estimated Proportion"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
cat("Cell proportions plot saved\n")

# --- 5. ComBat batch correction ----------------------------------------------

cat("\n--- ComBat batch correction ---\n")
cat("Note: This is a first-pass batch correction.\n")
cat("With Project 120 data added (~500 more samples),\n")
cat("we will revisit this with a more balanced design.\n")

# Sample order consistency check before ComBat
cat("\nFinal sample order verification before ComBat:\n")
cat("colnames(mVals) match Basename in targets_clean:", 
    all(colnames(mVals) == basename(targets_clean$Basename)), "\n")

# Build model matrix protecting biological variable (GP2_phenotype)
# This tells ComBat to preserve phenotype-related variation
mod <- model.matrix(~ GP2_phenotype, data = targets_clean)

# Apply ComBat to M values
# Following Gonzalez-Latapi et al. 2023 PPMI methodology
cat("Running ComBat...\n")
combat_mVals <- ComBat(
  dat         = as.matrix(mVals),
  batch       = targets_clean$Batch,
  mod         = mod,
  par.prior   = TRUE,
  prior.plots = FALSE
)

cat("ComBat complete\n")
cat("Corrected M values dimensions:", dim(combat_mVals), "\n")

# --- 6. PCA exploration AFTER batch correction -------------------------------

cat("\n--- PCA exploration (after correction) ---\n")

# Recompute PCA on corrected M values
top_var_probes_post <- order(rowVars(combat_mVals), 
                              decreasing = TRUE)[1:n_top_probes]
pca_post <- prcomp(t(combat_mVals[top_var_probes_post, ]), 
                    scale. = TRUE, center = TRUE)
var_explained_post <- (pca_post$sdev^2 / sum(pca_post$sdev^2)) * 100

cat("Variance explained by PC1 (after):", 
    round(var_explained_post[1], 2), "%\n")
cat("Variance explained by PC2 (after):", 
    round(var_explained_post[2], 2), "%\n")

pca_df_post <- data.frame(
  PC1       = pca_post$x[, 1],
  PC2       = pca_post$x[, 2],
  Batch         = targets_clean$Batch,
  GP2_phenotype = targets_clean$GP2_phenotype,
  sex           = targets_clean$sex,
  age           = targets_clean$age
)

# 6a. PCA after - by Batch (should show LESS clustering by batch)
png(file.path(DIR_RESULTS, "qc_11_pca_after_by_batch.png"),
    width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
ggplot(pca_df_post, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    title = "PCA After Batch Correction - Colored by Batch",
    subtitle = "Successful correction = reduced clustering by batch",
    x     = paste0("PC1 (", round(var_explained_post[1], 2), "%)"),
    y     = paste0("PC2 (", round(var_explained_post[2], 2), "%)")
  ) +
  theme_bw() +
  theme(legend.position = "right",
        legend.text = element_text(size = 7))
dev.off()

# 6b. PCA after - by Phenotype (should preserve biology)
png(file.path(DIR_RESULTS, "qc_12_pca_after_by_phenotype.png"),
    width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
ggplot(pca_df_post, aes(x = PC1, y = PC2, color = GP2_phenotype)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    title = "PCA After Batch Correction - Colored by Phenotype",
    color = "Phenotype",
    subtitle = "Successful correction = biological signal preserved",
    x     = paste0("PC1 (", round(var_explained_post[1], 2), "%)"),
    y     = paste0("PC2 (", round(var_explained_post[2], 2), "%)")
  ) +
  theme_bw()
dev.off()

cat("Post-correction PCA plots saved\n")

# --- 7. Save outputs ---------------------------------------------------------

cat("\n--- Saving outputs ---\n")

# Save batch-corrected M values
saveRDS(combat_mVals, file.path(DIR_RESULTS, "combat_mVals.rds"))
cat("Batch-corrected M values saved\n")

# Save updated sample sheet with cell counts and batch
write.csv(targets_clean,
          file.path(DIR_RESULTS, "sample_sheet_final.csv"),
          row.names = FALSE)
cat("Final sample sheet saved (with cell counts and batch)\n")

# Save analysis summary
analysis_summary <- data.frame(
  metric = c("Final samples",
             "Final probes (autosomal)",
             "Number of batches",
             "Smallest batch size",
             "Age by phenotype p-value",
             "Sex by phenotype p-value",
             "Batch by phenotype p-value",
             "PC1 variance explained (before)",
             "PC1 variance explained (after)"),
  value  = c(nrow(targets_clean),
             nrow(combat_mVals),
             length(unique(targets_clean$Batch)),
             min(table(targets_clean$Batch)),
             format(age_test$p.value, scientific = TRUE, digits = 3),
             format(sex_test$p.value, scientific = TRUE, digits = 3),
             format(batch_test$p.value, scientific = TRUE, digits = 3),
             round(var_explained[1], 2),
             round(var_explained_post[1], 2))
)
write.csv(analysis_summary,
          file.path(DIR_RESULTS, "batch_correction_summary.csv"),
          row.names = FALSE)

cat("\nScript 03 complete!\n")
cat("Outputs saved to:", DIR_RESULTS, "\n")
cat("\nKey files for Script 04:\n")
cat("  - combat_mVals.rds: batch-corrected M values\n")
cat("  - sample_sheet_final.csv: sample sheet with cell counts and batch\n")