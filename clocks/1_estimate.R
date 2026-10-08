# =============================================================================
# Clocks 1: estimate epigenetic age
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Sept 23, 2026
# Updated: Oct 8, 2026
# Description: Loads normalized, unfiltered beta values from qc/02_qc.R,
#              collapses EPICv2 probe names to 450K/EPIC IDs, records each clock's
#              CpG coverage, estimates epigenetic age (methylclock) and pace
#              of aging (DunedinPACE), and computes age acceleration vs
#              chronological age. Saves per-sample estimates with the model
#              covariates for 2_compare.R, and the accuracy and coverage
#              tables for 3_plots.R. Needs the preprocessing outputs (cell
#              proportions, ComBat batch in SAMPLE_SHEET_FINAL)
# Usage:       Rscript clocks/1_estimate.R <data source>   (see config.R)
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------
library(tidyverse)
library(methylclock)
library(DunedinPACE)

# Load shared configuration
source("config.R")
source("clocks/config.R")

# --- 1. Load unfiltered betas and sample sheet -------------------------------

cat("Loading unfiltered beta values...\n")
bVals <- readRDS(BVALS_UNFILTERED)
cat("Beta values dimensions:", dim(bVals), "\n")

# Final sample sheet: QC-passed samples with cell proportions and the ComBat
# batch (preprocessing/3_combat.R)
targets <- read.csv(SAMPLE_SHEET_FINAL,
                    colClasses = c(GP2ID            = "character",
                                   clinical_id      = "character",
                                   Sentrix_ID       = "character",
                                   Sentrix_Position = "character",
                                   Basename         = "character",
                                   combat_batch     = "character")) %>%
  mutate(Sample   = basename(Basename),               # minfi sample names: <Sentrix_ID>_<Sentrix_Position>
         chip_row = substr(Sentrix_Position, 1, 3))   # R01-R08

# Sample sheet in beta column order
targets <- targets[match(colnames(bVals), targets$Sample), ]
stopifnot(all(colnames(bVals) == targets$Sample))

# --- 2. Collapse EPICv2 probe names to 450K/EPIC IDs -------------------------

cat("\nCollapsing EPICv2 probe names (e.g. cg00000029_TC21 -> cg00000029)...\n")

# EPICv2 has replicate probes for some CpGs; average them into one row
collapse_epicv2 <- function(b) {
  cg_id <- sub("_.*$", "", rownames(b))
  keep  <- startsWith(cg_id, "cg") | startsWith(cg_id, "ch.")   # clocks use cg and ch. probes
  b     <- b[keep, , drop = FALSE]
  cg_id <- cg_id[keep]

  dup    <- cg_id %in% cg_id[duplicated(cg_id)]
  single <- b[!dup, , drop = FALSE]
  rownames(single) <- cg_id[!dup]

  sums  <- rowsum(b[dup, , drop = FALSE], cg_id[dup], na.rm = TRUE)
  n     <- rowsum((!is.na(b[dup, , drop = FALSE])) * 1, cg_id[dup])
  multi <- sums / n

  rbind(single, multi)
}

bVals_cg <- collapse_epicv2(bVals)
rm(bVals); invisible(gc())
cat("Collapsed beta values dimensions:", dim(bVals_cg), "\n")

# --- 3. Clock CpG coverage ---------------------------------------------------

# EPICv2 dropped some CpGs used by older clocks. methylclock leaves CpGs
# missing from the array out of the weighted sum (values missing for some
# samples are imputed), which shifts clock age: a constant offset that age
# acceleration removes, but also lost information, so clocks with low coverage
# should be interpreted with caution
cat("\nChecking clock CpG coverage...\n")

# CpGs per clock (methylclockData coefficient tables; first row is the
# intercept except for Hannum) and DunedinPACE (model probes, and the probes
# PACEProjector uses to normalize; it needs PACE_COV_THRESH of both)
clock_cpgs <- list(
  Horvath     = methylclockData::get_coefHorvath()$CpGmarker[-1],
  Hannum      = methylclockData::get_coefHannum()$CpGmarker,
  Levine      = methylclockData::get_coefLevine()$CpGmarker[-1],
  skinHorvath = methylclockData::get_coefSkin()$CpGmarker[-1],
  PedBE       = methylclockData::get_coefPedBE()$CpGmarker[-1],
  Wu          = methylclockData::get_coefWu()$CpGmarker[-1],
  TL          = methylclockData::get_coefTL()$CpGmarker[-1],
  BLUP        = methylclockData::get_coefBLUP()$CpGmarker[-1],
  EN          = methylclockData::get_coefEN()$CpGmarker[-1],
  DunedinPACE = getRequiredProbes()$DunedinPACE,
  DunedinPACE_normalization = getRequiredProbes(backgroundList = TRUE)$DunedinPACE
)

# pct_on_array: share of the clock's CpGs on the array, what the coverage
# thresholds are applied to
clock_coverage <- imap_dfr(clock_cpgs, function(cpgs, clock) {
  n_on_array <- sum(cpgs %in% rownames(bVals_cg))
  tibble(clock        = clock,
         n_cpgs       = length(cpgs),
         n_on_array   = n_on_array,
         pct_on_array = 100 * n_on_array / length(cpgs))
})

# --- 4. Estimate epigenetic age ----------------------------------------------

cat("\nEstimating epigenetic age with methylclock...\n")
# normalize = FALSE: data are already Funnorm-normalized in qc/02_qc.R
# cell.count = FALSE: cell proportions come from preprocessing/1_cell_counts.R
dnam_age <- DNAmAge(bVals_cg,
                    clocks     = "all",
                    normalize  = FALSE,
                    cell.count = FALSE,
                    min.perc   = CLOCK_COV_THRESH)

cat("\nEstimating DunedinPACE...\n")
pace <- PACEProjector(bVals_cg, proportionOfProbesRequired = PACE_COV_THRESH)

# Record which clocks were estimated (below the coverage threshold a clock is
# all NA); the normalization probes row has no estimate of its own
estimated <- c(sapply(setdiff(names(clock_cpgs), c("DunedinPACE", "DunedinPACE_normalization")),
                      function(clock) clock %in% names(dnam_age) && !all(is.na(dnam_age[[clock]]))),
               DunedinPACE = !all(is.na(pace$DunedinPACE)))
clock_coverage <- clock_coverage %>%
  mutate(estimated = unname(estimated[clock]))
cat("\nClock CpG coverage:\n")
print(clock_coverage %>% mutate(pct_on_array = round(pct_on_array, 1)))

# Attach phenotype and the model covariates from the final sample sheet
covariate_cols <- setdiff(CLOCK_COVARIATES, c("age", "sex"))
clocks <- dnam_age %>%
  rename(Sample = id) %>%
  mutate(DunedinPACE = pace$DunedinPACE[Sample]) %>%
  left_join(targets %>% select(Sample, GP2ID, phenotype, age, sex, all_of(covariate_cols)),
            by = "Sample")

# --- 5. Age acceleration and clock accuracy ----------------------------------

cat("\nComputing age acceleration...\n")
cat("Samples with age available:", sum(!is.na(clocks$age)), "/", nrow(clocks), "\n")

# Age acceleration = residual of clock age regressed on chronological age
for (clock in ADULT_CLOCKS) {
  fit <- lm(clocks[[clock]] ~ clocks$age, na.action = na.exclude)
  clocks[[paste0("AgeAccel_", clock)]] <- residuals(fit)
}

# Accuracy of each clock vs chronological age
# MAE = median absolute error; mean_offset = mean(clock age - chronological age)
clock_accuracy <- map_dfr(ADULT_CLOCKS, function(clock) {
  diff <- clocks[[clock]] - clocks$age
  tibble(clock       = clock,
         r           = cor(clocks[[clock]], clocks$age, use = "complete.obs"),
         MAE         = median(abs(diff), na.rm = TRUE),
         mean_offset = mean(diff, na.rm = TRUE))
})
cat("\nClock accuracy vs chronological age:\n")
print(clock_accuracy %>% mutate(across(where(is.numeric), ~ round(.x, 2))))

write.csv(clocks, CLOCK_RESULTS, row.names = FALSE)
write.csv(clock_accuracy, CLOCK_ACCURACY, row.names = FALSE)
write.csv(clock_coverage, CLOCK_COVERAGE, row.names = FALSE)
cat("Clock estimates, age acceleration, accuracy and coverage saved\n")

cat("\nMethylation clock estimation complete!\n")
cat("Results saved to:", DIR_RESULTS, "\n")
cat("Next: Rscript clocks/2_compare.R", DATA_SOURCE, "\n")
