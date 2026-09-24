# =============================================================================
# Methylation Clocks
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Sept 23, 2026
# Description: Loads normalized, unfiltered beta values from Script 02,
#              collapses EPICv2 probe names to cg IDs, estimates epigenetic
#              age (methylclock) and pace of aging (DunedinPACE), and
#              computes age acceleration vs chronological age at collection
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------
suppressWarnings({
  library(tidyverse)
  library(methylclock)
  library(DunedinPACE)
})

# Load shared configuration
source("~/methylation/scripts/00_config.R")

# Adult clocks used for age acceleration (pediatric/telomere clocks are skipped)
ADULT_CLOCKS <- c("Horvath", "Hannum", "Levine", "skinHorvath", "EN", "BLUP")

# --- 1. Load unfiltered betas and sample sheet -------------------------------

cat("Loading unfiltered beta values...\n")
bVals <- readRDS(BVALS_UNFILTERED)
cat("Beta values dimensions:", dim(bVals), "\n")

targets_clean <- read.csv(SAMPLE_SHEET_QC,
                          colClasses = c(Sentrix_ID       = "character",
                                         Sentrix_Position = "character",
                                         Basename         = "character"))

# minfi names samples by the idat basename (<Sentrix_ID>_<Sentrix_Position>)
targets_clean <- targets_clean[match(colnames(bVals),
                                     basename(targets_clean$Basename)), ]
stopifnot(all(colnames(bVals) == basename(targets_clean$Basename)))

# --- 2. Collapse EPICv2 probe names to cg IDs --------------------------------

cat("\nCollapsing EPICv2 probe names (e.g. cg00000029_TC21 -> cg00000029)...\n")

# EPICv2 has replicate probes for some CpGs; average them into one row
collapse_epicv2 <- function(b) {
  cg_id <- sub("_.*$", "", rownames(b))
  keep  <- startsWith(cg_id, "cg")
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
rm(bVals); gc()
cat("Collapsed beta values dimensions:", dim(bVals_cg), "\n")

# --- 3. Check clock CpG coverage ---------------------------------------------

# EPICv2 dropped some CpGs used by older clocks; missing ones get imputed,
# so clocks with low coverage should be interpreted with caution
cat("\nChecking clock CpG coverage...\n")
invisible(checkClocks(bVals_cg))

# --- 4. Estimate epigenetic age ----------------------------------------------

cat("\nEstimating epigenetic age with methylclock...\n")
# normalize = FALSE: data are already Funnorm-normalized in Script 02
# cell.count = FALSE: cell proportions come from the deconvolution step
dnam_age <- DNAmAge(bVals_cg,
                    clocks     = "all",
                    normalize  = FALSE,
                    cell.count = FALSE,
                    min.perc   = 0.8)

cat("\nEstimating DunedinPACE...\n")
pace <- PACEProjector(bVals_cg, proportionOfProbesRequired = 0.8)

# Attach phenotype, sex and age from the QC-passed sample sheet (merged from R12 in Script 01)
clocks <- dnam_age %>%
  rename(Sample = id) %>%
  mutate(DunedinPACE = pace$DunedinPACE[Sample]) %>%
  left_join(targets_clean %>%
              mutate(Sample = basename(Basename)) %>%
              select(Sample, GP2ID, GP2_phenotype, sex, age),
            by = "Sample")

write.csv(clocks, file.path(DIR_RESULTS, "clock_estimates.csv"), row.names = FALSE)
cat("Clock estimates saved\n")

# --- 5. Age acceleration -----------------------------------------------------

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
write.csv(clock_accuracy, file.path(DIR_RESULTS, "clock_accuracy.csv"), row.names = FALSE)

write.csv(clocks, file.path(DIR_RESULTS, "clock_age_acceleration.csv"), row.names = FALSE)

# Clock age vs chronological age (shared square axes so distance from y = x is comparable)
clocks_long <- clocks %>%
  pivot_longer(all_of(ADULT_CLOCKS), names_to = "clock", values_to = "DNAm_age")
axis_lim <- range(c(clocks_long$age, clocks_long$DNAm_age), na.rm = TRUE)

accuracy_labels <- clock_accuracy %>%
  mutate(label = sprintf("r = %.2f\nMAE = %.1f", r, MAE))

png(file.path(DIR_RESULTS, "clock_01_vs_chronological_age.png"), width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
p <- ggplot(clocks_long, aes(age, DNAm_age)) +
  geom_abline(linetype = "dashed", colour = "grey50") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.6) +
  geom_point(aes(colour = GP2_phenotype), alpha = 0.7) +
  geom_text(data = accuracy_labels, aes(label = label),
            x = axis_lim[1], y = axis_lim[2], hjust = 0, vjust = 1, size = 3) +
  facet_wrap(~ clock) +
  coord_equal(xlim = axis_lim, ylim = axis_lim) +
  labs(x = "Chronological age", y = "Epigenetic age", colour = "Phenotype",
       caption = "Dashed: epigenetic = chronological age. Solid: linear fit. MAE: median absolute error (years).") +
  theme_bw()
print(p)
dev.off()

# --- 6. Compare age acceleration between PD and controls ---------------------

# Restrict to PD vs Control, with Control as the reference level
clocks_cc <- clocks %>%
  filter(GP2_phenotype %in% c("PD", "Control")) %>%
  mutate(GP2_phenotype = factor(GP2_phenotype, levels = c("Control", "PD")))
cat("\nSamples in PD vs Control comparison:\n")
print(table(clocks_cc$GP2_phenotype))

# Skip the comparison if either group is missing (e.g. small test runs)
if (any(table(clocks_cc$GP2_phenotype) < 2)) {
  cat("\nFewer than 2 samples in PD or Control; skipping PD vs Control comparison\n")
  cat("\nMethylation clock pipeline complete!\n")
  quit(save = "no")
}

# NOTE: add cell proportions as covariates once deconvolution is run
cat("\nAge acceleration, PD vs Control (adjusted for sex):\n")
accel_results <- map_dfr(c(paste0("AgeAccel_", ADULT_CLOCKS), "DunedinPACE"), function(outcome) {
  covars <- if (outcome == "DunedinPACE") "GP2_phenotype + sex + age" else "GP2_phenotype + sex"
  fit    <- lm(as.formula(paste(outcome, "~", covars)), data = clocks_cc)
  broom::tidy(fit) %>%
    filter(str_starts(term, "GP2_phenotype")) %>%
    mutate(outcome = outcome, .before = 1)
})
print(accel_results)

write.csv(accel_results,
          file.path(DIR_RESULTS, "clock_acceleration_by_diagnosis.csv"),
          row.names = FALSE)

cat("\nMethylation clock pipeline complete!\n")
cat("Results saved to:", DIR_RESULTS, "\n")
