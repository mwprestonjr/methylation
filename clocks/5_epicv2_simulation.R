# =============================================================================
# Clocks 5: EPICv2 simulation
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 29, 2026
# Description: Re-estimates the clocks of an EPICv1 data source using only
#              the CpGs that are also on EPICv2, as if the samples had been
#              run on EPICv2. Compares clock accuracy with the full-array
#              estimates (1_estimate.R), the same with age in whole years,
#              and the EPICv2 data source EPICV2_SIM_SOURCE; how closely the
#              EPICv2-CpG estimates agree with the full-array ones; and
#              refits the group comparison (2_compare.R main model) on them.
#              Separates the effect of the CpGs missing from EPICv2 from
#              other differences between the cohorts. Not part of
#              `make clocks`; run after 1_estimate.R and 2_compare.R for the
#              EPICv1 source, and 1_estimate.R for EPICV2_SIM_SOURCE
# Usage:       Rscript clocks/5_epicv2_simulation.R ppmi_p140
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------
# EPICv2 annotation first: it attaches minfi's dependencies (S4Vectors, ...),
# whose rename/select etc. must not mask dplyr's, which is attached after
source("R/array_profiles.R")
epicv2_pkg <- ARRAY_PROFILES$EPICv2$anno_pkg
library(epicv2_pkg, character.only = TRUE)
library(tidyverse)
library(methylclock)
library(DunedinPACE)

# Load shared configuration
source("config.R")
source("clocks/config.R")
source("R/clock_estimate.R")
source("R/clock_models.R")

# Full-array estimates and model covariates from 1_estimate.R
clocks_full <- read.csv(CLOCK_RESULTS,
                        colClasses = c(GP2ID = "character", combat_batch = "character"))
cat("Samples:", nrow(clocks_full), "\n")

# --- 1. Load unfiltered betas ------------------------------------------------

cat("\nLoading unfiltered beta values...\n")
bVals <- readRDS(BVALS_UNFILTERED)
cat("Beta values dimensions:", dim(bVals), "\n")
# EPICv2 probe names carry a suffix (cg00000029_TC21); EPICv1 names don't
if (any(grepl("_", rownames(bVals)))) {
  stop("5_epicv2_simulation.R needs an EPICv1 data source; ", DATA_SOURCE, " looks like EPICv2")
}
stopifnot(all(clocks_full$Sample %in% colnames(bVals)))

bVals_cg <- collapse_epicv2(bVals[, clocks_full$Sample])
rm(bVals); invisible(gc())

# --- 2. Restrict to the CpGs on EPICv2 ---------------------------------------

# Every CpG on the EPICv2 manifest, replicate-probe suffixes removed
epicv2_cpgs <- unique(sub("_.*$", "", rownames(minfi::getAnnotation(get(epicv2_pkg), what = "Locations"))))
cat("\nCpGs on EPICv2:", length(epicv2_cpgs), "\n")

clock_cpgs <- clock_cpg_list()

# Check: clock coverage from the EPICv2 manifest vs what the EPICv2 data
# source measured (its probes, unfiltered). Small differences are probes
# missing from its betas
comparison_coverage <- read.csv(file.path(EPICV2_SIM_COMPARISON, "clock_cpg_coverage.csv"))
manifest_check <- clock_coverage_table(clock_cpgs, epicv2_cpgs) %>%
  select(clock, n_manifest = n_on_array) %>%
  left_join(comparison_coverage %>% select(clock, n_data = n_on_array), by = "clock")
cat("\nClock CpGs on the EPICv2 manifest vs in", EPICV2_SIM_SOURCE, "betas:\n")
print(manifest_check)

# CpGs on both arrays: EPICv2 can't add the CpGs EPICv1 doesn't have, so
# coverage can be lower than on either array alone
bVals_sim <- bVals_cg[rownames(bVals_cg) %in% epicv2_cpgs, , drop = FALSE]
cat("\nCpGs kept (on EPICv1 data and EPICv2):", nrow(bVals_sim), "of", nrow(bVals_cg), "\n")
rm(bVals_cg); invisible(gc())

sim_coverage <- clock_coverage_table(clock_cpgs, rownames(bVals_sim))

# --- 3. Estimate the clocks from the EPICv2 CpGs -----------------------------

# Same settings as 1_estimate.R, including the coverage thresholds, so a clock
# below them is skipped here as it would be on EPICv2
cat("\nEstimating epigenetic age with methylclock (EPICv2 CpGs)...\n")
dnam_age <- DNAmAge(bVals_sim,
                    clocks     = "all",
                    normalize  = FALSE,
                    cell.count = FALSE,
                    min.perc   = CLOCK_COV_THRESH)

cat("\nEstimating DunedinPACE (EPICv2 CpGs)...\n")
pace <- PACEProjector(bVals_sim, proportionOfProbesRequired = PACE_COV_THRESH)

estimated <- clocks_estimated(clock_cpgs, dnam_age, pace)
sim_coverage <- sim_coverage %>% mutate(estimated = unname(estimated[clock]))
cat("\nClock CpG coverage, EPICv2 CpGs only:\n")
print(sim_coverage %>% mutate(pct_on_array = round(pct_on_array, 1)))

# Same table layout as CLOCK_RESULTS: clock ages and pace replaced, covariates
# kept, age acceleration recomputed
sim_outcomes <- intersect(ADULT_CLOCKS, names(dnam_age)[sapply(dnam_age, function(x) !all(is.na(x)))])
clocks_sim <- clocks_full %>%
  select(-any_of(c(names(dnam_age), "DunedinPACE")), -starts_with("AgeAccel_"), Sample) %>%
  left_join(dnam_age %>% rename(Sample = id), by = "Sample") %>%
  mutate(DunedinPACE = pace$DunedinPACE[Sample])
for (clock in sim_outcomes) {
  fit <- lm(clocks_sim[[clock]] ~ clocks_sim$age, na.action = na.exclude)
  clocks_sim[[paste0("AgeAccel_", clock)]] <- residuals(fit)
}

# --- 4. Accuracy: full array, EPICv2 CpGs, EPICv2 data source ----------------

# r, median absolute error and mean offset vs chronological age, as in
# 1_estimate.R. "whole-year age": the EPICv2-CpG estimates against age
# truncated to whole years, as recorded for the EPICv2 data source, to see how
# much age precision alone lowers the accuracy. age_sd is shown because r
# depends on the spread of ages, MAE much less
accuracy_of <- function(data, setting, age = data$age) {
  map_dfr(intersect(ADULT_CLOCKS, names(data)), function(clock) {
    diff <- data[[clock]] - age
    tibble(clock       = clock,
           setting     = setting,
           n           = sum(!is.na(diff)),
           age_sd      = sd(age[!is.na(diff)]),
           r           = cor(data[[clock]], age, use = "complete.obs"),
           MAE         = median(abs(diff), na.rm = TRUE),
           mean_offset = mean(diff, na.rm = TRUE))
  }) %>% filter(n > 0)
}

# The EPICv2 source's own analyses leave out the clocks not in its
# ADULT_CLOCKS, but all clocks were estimated and saved, so all are compared
comparison <- read.csv(file.path(EPICV2_SIM_COMPARISON, "clock_age_acceleration.csv"))
settings <- c(full = paste(DATA_SOURCE, "- all CpGs"),
              sim  = paste(DATA_SOURCE, "- EPICv2 CpGs"),
              sim_whole_years = paste(DATA_SOURCE, "- EPICv2 CpGs, whole-year age"),
              comparison = paste(EPICV2_SIM_SOURCE, "(EPICv2)"))
sim_accuracy <- bind_rows(
  accuracy_of(clocks_full, settings[["full"]]),
  accuracy_of(clocks_sim,  settings[["sim"]]),
  accuracy_of(clocks_sim,  settings[["sim_whole_years"]], age = floor(clocks_sim$age)),
  accuracy_of(comparison,  settings[["comparison"]])
) %>%
  left_join(bind_rows(
    read.csv(CLOCK_COVERAGE) %>% mutate(setting = settings[["full"]]),
    sim_coverage %>% mutate(setting = settings[["sim"]]),
    sim_coverage %>% mutate(setting = settings[["sim_whole_years"]]),
    comparison_coverage %>% mutate(setting = settings[["comparison"]])
  ) %>% select(clock, setting, pct_on_array), by = c("clock", "setting")) %>%
  mutate(setting = factor(setting, levels = settings), clock = factor(clock, levels = ADULT_CLOCKS)) %>%
  arrange(clock, setting)
cat("\nClock accuracy vs chronological age:\n")
print(sim_accuracy %>% mutate(across(where(is.numeric), ~ round(.x, 2))), n = Inf)

# --- 5. Agreement with the full-array estimates ------------------------------

# Per outcome: correlation of the EPICv2-CpG and full-array values across
# samples (do they rank people the same way?) and the mean and SD of the
# difference. For AgeAccel the mean difference is ~0 by construction (both are
# residuals); the shift in clock age shows in the accuracy offsets above
agreement_outcomes <- c(paste0("AgeAccel_", sim_outcomes), "DunedinPACE")
sim_agreement <- map_dfr(agreement_outcomes, function(outcome) {
  full <- clocks_full[[outcome]]
  sim  <- clocks_sim[[outcome]][match(clocks_full$Sample, clocks_sim$Sample)]
  tibble(outcome   = outcome,
         n         = sum(!is.na(full) & !is.na(sim)),
         r         = cor(full, sim, use = "complete.obs"),
         mean_diff = mean(sim - full, na.rm = TRUE),
         sd_diff   = sd(sim - full, na.rm = TRUE),
         sd_full   = sd(full, na.rm = TRUE))
})
cat("\nEPICv2-CpG vs full-array estimates:\n")
print(sim_agreement %>% mutate(across(where(is.numeric), ~ round(.x, 3))))

# --- 6. Group comparison on the EPICv2-CpG estimates -------------------------

# 2_compare.R main model, so the effects can be set against CLOCK_GROUP_EFFECTS
clocks_cmp <- clocks_sim %>%
  filter(phenotype %in% CLOCK_GROUPS) %>%
  mutate(phenotype = factor(phenotype, levels = CLOCK_GROUPS))
group_sizes <- table(clocks_cmp$phenotype)
if (all(group_sizes >= CLOCK_MIN_GROUP_N)) {
  usable <- usable_covariates(clocks_cmp, CLOCK_COVARIATES)
  cat("\nGroup comparison, EPICv2 CpGs - model: <outcome> ~",
      paste(c("phenotype", usable), collapse = " + "), "\n")
  sim_group_effects <- fit_predictor(clocks_cmp, estimated_outcomes(clocks_cmp, agreement_outcomes),
                                     "phenotype", usable) %>%
    mutate(fdr = p.adjust(p.value, method = "BH"))
  print(sim_group_effects %>% select(outcome, n, term, estimate, std.error, p.value, fdr))
  write.csv(sim_group_effects, EPICV2_SIM_GROUP_EFFECTS, row.names = FALSE)
} else {
  cat("\nFewer than", CLOCK_MIN_GROUP_N, "samples in a group; skipping the group comparison\n")
  if (file.exists(EPICV2_SIM_GROUP_EFFECTS)) file.remove(EPICV2_SIM_GROUP_EFFECTS)
}

write.csv(clocks_sim, EPICV2_SIM_RESULTS, row.names = FALSE)
write.csv(sim_coverage, EPICV2_SIM_COVERAGE, row.names = FALSE)
write.csv(sim_accuracy, EPICV2_SIM_ACCURACY, row.names = FALSE)
write.csv(sim_agreement, EPICV2_SIM_AGREEMENT, row.names = FALSE)
cat("\nEPICv2 simulation results saved to:", DIR_RESULTS_CLOCKS, "\n")

# --- 7. Figure: EPICv2-CpG vs full-array estimates ---------------------------

agreement_long <- map_dfr(agreement_outcomes, function(o) {
  tibble(outcome   = sub("^AgeAccel_", "", o),
         full      = clocks_full[[o]],
         sim       = clocks_sim[[o]][match(clocks_full$Sample, clocks_sim$Sample)],
         phenotype = clocks_full$phenotype)
}) %>%
  mutate(outcome = factor(outcome, levels = c(sim_outcomes, "DunedinPACE")))
agreement_labels <- sim_agreement %>%
  mutate(outcome = factor(sub("^AgeAccel_", "", outcome), levels = levels(agreement_long$outcome)),
         label   = sprintf("r = %.2f", r))

png(file.path(DIR_FIGURES_CLOCKS, "clock_05_epicv2_simulation.png"),
    width = FIG_WIDTH * 1.3, height = FIG_HEIGHT * 1.3, units = "in", res = FIG_RES)
print(ggplot(agreement_long, aes(full, sim)) +
  geom_abline(linetype = "dashed", color = "gray50") +
  geom_point(aes(color = phenotype), alpha = 0.6, size = 0.8) +
  geom_text(data = agreement_labels, aes(label = label), x = -Inf, y = Inf,
            hjust = -0.1, vjust = 1.3, size = 3) +
  facet_wrap(~ outcome, scales = "free") +
  labs(x = "All CpGs", y = "EPICv2 CpGs only", color = "Phenotype",
       title = paste("Clock estimates from EPICv2 CpGs only -", DATA_SOURCE),
       caption = paste("Age acceleration (years) and DunedinPACE (years per year).",
                       "Dashed: identical estimates.", sep = "\n")) +
  theme_bw(base_size = 9) +
  theme(plot.caption = element_text(hjust = 0), plot.caption.position = "plot"))
dev.off()

cat("\nEPICv2 simulation complete!\n")
