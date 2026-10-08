# =============================================================================
# Clocks 3: plots
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 8, 2026
# Description: Figures from the tables saved by 1_estimate.R and
#              2_compare.R, in DIR_FIGURES:
#   clock_01_vs_chronological_age   clock age vs chronological age per clock,
#                                   with r, median absolute error and the
#                                   share of the clock's CpGs on the array
#   clock_02_acceleration_by_group  age acceleration and DunedinPACE by
#                                   phenotype (all phenotypes in the data)
#   clock_03_group_effects          adjusted group difference (CLOCK_GROUPS)
#                                   per outcome with 95% CI; only if
#                                   2_compare.R ran the comparison
# Usage:       Rscript clocks/3_plots.R <data source>   (see config.R)
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------
library(tidyverse)

# Load shared configuration
source("config.R")
source("clocks/config.R")

fig <- function(name) file.path(DIR_FIGURES, paste0("clock_", name, ".png"))

clocks         <- read.csv(CLOCK_RESULTS)
clock_accuracy <- read.csv(CLOCK_ACCURACY)
clock_coverage <- read.csv(CLOCK_COVERAGE)
cat("Samples:", nrow(clocks), "\n")

# --- 1. Clock age vs chronological age ---------------------------------------

# Shared square axes so distance from y = x is comparable across clocks
clocks_long <- clocks %>%
  pivot_longer(all_of(ADULT_CLOCKS), names_to = "clock", values_to = "DNAm_age")
axis_lim <- range(c(clocks_long$age, clocks_long$DNAm_age), na.rm = TRUE)

accuracy_labels <- clock_accuracy %>%
  left_join(clock_coverage %>% select(clock, pct_on_array), by = "clock") %>%
  mutate(label = sprintf("r = %.2f\nMAE = %.1f\nCpGs = %.0f%%", r, MAE, pct_on_array))

png(fig("01_vs_chronological_age"), width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
print(ggplot(clocks_long, aes(age, DNAm_age)) +
  geom_abline(linetype = "dashed", colour = "grey50") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.6) +
  geom_point(aes(colour = phenotype), alpha = 0.7) +
  geom_text(data = accuracy_labels, aes(label = label),
            x = axis_lim[1], y = axis_lim[2], hjust = 0, vjust = 1, size = 3) +
  facet_wrap(~ clock) +
  coord_equal(xlim = axis_lim, ylim = axis_lim) +
  labs(x = "Chronological age", y = "Epigenetic age", colour = "Phenotype",
       caption = "Dashed: epigenetic = chronological age. Solid: linear fit. MAE: median absolute error (years).\nCpGs: share of the clock's CpGs on the array; missing ones are left out of the clock.") +
  theme_bw() +
  theme(plot.caption = element_text(hjust = 0), plot.caption.position = "plot"))
dev.off()

# --- 2. Age acceleration by phenotype ----------------------------------------

# Unadjusted: AgeAccel (years) and DunedinPACE per phenotype, so groups not in
# CLOCK_GROUPS (e.g. DLB, MSA, PSP) are visible too
accel_long <- clocks %>%
  select(phenotype, all_of(paste0("AgeAccel_", ADULT_CLOCKS)), DunedinPACE) %>%
  pivot_longer(-phenotype, names_to = "outcome", values_to = "value") %>%
  mutate(outcome = factor(sub("^AgeAccel_", "", outcome),
                          levels = c(ADULT_CLOCKS, "DunedinPACE")))

png(fig("02_acceleration_by_group"), width = FIG_WIDTH * 1.5, height = FIG_HEIGHT * 1.3,
    units = "in", res = FIG_RES)
print(ggplot(accel_long, aes(phenotype, value)) +
  geom_hline(data = data.frame(outcome = factor(c(ADULT_CLOCKS, "DunedinPACE"),
                                                levels = c(ADULT_CLOCKS, "DunedinPACE")),
                               ref = c(rep(0, length(ADULT_CLOCKS)), 1)),
             aes(yintercept = ref), linetype = "dashed", colour = "grey50") +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.4) +
  facet_wrap(~ outcome, scales = "free_y") +
  labs(title = paste("Age acceleration by phenotype -", DATA_SOURCE),
       subtitle = "Clocks: years older (+) or younger (-) than expected for age. DunedinPACE: years of aging per year",
       x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()

# --- 3. Adjusted group effects -----------------------------------------------

if (file.exists(CLOCK_GROUP_EFFECTS)) {
  effects <- read.csv(CLOCK_GROUP_EFFECTS) %>%
    mutate(outcome = factor(sub("^AgeAccel_", "", outcome),
                            levels = rev(c(ADULT_CLOCKS, "DunedinPACE"))),
           unit    = ifelse(outcome == "DunedinPACE", "DunedinPACE (years per year)",
                            "Age acceleration (years)"))
  png(fig("03_group_effects"), width = FIG_WIDTH * 1.4, height = FIG_HEIGHT,
      units = "in", res = FIG_RES)
  print(ggplot(effects, aes(estimate, outcome)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
    geom_point(size = 2) +
    geom_text(aes(label = sprintf("p = %.2g", p.value)), vjust = -0.8, size = 2.5) +
    facet_wrap(~ unit, scales = "free") +
    labs(title = paste0(paste(rev(CLOCK_GROUPS), collapse = " vs "), " - ", DATA_SOURCE),
         subtitle = paste("Adjusted difference (95% CI); model:", effects$covariates[1]),
         x = paste("Difference,", CLOCK_GROUPS[2], "minus", CLOCK_GROUPS[1]), y = NULL) +
    theme_bw(base_size = 8))
  dev.off()
} else {
  cat("No group comparison results (", basename(CLOCK_GROUP_EFFECTS),
      "); skipping the group effects plot\n")
}

cat("\nClock plots saved to:", DIR_FIGURES, "\n")
