# =============================================================================
# Clocks 3: associations with clinical traits in cases
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 9, 2026
# Description: Within CLOCK_CASE_GROUP, tests each trait in CLOCK_TRAITS
#              (clocks/config.R) against age acceleration and DunedinPACE:
#                <outcome> ~ trait + CLOCK_COVARIATES
#              one model per outcome and trait. Traits are measured at the
#              blood draw (R/clock_traits.R). Traits with fewer than
#              CLOCK_MIN_TRAIT_N cases are skipped; categorical levels with
#              fewer than CLOCK_MIN_LEVEL_N cases are left out. FDR (BH) is
#              across all outcome x trait tests. Saves the trait values, a
#              per-trait summary and the effects for 4_plots.R
# Usage:       Rscript clocks/3_associations.R <data source>   (see config.R)
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------
library(tidyverse)

# Load shared configuration
source("config.R")
source("clocks/config.R")
source("R/clock_models.R")
source("R/clock_traits.R")

# --- 1. Cases and their traits -----------------------------------------------

# combat_batch is a label (plate or chip number), not a number
clocks <- read.csv(CLOCK_RESULTS,
                   colClasses = c(GP2ID = "character", combat_batch = "character")) %>%
  filter(phenotype == CLOCK_CASE_GROUP)
cat("Cases (", CLOCK_CASE_GROUP, "):", nrow(clocks), "\n")

# clinical_id (e.g. PPMI PATNO) is needed for source-specific trait files
ids <- read.csv(SAMPLE_SHEET_FINAL, colClasses = "character") %>%
  distinct(GP2ID, clinical_id)

cat("\nBuilding trait table...\n")
traits <- build_clock_traits(clocks %>% select(GP2ID, age) %>% left_join(ids, by = "GP2ID"))
write.csv(traits, CLOCK_TRAIT_VALUES, row.names = FALSE)

# --- 2. Prepare each trait ---------------------------------------------------

# Returns the trait as a model variable (numeric, or factor with the reference
# level first and small levels set to NA) and a summary row saying whether it
# is tested and why not
prepare_trait <- function(x, type, ref) {
  if (type == "continuous") {
    x <- suppressWarnings(as.numeric(x))
    n <- sum(!is.na(x))
    status <- case_when(n < CLOCK_MIN_TRAIT_N ~ paste("skipped: fewer than", CLOCK_MIN_TRAIT_N, "cases"),
                        sd(x, na.rm = TRUE) %in% c(0, NA) ~ "skipped: no variation",
                        TRUE ~ "tested")
    distribution <- if (n == 0) "no values" else
      sprintf("mean %.2f, SD %.2f, range %.1f-%.1f", mean(x, na.rm = TRUE),
              sd(x, na.rm = TRUE), min(x, na.rm = TRUE), max(x, na.rm = TRUE))
    summary <- tibble(n = n, distribution = distribution,
                      trait_sd = sd(x, na.rm = TRUE), status = status)
    return(list(x = x, summary = summary))
  }

  counts  <- table(x)
  dropped <- names(counts)[counts < CLOCK_MIN_LEVEL_N]
  x[x %in% dropped] <- NA
  levels  <- setdiff(names(counts), dropped)
  levels  <- c(intersect(ref, levels), setdiff(levels, ref))   # reference first
  x <- factor(x, levels = levels)
  n <- sum(!is.na(x))
  status <- case_when(n < CLOCK_MIN_TRAIT_N ~ paste("skipped: fewer than", CLOCK_MIN_TRAIT_N, "cases"),
                      !ref %in% levels ~ paste0("skipped: reference level '", ref, "' missing or too small"),
                      length(levels) < 2 ~ "skipped: fewer than 2 levels",
                      TRUE ~ "tested")
  distribution <- if (length(counts) == 0) "no values" else
    paste(paste0(names(counts), " ", counts), collapse = ", ")
  if (length(dropped)) {
    distribution <- paste0(distribution, " (left out: ", paste(dropped, collapse = ", "), ")")
  }
  list(x = x, summary = tibble(n = n, distribution = distribution, trait_sd = NA_real_, status = status))
}

data     <- clocks %>% left_join(traits %>% select(-age, -clinical_id), by = "GP2ID")
outcomes <- estimated_outcomes(data, c(paste0("AgeAccel_", ADULT_CLOCKS), "DunedinPACE"))

trait_summary <- tibble()
for (i in seq_len(nrow(CLOCK_TRAITS))) {
  t <- CLOCK_TRAITS[i, ]
  p <- prepare_trait(data[[t$trait]], t$type, t$ref)
  data[[t$trait]] <- p$x
  trait_summary <- bind_rows(trait_summary,
                             p$summary %>% mutate(trait = t$trait, type = t$type, label = t$label, .before = 1))
}
cat("\nTraits:\n")
print(as.data.frame(trait_summary %>% select(trait, n, status, distribution)), right = FALSE)
write.csv(trait_summary, CLOCK_TRAIT_SUMMARY, row.names = FALSE)

# --- 3. Models ---------------------------------------------------------------

tested <- trait_summary %>% filter(status == "tested")

if (nrow(tested) > 0) {
  effects <- pmap_dfr(tested %>% select(trait, label, type, trait_sd), function(trait, label, type, trait_sd) {
    # Covariates are chosen among the cases with this trait
    covars <- usable_covariates(data %>% filter(!is.na(.data[[trait]])), CLOCK_COVARIATES)
    fit_predictor(data, outcomes, trait, covars) %>%
      mutate(trait = trait, label = label, type = type, trait_sd = trait_sd,
             level = ifelse(type == "categorical", str_remove(term, fixed(trait)), NA),
             .after = outcome)
  }) %>%
    mutate(fdr = p.adjust(p.value, method = "BH"))

  cat("\nAssociations (", nrow(effects), " tests; FDR across all):\n", sep = "")
  print(as.data.frame(effects %>%
                        select(outcome, trait, level, n, estimate, p.value, fdr) %>%
                        mutate(across(c(estimate), ~ signif(.x, 3)),
                               across(c(p.value, fdr), ~ signif(.x, 2)))),
        right = FALSE)
  cat("\nFDR < 0.05:", sum(effects$fdr < 0.05), "of", nrow(effects), "\n")
  write.csv(effects, CLOCK_TRAIT_EFFECTS, row.names = FALSE)
  cat("Associations saved to:", CLOCK_TRAIT_EFFECTS, "\n")
} else {
  cat("\nNo traits with enough data; skipping the associations\n")
  # Remove results from an earlier run so 4_plots.R doesn't draw stale ones
  if (file.exists(CLOCK_TRAIT_EFFECTS)) file.remove(CLOCK_TRAIT_EFFECTS)
}

cat("\nNext: Rscript clocks/4_plots.R", DATA_SOURCE, "\n")
