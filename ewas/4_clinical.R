# =============================================================================
# EWAS 4: clinical traits
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Nov 4, 2026
# Description: CpG-wise limma models of the ComBat-corrected M values
#              (preprocessing/3_combat.R) on each trait in EWAS_CLINICAL_TRAITS
#              (ewas/config.R), within cases (EWAS_GROUPS[2]), adjusted for
#              EWAS_COVARIATES and disease duration. Traits are from the GP2
#              clinical data at the blood-draw visit, or the GP2 master key. A
#              continuous trait gives one result set (clinical_<name>); a
#              categorical one gives one per level vs its reference
#              (clinical_<name>_<level>), all from the same model. Saved in
#              DIR_RESULTS_EWAS, with the per-sample values
#              (EWAS_CLINICAL_VALUES). A trait is skipped if it has too little
#              data
# Usage:       Rscript ewas/4_clinical.R <data source>   (see config.R)
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

library(tidyverse)
library(limma)

# Load shared configuration and EWAS helpers
source("config.R")
source("ewas/config.R")
source("R/ewas_functions.R")

# Clinical outputs from earlier runs (possibly of another trait or level) would
# otherwise be picked up by 5_plots.R
remove_ewas_outputs("clinical_.*")

inp <- load_ewas_inputs()

# --- 1. Trait values and disease duration ------------------------------------

clinical_items <- unique(EWAS_CLINICAL_TRAITS$column[EWAS_CLINICAL_TRAITS$file == "clinical"])
master_items   <- unique(c(EWAS_CLINICAL_TRAITS$column[EWAS_CLINICAL_TRAITS$file == "master"],
                           DISEASE_DURATION_FROM))

cat("\nLoading", paste(clinical_items, collapse = ", "), "from", CLINICAL_DATA, "...\n")
cat("Loading", paste(master_items, collapse = ", "), "from", GP2_MASTER_KEY, "...\n")
inp$targets <- inp$targets %>%
  left_join(load_clinical_items(inp$targets$GP2ID, clinical_items), by = "GP2ID") %>%
  left_join(load_master_items(inp$targets$GP2ID, master_items), by = "GP2ID") %>%
  mutate(disease_duration = age - as.numeric(.data[[DISEASE_DURATION_FROM]]))
stopifnot(identical(inp$targets$Sample, colnames(inp$mVals)))

is_case <- inp$targets$phenotype %in% EWAS_GROUPS[2]
cat("\nCases (", EWAS_GROUPS[2], "):", sum(is_case), "\n")
cat("Disease duration in cases (years, from ", DISEASE_DURATION_FROM, "):\n", sep = "")
print(summary(inp$targets$disease_duration[is_case]))

write.csv(inp$targets %>%
            filter(is_case) %>%
            select(Sample, GP2ID, phenotype, disease_duration,
                   all_of(unique(EWAS_CLINICAL_TRAITS$column))),
          EWAS_CLINICAL_VALUES, row.names = FALSE)
cat("Clinical values saved to:", EWAS_CLINICAL_VALUES, "\n")

# --- 2. Models ---------------------------------------------------------------

# Within cases, so phenotype isn't a covariate
for (i in seq_len(nrow(EWAS_CLINICAL_TRAITS))) {
  spec   <- EWAS_CLINICAL_TRAITS[i, ]
  values <- replace(inp$targets[[spec$column]], !is_case, NA)
  fit_trait_model(inp, spec, values, "clinical", extra = "disease_duration")
}

cat("\nNext: Rscript ewas/5_plots.R", DATA_SOURCE, "\n")
