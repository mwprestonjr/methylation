# =============================================================================
# Monogenic PD variants
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 30, 2026
# Description: Monogenic PD variant status per participant, shared by
#              clocks/ (pd_gene trait, R/clock_traits.R) and ewas/ (carrier
#              models, ewas/1_case_control.R). Source after config.R
# =============================================================================

# One row per participant in PPMI_GENETICS (config.R), by clinical_id (PPMI
# PATNO): pd_gene None (no variant), LRRK2, GBA1, PRKN, LRRK2+GBA1, ... Data
# sources without variant data (PPMI_GENETICS NULL) give no rows, so a join
# leaves pd_gene NA
load_pd_gene <- function() {
  if (is.null(PPMI_GENETICS)) {
    return(tibble(clinical_id = character(), pd_gene = character()))
  }
  if (!file.exists(PPMI_GENETICS)) stop(PPMI_GENETICS, " not found (PPMI_GENETICS in config.R)")
  read.csv(PPMI_GENETICS, colClasses = "character") %>%
    transmute(clinical_id = PATNO,
              pd_gene = recode(VAR_GENE, "0" = "None", "GBA" = "GBA1",
                               "LRRK2, GBA" = "LRRK2+GBA1"))
}
