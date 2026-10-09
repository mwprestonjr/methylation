# =============================================================================
# Clinical traits for clock associations
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 9, 2026
# Description: Builds one row per case with the traits in CLOCK_TRAITS
#              (clocks/config.R):
#                disease_duration   age at the blood draw - CLOCK_DURATION_FROM
#                                   (R12 master key)
#                family_history_pd  R12 master key (Yes/No)
#                pd_gene            PPMI only: monogenic PD gene from
#                                   PPMI_GENETICS (Idiopathic, LRRK2, GBA1,
#                                   PRKN, LRRK2+GBA1)
#                other traits       CLINICAL_DATA columns at the first of
#                                   CLOCK_CLINICAL_VISIT_MONTHS with a value
#              Traits a data source doesn't have are NA. Source after
#              config.R and clocks/config.R
# =============================================================================

# samples: GP2ID, clinical_id, age (age at the blood draw)
build_clock_traits <- function(samples) {
  traits <- samples %>% select(GP2ID, clinical_id, age)

  # --- R12 master key: disease duration, family history ---
  r12 <- read.csv(FNAME_METADATA, colClasses = "character", na.strings = c("", "NA")) %>%
    filter(GP2ID %in% traits$GP2ID) %>%
    distinct(GP2ID, .keep_all = TRUE) %>%
    transmute(GP2ID,
              duration_from     = as.numeric(.data[[CLOCK_DURATION_FROM]]),
              family_history_pd = ifelse(family_history_pd %in% c("Yes", "No"),
                                         family_history_pd, NA))
  traits <- traits %>%
    left_join(r12, by = "GP2ID") %>%
    mutate(disease_duration = age - duration_from, .after = age) %>%
    select(-duration_from)

  # --- PPMI monogenic PD status ---
  traits$pd_gene <- NA_character_
  if (!is.null(PPMI_GENETICS)) {
    genetics <- read.csv(PPMI_GENETICS, colClasses = "character") %>%
      transmute(clinical_id = PATNO,
                pd_gene = recode(VAR_GENE, "0" = "Idiopathic", "GBA" = "GBA1",
                                 "LRRK2, GBA" = "LRRK2+GBA1"))
    traits <- traits %>%
      select(-pd_gene) %>%
      left_join(genetics, by = "clinical_id")
  }

  # --- GP2 extended clinical data at the blood-draw visit ---
  clinical_cols <- setdiff(CLOCK_TRAITS$trait, names(traits))
  clinical <- readr::read_csv(CLINICAL_DATA,
                              col_types = readr::cols(.default = "c"),
                              col_select = any_of(c("GP2ID", "visit_month", clinical_cols)),
                              progress = FALSE) %>%
    filter(GP2ID %in% traits$GP2ID) %>%
    mutate(visit_month = as.numeric(visit_month),
           visit_rank  = match(visit_month, CLOCK_CLINICAL_VISIT_MONTHS)) %>%
    filter(!is.na(visit_rank)) %>%
    arrange(GP2ID, visit_rank) %>%
    group_by(GP2ID) %>%
    summarise(across(any_of(clinical_cols), ~ as.numeric(first(na.omit(.x)))), .groups = "drop")

  # Traits in CLOCK_TRAITS that aren't columns of the clinical data are NA
  missing <- setdiff(clinical_cols, names(clinical))
  if (length(missing)) cat("Traits not in the clinical data:", missing, "\n")

  traits %>%
    left_join(clinical, by = "GP2ID") %>%
    mutate(!!!setNames(rep(list(NA_real_), length(missing)), missing))
}
