# =============================================================================
# Overview: samples per data source
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Nov 5, 2026
# Description: Counts samples and participants per data source before and after
#              QC (OVERVIEW_STAGES), overall and by phenotype: phenotype, sex,
#              GP2 genetic ancestry (GP2_ANCESTRY_COL), self-reported race,
#              monogenic PD variant (where known), array, dataset, chips and
#              batches; and age. Fits no models and reads no M values. Not part
#              of `make`. Needs the master key and PPMI_GENETICS, and `make qc`
#              for every data source. Saves OVERVIEW_COUNTS, OVERVIEW_TABLE1 and
#              these figures in DIR_OVERVIEW, from OVERVIEW_MAIN_STAGE:
#   overview_01_phenotype   samples by phenotype
#   overview_02_sex         share of each sex by phenotype, with counts
#   overview_03_ancestry    samples by genetic ancestry
#   overview_04_age         age at sample collection by phenotype
# Usage:       Rscript overview/sample_overview.R
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

library(tidyverse)

# config.R needs a data source; all of them are loaded below
DATA_SOURCE <- "ppmi_p140"
source("config.R")
source("overview/config.R")
source("R/pd_genetics.R")

dir.create(DIR_OVERVIEW, showWarnings = FALSE, recursive = TRUE)
fig <- function(name) file.path(DIR_OVERVIEW, paste0("overview_", name, ".png"))

# --- 1. Samples --------------------------------------------------------------

if (!file.exists(GP2_MASTER_KEY)) {
  stop(GP2_MASTER_KEY, " not found (GP2_MASTER_KEY in config.R) - is the gp2_release12 bucket mounted?")
}
ancestry <- read.csv(GP2_MASTER_KEY, colClasses = "character", na.strings = c("", "NA")) %>%
  distinct(GP2ID, .keep_all = TRUE) %>%
  transmute(GP2ID, ancestry = .data[[GP2_ANCESTRY_COL]])

# One sample sheet of the current DATA_SOURCE, with ancestry and variant status
load_stage <- function(file, stage, pd_gene) {
  if (!file.exists(file)) stop(file, " not found - run make qc SOURCE=", DATA_SOURCE)
  sheet <- read.csv(file, colClasses = "character", na.strings = c("", "NA")) %>%
    select(all_of(SAMPLE_SHEET_COLUMNS)) %>%
    left_join(ancestry, by = "GP2ID") %>%
    left_join(pd_gene, by = "clinical_id")
  cat(DATA_SOURCE, stage, ":", nrow(sheet), "samples\n")
  sheet %>% mutate(data_source = DATA_SOURCE, stage = stage, age = as.numeric(age))
}

# Source config.R per data source for its own paths and PPMI_GENETICS
samples <- list()
for (data_source in DATA_SOURCES) {
  DATA_SOURCE <- data_source
  source("config.R")
  stage_files <- setNames(c(SAMPLE_SHEET, SAMPLE_SHEET_QC), c("before_qc", "after_qc"))
  stopifnot(all(OVERVIEW_STAGES %in% names(stage_files)))
  pd_gene <- load_pd_gene()
  for (stage in OVERVIEW_STAGES) {
    samples[[paste(data_source, stage)]] <- load_stage(stage_files[[stage]], stage, pd_gene)
  }
}

# Missing values get a label, so they are counted. Sources without variant data
# keep pd_gene NA and have no monogenic rows
samples <- bind_rows(samples) %>%
  group_by(data_source) %>%
  mutate(has_variants = any(!is.na(pd_gene))) %>%
  ungroup() %>%
  mutate(across(c(phenotype, sex, race, ancestry), ~ replace_na(.x, "Unknown")),
         pd_gene = ifelse(has_variants, replace_na(pd_gene, "No variant data"), NA))

if (anyNA(samples$age)) {
  cat("\nSamples without an age:\n")
  print(samples %>% filter(is.na(age)) %>% count(data_source, stage))
}

# --- 2. Counts ---------------------------------------------------------------

# Every row twice: under its phenotype and under "All"
by_phenotype <- bind_rows(samples, samples %>% mutate(phenotype = "All"))
groups <- c("data_source", "stage", "phenotype")

# Samples per level of a column, with the share of the group's known samples
count_levels <- function(data, variable, column) {
  data %>%
    filter(!is.na(.data[[column]])) %>%
    group_by(across(all_of(groups))) %>%
    mutate(n_group = n()) %>%
    group_by(across(all_of(c(groups, column))), n_group) %>%
    summarize(n = n(), .groups = "drop") %>%
    transmute(across(all_of(groups)), variable = variable, level = .data[[column]],
              n, pct = round(100 * n / n_group, 1))
}

# One number per group
count_total <- function(variable, value) {
  by_phenotype %>%
    group_by(across(all_of(groups))) %>%
    summarize(n = {{ value }}, .groups = "drop") %>%
    mutate(variable = variable, level = "Total", pct = NA_real_)
}

counts <- bind_rows(
  count_total("Samples", n()),
  count_total("Participants", n_distinct(GP2ID)),
  count_total("Chips", n_distinct(Sentrix_ID)),
  count_total("Batches", n_distinct(Batch)),
  count_levels(samples %>% mutate(label = phenotype, phenotype = "All"), "Phenotype", "label"),
  count_levels(by_phenotype, "Sex", "sex"),
  count_levels(by_phenotype, "Ancestry", "ancestry"),
  count_levels(by_phenotype, "Race", "race"),
  count_levels(by_phenotype, "Monogenic variant", "pd_gene"),
  count_levels(by_phenotype, "Array", "Array"),
  count_levels(by_phenotype, "Dataset", "Dataset")
) %>%
  arrange(data_source, match(stage, OVERVIEW_STAGES), phenotype != "All", phenotype) %>%
  select(all_of(groups), variable, level, n, pct)

write.csv(counts, OVERVIEW_COUNTS, row.names = FALSE)
cat("\nCounts saved to:", OVERVIEW_COUNTS, "\n")

cat("\nSamples before and after QC:\n")
counts %>%
  filter(phenotype == "All", variable %in% c("Samples", "Participants", "Phenotype", "Sex", "Ancestry")) %>%
  select(data_source, stage, variable, level, n) %>%
  pivot_wider(names_from = stage, values_from = n, values_fill = 0L) %>%
  mutate(removed = before_qc - after_qc) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# --- 3. Table 1 --------------------------------------------------------------

main <- by_phenotype %>% filter(stage == OVERVIEW_MAIN_STAGE)
age_rows <- main %>%
  group_by(data_source, phenotype) %>%
  summarize(`Mean (SD)` = sprintf("%.1f (%.1f)", mean(age, na.rm = TRUE), sd(age, na.rm = TRUE)),
            Range       = sprintf("%.1f to %.1f", min(age, na.rm = TRUE), max(age, na.rm = TRUE)),
            .groups = "drop") %>%
  pivot_longer(c(`Mean (SD)`, Range), names_to = "level", values_to = "cell") %>%
  mutate(variable = "Age")

table_variables <- c("Samples", "Participants", "Age", "Sex", "Ancestry", "Race", "Monogenic variant")
level_sizes <- counts %>%
  filter(stage == OVERVIEW_MAIN_STAGE, phenotype == "All") %>%
  group_by(variable, level) %>%
  summarize(size = sum(n), .groups = "drop")

# Columns: per data source, all samples, then phenotypes from largest to smallest
columns <- counts %>%
  filter(stage == OVERVIEW_MAIN_STAGE, variable == "Samples") %>%
  arrange(data_source, phenotype != "All", desc(n)) %>%
  transmute(column = paste(data_source, phenotype)) %>%
  pull(column)

# "-": none in that group, or not known for that data source
table1 <- counts %>%
  filter(stage == OVERVIEW_MAIN_STAGE, variable %in% table_variables) %>%
  mutate(cell = ifelse(is.na(pct), as.character(n), sprintf("%d (%.1f%%)", n, pct))) %>%
  select(data_source, phenotype, variable, level, cell) %>%
  bind_rows(age_rows) %>%
  left_join(level_sizes, by = c("variable", "level")) %>%
  arrange(match(variable, table_variables), desc(size)) %>%
  mutate(column = paste(data_source, phenotype)) %>%
  select(variable, level, column, cell) %>%
  pivot_wider(names_from = column, values_from = cell, values_fill = "-") %>%
  select(variable, level, all_of(columns))

write.csv(table1, OVERVIEW_TABLE1, row.names = FALSE)
cat("Table 1 saved to:", OVERVIEW_TABLE1, "\n")

# --- 4. Figures --------------------------------------------------------------

main_samples <- samples %>% filter(stage == OVERVIEW_MAIN_STAGE)

overview_theme <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

# One row per data source; bars keep one thickness however many levels it has
source_facets <- facet_grid(rows = vars(data_source), scales = "free_y", space = "free_y")

save_fig <- function(plot, name) {
  ggsave(fig(name), plot, width = FIG_WIDTH, height = FIG_HEIGHT, dpi = FIG_RES)
  cat("Saved:", fig(name), "\n")
}

# Counts of one column, as labeled horizontal bars per data source
count_bars <- function(column, x_label) {
  main_samples %>%
    count(data_source, level = .data[[column]]) %>%
    mutate(level = fct_reorder(level, n, .fun = sum)) %>%
    ggplot(aes(n, level)) +
    geom_col(fill = OVERVIEW_BAR_COLOR, width = 0.7) +
    geom_text(aes(label = n), hjust = -0.2, size = 2.8) +
    source_facets +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = x_label, y = NULL) +
    overview_theme +
    theme(panel.grid.major.y = element_blank())
}

save_fig(count_bars("phenotype", "Samples"), "01_phenotype")

# Phenotype labels with the group size, ordered by size
with_group_size <- function(data) {
  data %>%
    add_count(data_source, phenotype, name = "n_group") %>%
    mutate(group = fct_reorder(paste0(phenotype, " (n=", n_group, ")"), n_group))
}

sex_counts <- main_samples %>%
  with_group_size() %>%
  count(data_source, group, sex) %>%
  group_by(data_source, group) %>%
  mutate(share = n / sum(n)) %>%
  ungroup()
sex_plot <- ggplot(sex_counts, aes(share, group, fill = sex)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(share >= OVERVIEW_MIN_LABEL_SHARE, n, ""), color = sex),
            position = position_stack(vjust = 0.5), size = 2.8, show.legend = FALSE) +
  source_facets +
  scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.03))) +
  scale_fill_manual(values = OVERVIEW_SEX_COLORS) +
  scale_color_manual(values = OVERVIEW_SEX_TEXT) +
  labs(x = "Share of samples (count on the bar)", y = NULL, fill = NULL) +
  overview_theme +
  theme(panel.grid.major.y = element_blank(), legend.position = "bottom")
save_fig(sex_plot, "02_sex")

save_fig(count_bars("ancestry", "Samples"), "03_ancestry")

age_plot <- main_samples %>%
  filter(!is.na(age)) %>%
  with_group_size() %>%
  ggplot(aes(age, group)) +
  geom_jitter(height = 0.2, width = 0, size = 0.7, alpha = 0.4, color = OVERVIEW_BAR_COLOR) +
  geom_boxplot(fill = NA, outlier.shape = NA, width = 0.6, linewidth = 0.4) +
  source_facets +
  labs(x = "Age at sample collection (years)", y = NULL) +
  overview_theme +
  theme(panel.grid.major.y = element_blank())
save_fig(age_plot, "04_age")

cat("\nNext: Rscript overview/check_trait_data.R\n")
