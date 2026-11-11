# =============================================================================
# Configuration - overview module (overview/ scripts)
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Nov 5, 2026
# Updated: Nov 11, 2026
# Description: Settings used only by overview/sample_overview.R and
#              overview/array_comparison.R. Source after config.R. Outputs
#              cover all data sources, so they sit next to the per-source
#              output folders rather than in one of them
# =============================================================================

# --- Stages ------------------------------------------------------------------

# Sample sheets counted, in pipeline order. Preprocessing removes no samples, so
# SAMPLE_SHEET_FINAL would repeat the second one. Table 1 and the figures use
# OVERVIEW_MAIN_STAGE, the samples later analyses start from
OVERVIEW_STAGES     <- c("before_qc", "after_qc")
OVERVIEW_MAIN_STAGE <- "after_qc"

# --- Figures -----------------------------------------------------------------

# Bars of one series, and sex. Colorblind-safe pair; labels printed on a
# segment use OVERVIEW_SEX_TEXT, the readable ink for that fill
OVERVIEW_BAR_COLOR <- "#2a78d6"
OVERVIEW_SEX_COLORS <- c(Female = "#2a78d6", Male = "#eb6834", Unknown = "#b3b2ab")
OVERVIEW_SEX_TEXT   <- c(Female = "white", Male = "black", Unknown = "black")

# Segments below this share of a bar get no count: it wouldn't fit. The count
# is still in OVERVIEW_COUNTS
OVERVIEW_MIN_LABEL_SHARE <- 0.1

# --- Array comparison (array_comparison.R) -----------------------------------

# Island context colors in overview_06 (Island, Shore, Shelf, OpenSea), from the
# same palette; segments are labeled, since the green and yellow are light
OVERVIEW_ISLAND_COLORS <- c(Island = "#2a78d6", Shore = "#eb6834", Shelf = "#1baf7a", OpenSea = "#eda100")

# CpGs the mQTL scan tested, per ancestry folder of DIR_RESULTS_MQTL (written by
# mqtl/3_mqtl.R); ancestries without it are skipped
OVERVIEW_MQTL_TESTED_FILE <- "mqtl_cpgs_tested.tsv.gz"

# --- Outputs -----------------------------------------------------------------

DIR_OVERVIEW     <- file.path(dirname(DIR_OUTPUT), "reports", "overview")
OVERVIEW_COUNTS  <- file.path(DIR_OVERVIEW, "overview_counts.csv")   # long: data source x stage x phenotype x variable x level
OVERVIEW_TABLE1  <- file.path(DIR_OVERVIEW, "overview_table1.csv")   # formatted summary of OVERVIEW_MAIN_STAGE, one column per data source x phenotype
# array_comparison.R
OVERVIEW_ARRAY_SUMMARY <- file.path(DIR_OVERVIEW, "array_comparison_summary.csv")   # section x metric: probe and site counts, overlap, replicates, QC and mQTL overlap
OVERVIEW_ARRAY_CONTEXT <- file.path(DIR_OVERVIEW, "array_comparison_context.csv")   # overlap category x island context / gene region: n and percent
OVERVIEW_ARRAY_SITES   <- file.path(DIR_OVERVIEW, "array_comparison_sites.csv.gz")  # one row per CpG site on either array: presence, probes, design, position, QC and mQTL status
