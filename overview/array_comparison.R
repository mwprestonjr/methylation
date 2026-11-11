# =============================================================================
# Overview: EPIC array comparison
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Nov 11, 2026
# Description: Compares the CpG content of the two arrays used by the data
#              sources (EPICv1 and EPICv2), for planning analyses that combine
#              them (meta-analysis, cross-ancestry mQTL). A site is a probe
#              name without its EPICv2 suffix (cg00002033_TC11 -> cg00002033),
#              as QC averages replicate probes (R/mvalues.R). Counts, from the
#              array manifests (minfi annotation packages, ARRAY_PROFILES):
#              probes by type, sites, EPICv2 sites with several probes, sites
#              on both arrays or one, other ways of matching them (EPICv1_Loci
#              column, hg38 position), shared sites with a different probe
#              design or position, and 450K sites. Then, per data source,
#              sites passing QC (MVALS_QC) and, where mqtl/ has run, sites the
#              mQTL scan tested; and the island and gene context of shared and
#              array-specific sites (CPG_ANNOTATION). Needs `make
#              preprocessing` for every data source. Not part of `make`. Saves
#              OVERVIEW_ARRAY_SUMMARY, OVERVIEW_ARRAY_CONTEXT,
#              OVERVIEW_ARRAY_SITES and these figures in DIR_OVERVIEW:
#   overview_05_array_overlap   sites on both arrays or one: manifests, after
#                               QC, tested in the mQTL scan
#   overview_06_array_context   island context of shared and array-specific
#                               sites
# Usage:       Rscript overview/array_comparison.R
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

# Annotation packages are loaded before tidyverse (they mask dplyr verbs)
suppressMessages(library(minfi))
library(tidyverse)

# config.R needs a data source; all of them are loaded below
DATA_SOURCE <- "ppmi_p140"
source("config.R")
source("overview/config.R")
source("R/array_profiles.R")

dir.create(DIR_OVERVIEW, showWarnings = FALSE, recursive = TRUE)
fig <- function(name) file.path(DIR_OVERVIEW, paste0("overview_", name, ".png"))
to_site <- function(probe) sub("_.*$", "", probe)

# --- 1. Data sources and their arrays ----------------------------------------

# Source config.R per data source for its own paths
sources <- list()
for (data_source in DATA_SOURCES) {
  DATA_SOURCE <- data_source
  source("config.R")
  for (f in c(SAMPLE_SHEET_QC, MVALS_QC, CPG_ANNOTATION)) {
    if (!file.exists(f)) stop(f, " not found - run make preprocessing SOURCE=", data_source)
  }
  array <- unique(read.csv(SAMPLE_SHEET_QC, colClasses = "character")$Array)
  if (length(array) != 1) stop(data_source, " has more than one array: ", paste(array, collapse = ", "))
  mqtl_dirs <- list.dirs(DIR_RESULTS_MQTL, recursive = FALSE)
  mqtl_dirs <- mqtl_dirs[file.exists(file.path(mqtl_dirs, OVERVIEW_MQTL_TESTED_FILE))]
  sources[[data_source]] <- list(array = array, mvals = MVALS_QC, annotation = CPG_ANNOTATION,
                                 mqtl = setNames(file.path(mqtl_dirs, OVERVIEW_MQTL_TESTED_FILE),
                                                 basename(mqtl_dirs)))
  cat(data_source, ":", array, "\n")
}

# Arrays in data source order (the first is the reference, e.g. EPICv1)
arrays <- unique(map_chr(sources, "array"))
if (length(arrays) != 2) stop("Expects two arrays across DATA_SOURCES, found: ", paste(arrays, collapse = ", "))
a1 <- arrays[1]
a2 <- arrays[2]
# One data source per array supplies its CPG_ANNOTATION
annotation_source <- map_chr(set_names(arrays), ~ names(sources)[map_chr(sources, "array") == .x][1])

# --- 2. Array manifests ------------------------------------------------------

# Every probe of an array's minfi annotation (SNP rs probes aren't in it)
load_manifest <- function(array) {
  pkg <- ARRAY_PROFILES[[array]]$anno_pkg
  cat("Loading", pkg, "...\n")
  suppressMessages(library(pkg, character.only = TRUE))
  ann <- getAnnotation(get(pkg))
  tibble(probe     = ann$Name,
         site      = to_site(ann$Name),
         probe_type = substr(ann$Name, 1, 2),          # cg, ch, nv
         design    = ann$Type,                         # Infinium I or II
         on_450k   = !is.na(ann$Methyl450_Loci) & !ann$Methyl450_Loci %in% c("", "FALSE"),
         epicv1_loci = if ("EPICv1_Loci" %in% colnames(ann)) ann$EPICv1_Loci else NA_character_)
}
manifests <- map(set_names(arrays), load_manifest)

# One row per site: probe count, design(s), 450K
per_site <- function(manifest, array) {
  manifest %>%
    group_by(site) %>%
    summarize(probe_type = first(probe_type),
              n_probes   = n(),
              design     = paste(sort(unique(design)), collapse = "/"),
              on_450k    = any(on_450k), .groups = "drop") %>%
    rename_with(~ paste0(.x, "_", array), -site)
}

# hg38 position per site (CPG_ANNOTATION; the first probe of an EPICv2 site)
load_positions <- function(array) {
  read.csv(sources[[annotation_source[[array]]]]$annotation, colClasses = c(chr = "character")) %>%
    mutate(site = to_site(cpg)) %>%
    distinct(site, .keep_all = TRUE) %>%
    transmute(site, chr, pos_hg38, island, gene_region) %>%
    rename_with(~ paste0(.x, "_", array), -site)
}

sites <- full_join(per_site(manifests[[a1]], a1), per_site(manifests[[a2]], a2), by = "site") %>%
  mutate(!!paste0("on_", a1) := !is.na(.data[[paste0("n_probes_", a1)]]),
         !!paste0("on_", a2) := !is.na(.data[[paste0("n_probes_", a2)]]),
         probe_type = coalesce(.data[[paste0("probe_type_", a2)]], .data[[paste0("probe_type_", a1)]]),
         category = case_when(.data[[paste0("on_", a1)]] & .data[[paste0("on_", a2)]] ~ "Both arrays",
                              .data[[paste0("on_", a1)]] ~ paste(a1, "only"),
                              TRUE ~ paste(a2, "only"))) %>%
  select(-starts_with("probe_type_")) %>%
  left_join(load_positions(a1), by = "site") %>%
  left_join(load_positions(a2), by = "site")

on1 <- sites[[paste0("on_", a1)]]
on2 <- sites[[paste0("on_", a2)]]
both <- on1 & on2
same_position <- both &
  sites[[paste0("chr_", a1)]] %in% sites[[paste0("chr_", a2)]] &
  sites[[paste0("chr_", a1)]] == sites[[paste0("chr_", a2)]] &
  sites[[paste0("pos_hg38_", a1)]] == sites[[paste0("pos_hg38_", a2)]]
positioned <- both & !is.na(sites[[paste0("pos_hg38_", a1)]]) & !is.na(sites[[paste0("pos_hg38_", a2)]])

# --- 3. Summary --------------------------------------------------------------

summary_rows <- list()
add <- function(section, metric, value, note = "") {
  summary_rows[[length(summary_rows) + 1]] <<- tibble(section = section, metric = metric,
                                                      value = value, note = note)
}

for (array in arrays) {
  m <- manifests[[array]]
  add("Manifest", paste(array, "probes"), nrow(m))
  for (type in sort(unique(m$probe_type))) add("Manifest", paste(array, type, "probes"), sum(m$probe_type == type))
  add("Manifest", paste(array, "sites"), n_distinct(m$site), "probe name without suffix")
  add("Manifest", paste(array, "sites on 450K"), sum(sites[[paste0("on_450k_", array)]] %in% TRUE))
}

# Sites measured by more than one probe (EPICv2 replicates)
for (array in arrays) {
  n_probes <- sites[[paste0("n_probes_", array)]]
  if (any(n_probes > 1, na.rm = TRUE)) {
    add("Replicate probes", paste(array, "sites with more than one probe"), sum(n_probes > 1, na.rm = TRUE),
        "averaged into one value per site in QC")
    add("Replicate probes", paste(array, "probes at those sites"), sum(n_probes[n_probes > 1], na.rm = TRUE))
    for (k in sort(unique(n_probes[n_probes > 1]))) {
      add("Replicate probes", paste(array, "sites with", k, "probes"), sum(n_probes == k, na.rm = TRUE))
    }
  }
}

add("Overlap", "Sites on both arrays", sum(both), "matched by site name")
add("Overlap", paste("Sites on", a1, "only"), sum(on1 & !on2))
add("Overlap", paste("Sites on", a2, "only"), sum(on2 & !on1))
add("Overlap", "Sites on 450K and both arrays",
    sum(both & sites[[paste0("on_450k_", a1)]] & sites[[paste0("on_450k_", a2)]]),
    "for comparison with 450K-based resources (e.g. GoDMC, older clocks)")
add("Overlap", "Shared sites with the same hg38 position", sum(same_position, na.rm = TRUE),
    paste("of", sum(positioned), "shared sites with a position on both"))
add("Overlap", "Shared sites with a different hg38 position", sum(positioned & !same_position, na.rm = TRUE),
    "check before combining these")
add("Overlap", "Shared sites with a different probe design", sum(both & sites[[paste0("design_", a1)]] != sites[[paste0("design_", a2)]]),
    "Infinium I vs II; design affects the measured value range")

# The EPICv2 manifest also maps probes to EPICv1 names (EPICv1_Loci); a few
# differ from the probe's own name
mapped <- manifests[[a2]] %>% filter(!is.na(epicv1_loci), epicv1_loci != "")
if (nrow(mapped)) {
  add("Overlap", paste(a2, "probes with an EPICv1_Loci name"), nrow(mapped))
  add("Overlap", paste(a2, "probes whose EPICv1_Loci differs from their site name"),
      sum(mapped$epicv1_loci != mapped$site),
      "name matching (used by the pipeline) misses these")
}

# --- 4. After QC and in the mQTL scan ----------------------------------------

# Sites passing QC per data source: row names of MVALS_QC
for (data_source in names(sources)) {
  cat("Reading the row names of", sources[[data_source]]$mvals, "...\n")
  passed <- unique(to_site(rownames(readRDS(sources[[data_source]]$mvals))))
  invisible(gc())
  sites[[paste0("qc_", data_source)]] <- sites$site %in% passed
}
qc_cols <- paste0("qc_", names(sources))
qc_both <- reduce(map(qc_cols, ~ sites[[.x]]), `&`)
add("After QC", "Sites passing QC in every data source", sum(qc_both))
for (data_source in names(sources)) {
  array <- sources[[data_source]]$array
  other <- setdiff(arrays, array)
  passed <- sites[[paste0("qc_", data_source)]]
  add("After QC", paste(data_source, "sites passing QC"), sum(passed))
  add("After QC", paste(data_source, "sites passing QC, failing in another data source"),
      sum(passed & !qc_both & sites[[paste0("on_", other)]]), paste("on", other, "but removed by its QC"))
  add("After QC", paste(data_source, "sites passing QC, not on", other), sum(passed & !sites[[paste0("on_", other)]]))
}

# Sites the mQTL scan tested (at least one cis SNP), per data source and ancestry
mqtl_cols <- character(0)
for (data_source in names(sources)) {
  for (ancestry in names(sources[[data_source]]$mqtl)) {
    tested <- data.table::fread(sources[[data_source]]$mqtl[[ancestry]], select = c("cpg", "n_snps"))
    col <- paste("mqtl", data_source, ancestry, sep = "_")
    sites[[col]] <- sites$site %in% to_site(tested$cpg[tested$n_snps > 0])
    mqtl_cols <- c(mqtl_cols, col)
    add("mQTL", paste(data_source, ancestry, "sites tested"), sum(sites[[col]]))
  }
}
if (length(mqtl_cols) > 1) {
  add("mQTL", paste("Sites tested in every mQTL run:", paste(sub("^mqtl_", "", mqtl_cols), collapse = ", ")),
      sum(reduce(map(mqtl_cols, ~ sites[[.x]]), `&`)), "the CpGs a cross-ancestry comparison can use")
}

summary <- bind_rows(summary_rows)
write.csv(summary, OVERVIEW_ARRAY_SUMMARY, row.names = FALSE)
cat("\n")
print(as.data.frame(summary), row.names = FALSE, right = FALSE)
cat("\nSummary saved to:", OVERVIEW_ARRAY_SUMMARY, "\n")

# --- 5. Genomic context ------------------------------------------------------

# Context from the array that has the site (EPICv2 for shared sites: its
# annotation is native hg38)
context <- sites %>%
  mutate(island      = coalesce(.data[[paste0("island_", a2)]], .data[[paste0("island_", a1)]]),
         gene_region = coalesce(.data[[paste0("gene_region_", a2)]], .data[[paste0("gene_region_", a1)]]))
context_counts <- bind_rows(
  context %>% count(category, level = island) %>% mutate(feature = "island"),
  context %>% count(category, level = gene_region) %>% mutate(feature = "gene_region")) %>%
  mutate(level = replace_na(level, "No annotation")) %>%
  group_by(feature, category) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup() %>%
  select(feature, category, level, n, pct)
write.csv(context_counts, OVERVIEW_ARRAY_CONTEXT, row.names = FALSE)
cat("Context saved to:", OVERVIEW_ARRAY_CONTEXT, "\n")

data.table::fwrite(sites %>% select(-category, everything(), category), OVERVIEW_ARRAY_SITES)
cat("Per-site table saved to:", OVERVIEW_ARRAY_SITES, "\n")

# --- 6. Figures --------------------------------------------------------------

categories <- c("Both arrays", paste(a1, "only"), paste(a2, "only"))
overview_theme <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        strip.background = element_blank(), strip.text = element_text(face = "bold"))
save_fig <- function(plot, name) {
  ggsave(fig(name), plot, width = FIG_WIDTH, height = FIG_HEIGHT, dpi = FIG_RES)
  cat("Saved:", fig(name), "\n")
}

# Sites in each category at each stage; after QC and in the mQTL scan, "only"
# means kept in that array's data sources but not the other's
kept_by_array <- function(cols) {
  map(set_names(arrays), function(array) {
    own <- cols[map_chr(sub("^(qc|mqtl)_", "", cols), ~ sources[[names(sources)[startsWith(.x, names(sources))][1]]]$array) == array]
    if (length(own) == 0) return(NULL)
    reduce(map(own, ~ sites[[.x]]), `|`)
  })
}
stage_counts <- function(stage, keep) {
  if (any(map_lgl(keep, is.null))) return(NULL)
  tibble(stage = stage, category = categories,
         n = c(sum(keep[[a1]] & keep[[a2]]), sum(keep[[a1]] & !keep[[a2]]), sum(keep[[a2]] & !keep[[a1]])))
}
overlap <- bind_rows(
  stage_counts("Array manifests", list(on1, on2) %>% set_names(arrays)),
  stage_counts("Passed QC", kept_by_array(qc_cols)),
  stage_counts("Tested in mQTL scan", kept_by_array(mqtl_cols))) %>%
  mutate(stage = factor(stage, levels = unique(stage)),
         category = factor(category, levels = rev(categories)))

overlap_plot <- ggplot(overlap, aes(n / 1000, category)) +
  geom_col(fill = OVERVIEW_BAR_COLOR, width = 0.7) +
  geom_text(aes(label = scales::comma(n)), hjust = -0.1, size = 2.6) +
  facet_grid(rows = vars(stage)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "CpG sites (thousands)", y = NULL) +
  overview_theme
save_fig(overlap_plot, "05_array_overlap")

# Shares of all sites in the category (sites without an annotation are left
# out of the bars but counted in n)
island_shares <- context_counts %>%
  filter(feature == "island") %>%
  group_by(category) %>%
  mutate(share = n / sum(n), total = sum(n)) %>%
  ungroup() %>%
  filter(level %in% names(OVERVIEW_ISLAND_COLORS)) %>%
  mutate(level = factor(level, levels = names(OVERVIEW_ISLAND_COLORS)),
         category = factor(paste0(category, " (n=", scales::comma(total), ")"),
                           levels = paste0(categories, " (n=", scales::comma(total[match(categories, category)]), ")")) %>%
           fct_rev())
context_plot <- ggplot(island_shares, aes(share, category, fill = level)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.4, position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = ifelse(share >= OVERVIEW_MIN_LABEL_SHARE, scales::percent(share, 1), "")),
            position = position_stack(vjust = 0.5, reverse = TRUE), size = 2.6) +
  scale_fill_manual(values = OVERVIEW_ISLAND_COLORS) +
  scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.03))) +
  labs(x = "Share of sites (CpG island context)", y = NULL, fill = NULL) +
  overview_theme +
  theme(legend.position = "bottom")
save_fig(context_plot, "06_array_context")

cat("\nNext: see overview/README.md for the outputs\n")
