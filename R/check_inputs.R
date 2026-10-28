# =============================================================================
# Step 0: Check input paths
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 28, 2026
# Description: Checks that the input files and folders a module reads (set in
#              config.R and mqtl/config.R) exist and are readable, that the
#              output folder is writable and, for mqtl, that plink2 is on the
#              PATH. Empty folders are reported too: a gcsfuse bucket that
#              isn't mounted shows up as an empty folder. Stops with an error
#              if anything is missing, so make stops before the module runs.
# Usage:       Rscript R/check_inputs.R <data source> [modules] [ancestry]
#              modules: comma-separated, default all (qc, preprocessing,
#              clocks, ewas, mqtl); ancestry: for mqtl, checks that
#              ancestry's genotype files, e.g.
#                Rscript R/check_inputs.R psomagen qc
#                Rscript R/check_inputs.R ppmi_p140 mqtl EUR
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

source("config.R")
source("mqtl/config.R")

args     <- commandArgs(trailingOnly = TRUE)
MODULES  <- c("qc", "preprocessing", "clocks", "ewas", "mqtl")
modules  <- if (length(args) >= 2 && !args[2] %in% c("", "all")) unlist(strsplit(args[2], ",")) else MODULES
ancestry <- if (length(args) >= 3 && args[3] != "") args[3] else NA
if (!all(modules %in% MODULES)) {
  stop("Unknown module: ", paste(setdiff(modules, MODULES), collapse = ", "),
       " - one of: ", paste(MODULES, collapse = ", "))
}
cat("Modules:", paste(modules, collapse = ", "), "\n")

# --- 1. Inputs, and the modules that read them -------------------------------

# One row per input: config name, path, file or folder, modules reading it
inputs <- list(
  list("GP2_MASTER_KEY", GP2_MASTER_KEY, "file", c("qc", "clocks", "mqtl")),
  list("CLINICAL_DATA",  CLINICAL_DATA,  "file", c("clocks", "ewas")),
  list("GENCODE_GTF",    GENCODE_GTF,    "file", "preprocessing")
)

if (DATA_SOURCE == "ppmi_p140") {
  inputs <- c(inputs, list(
    list("IDAT_DIR",           IDAT_DIR,           "dir",  "qc"),
    list("LINK_LIST",          LINK_LIST,          "file", "qc"),
    list("PARTICIPANT_STATUS", PARTICIPANT_STATUS, "file", "qc"),
    list("AGE_AT_VISIT",       AGE_AT_VISIT,       "file", "qc"),
    # EPICv1 annotation is hg19, lifted to hg38
    list("CHAIN_HG19_HG38",    CHAIN_HG19_HG38,    "file", c("preprocessing", "mqtl"))
  ))
}

if (DATA_SOURCE == "psomagen") {
  inputs <- c(inputs, list(
    list("GP2_INTERNAL_KEY", GP2_INTERNAL_KEY, "file", "qc")
  ))
  for (ds in DATASETS) {
    inputs <- c(inputs, list(
      list(paste0("DATASETS: ", ds),
           file.path(DIR_DELIVERY, ds, paste0(ds, "_QC_table.csv")), "file", "qc")
    ))
  }
}

# Genotypes: one file set per ancestry (GENO_PFILE_NAME); without an ancestry
# only the folder is checked
if (is.na(ancestry)) {
  inputs <- c(inputs, list(list("GENO_PFILE_PATH", GENO_PFILE_PATH, "dir", "mqtl")))
} else {
  pfile <- file.path(GENO_PFILE_PATH, gsub("{ANCESTRY}", ancestry, GENO_PFILE_NAME, fixed = TRUE))
  for (ext in c(".pgen", ".pvar", ".psam")) {
    inputs <- c(inputs, list(list(paste0("GENO_PFILE (", ancestry, ")"),
                                  paste0(pfile, ext), "file", "mqtl")))
  }
}

# --- 2. Check ----------------------------------------------------------------

check_path <- function(path, type) {
  path <- path.expand(path)
  if (type == "dir") {
    if (!dir.exists(path))                  return("missing")
    if (file.access(path, 4) != 0)          return("not readable")
    if (length(list.files(path)) == 0)      return("empty (bucket not mounted?)")
  } else {
    if (!file.exists(path))                 return("missing")
    if (dir.exists(path))                   return("a folder, expected a file")
    if (file.access(path, 4) != 0)          return("not readable")
    if (file.size(path) == 0)               return("empty")
  }
  "ok"
}

needed <- Filter(function(x) any(x[[4]] %in% modules), inputs)
status <- data.frame(
  name   = vapply(needed, `[[`, character(1), 1),
  path   = vapply(needed, `[[`, character(1), 2),
  status = vapply(needed, function(x) check_path(x[[2]], x[[3]]), character(1))
)

# Output folder (created by config.R) must be writable
status <- rbind(status, data.frame(
  name = "DIR_OUTPUT", path = DIR_OUTPUT,
  status = if (file.access(DIR_OUTPUT, 2) == 0) "ok" else "not writable"))

if ("mqtl" %in% modules) {
  status <- rbind(status, data.frame(
    name = "plink2", path = Sys.which("plink2"),
    status = if (nzchar(Sys.which("plink2"))) "ok" else "not on PATH"))
}

cat("\n")
for (i in seq_len(nrow(status))) {
  cat(sprintf("  %-4s %-28s %s\n", if (status$status[i] == "ok") "ok" else "FAIL",
              status$name[i], status$path[i]))
  if (status$status[i] != "ok") cat("       ->", status$status[i], "\n")
}

# Missing files usually mean a bucket isn't mounted: name the mount points
# (~/gp2_release12, /mnt/<bucket>) that don't exist or are empty
failed <- status[status$status != "ok", ]
if (nrow(failed) > 0) {
  roots <- unique(path.expand(sub("^((~|/home/[^/]+)/[^/]+|/mnt/[^/]+).*", "\\1", failed$path)))
  roots <- roots[!dir.exists(roots) | vapply(roots, function(r) length(list.files(r)) == 0, logical(1))]
  if (length(roots) > 0) {
    cat("\nNot mounted (missing or empty):", paste(roots, collapse = ", "),
        "\n  e.g. the gp2_release12 bucket at ~/gp2_release12 (see README.md)\n")
  }
  stop(nrow(failed), " input(s) missing or unreadable (above)", call. = FALSE)
}
if ("mqtl" %in% modules && is.na(ancestry)) {
  cat("\nGenotype files not checked: give an ancestry (e.g. EUR) to check them\n")
}
cat("\nAll inputs found.\n")
