# =============================================================================
# IDAT and sample QC helpers
# Author: GP2 Subtypes and Mechanisms - M.E., M.P.
# Date: Sept 30, 2026
# Description: QC steps used by the 01 and 02 scripts that work on idat files
#              or the raw minfi objects: checking idat files exist, SeSAMe
#              per-sample QC stats, sex check, and low bead count probes
# =============================================================================

# Adds both_exist (TRUE if both <Basename>_Red.idat and _Grn.idat exist) to
# the sample sheet and reports samples with missing files, showing id_cols
check_idats <- function(sample_sheet, id_cols = "Basename") {
  cat("\nVerifying idat files exist on disk...\n")
  sample_sheet$both_exist <- file.exists(paste0(sample_sheet$Basename, "_Red.idat")) &
                             file.exists(paste0(sample_sheet$Basename, "_Grn.idat"))

  cat("Samples with both Red and Green idat files:",
      sum(sample_sheet$both_exist), "/", nrow(sample_sheet), "\n")
  if (any(!sample_sheet$both_exist)) {
    cat("WARNING:", sum(!sample_sheet$both_exist), "samples missing idat files:\n")
    print(sample_sheet[!sample_sheet$both_exist, id_cols, drop = FALSE])
  } else {
    cat("All idat files found!\n")
  }
  sample_sheet
}

# SeSAMe per-sample QC stats (detection pOOBAH, intensity, dye bias, beta
# distribution), one row per basename in the same order. Stops if any sample
# fails, naming the samples.
# NOTE: sesame's bisulfite conversion score (bisConversionControl) fails on
# EPICv2 in sesame 1.24, so it is not included
sesame_qc_stats <- function(basenames, cores = max(1, parallel::detectCores() - 1)) {
  # The idats are on a gcsfuse mount, where reads occasionally fail with
  # "error reading from connection" under parallel load; retry those reads
  read_idat_pair <- function(b, attempts = 3) {
    for (i in seq_len(attempts)) {
      sdf <- try(sesame::readIDATpair(b), silent = TRUE)
      if (!inherits(sdf, "try-error")) return(sdf)
      Sys.sleep(5 * i)
    }
    stop(sprintf("reading %s failed after %d attempts: %s", b, attempts, sdf))
  }

  # mc.preschedule = FALSE runs each sample as its own job, so one error
  # only affects that sample instead of every sample on the same core
  stats <- parallel::mclapply(basenames, function(b) {
    sdf <- read_idat_pair(b)
    as.data.frame(sesame::sesameQC_getStats(sesame::sesameQC_calcStats(sdf)))
  }, mc.cores = cores, mc.preschedule = FALSE)

  # mclapply returns errors instead of stopping; report which samples failed
  failed <- vapply(stats, inherits, logical(1), "try-error")
  if (any(failed)) {
    stop("SeSAMe QC failed for: ", paste(basenames[failed], collapse = ", "),
         "\n", paste(unique(unlist(stats[failed])), collapse = "\n"))
  }
  dplyr::bind_rows(stats)
}

# Predicted vs reported sex per sample, in the column order of mSet (a raw
# MethylSet). reported_sex is "Female"/"Male" (R12); anything else is
# "Unknown" and never counts as discordant. xMed/yMed are the median X and Y
# chromosome intensities getSex uses (plotted in qc_03)
check_sex <- function(mSet, reported_sex, cutoff = -2) {
  predicted <- minfi::getSex(minfi::mapToGenome(mSet), cutoff = cutoff)
  reported_sex_label <- dplyr::case_when(reported_sex == "Female" ~ "F",
                                         reported_sex == "Male"   ~ "M",
                                         TRUE                     ~ "Unknown")
  data.frame(reported_sex       = reported_sex,
             predicted_sex      = predicted$predictedSex,
             reported_sex_label = reported_sex_label,
             sex_discordant     = reported_sex_label != predicted$predictedSex &
                                  reported_sex_label != "Unknown",
             xMed               = predicted$xMed,
             yMed               = predicted$yMed,
             row.names = NULL)
}

# Names of probes with fewer than min_beads beads in more than max_frac of
# the samples in keep. rgSet must be read with extended = TRUE; ann is
# getAnnotation(rgSet).
# getNBeads() rows are bead addresses, not probe names; map them to probes.
# Type I probes use two addresses (A and B), so take the lower count of the two
find_low_bead_probes <- function(rgSet, ann, keep, min_beads, max_frac = 0.05) {
  nbeads   <- minfi::getNBeads(rgSet)[, keep, drop = FALSE]
  nbeads_A <- nbeads[match(as.character(ann$AddressA), rownames(nbeads)), , drop = FALSE]
  nbeads_B <- nbeads[match(as.character(ann$AddressB), rownames(nbeads)), , drop = FALSE]
  nbeads_probe <- ifelse(is.na(nbeads_B), nbeads_A, pmin(nbeads_A, nbeads_B))
  low <- rowSums(nbeads_probe < min_beads, na.rm = TRUE) > (max_frac * ncol(nbeads_probe))
  ann$Name[low]
}
