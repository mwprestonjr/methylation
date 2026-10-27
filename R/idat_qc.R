# =============================================================================
# IDAT and sample QC helpers
# Author: GP2 Subtypes and Mechanisms - M.E., M.P.
# Date: Sept 30, 2026
# Description: QC steps used by qc/1 and qc/2 that work on idat files
#              or the raw minfi objects: checking idat files exist, SeSAMe
#              per-sample QC stats, sex check, duplicate arrays of one person,
#              and low bead count probes
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

# Duplicate arrays of one person among the samples in keep, and which one to
# keep. person identifies the participant per rgSet column (NA or "" = can't
# be matched, never a duplicate); quality ranks the arrays (higher is
# better), ties going to the higher sentrix_id (more recent chip). Sample
# swap check: if any two arrays of a person have rs SNP probe betas
# (getSnpBeta) correlating below min_cor, or a different predicted_sex, all
# of that person's arrays are removed. Returns one row per array of a
# duplicated person: index (rgSet column), person, n_arrays, min_snp_cor
# (lowest pairwise correlation), same_sex, and decision ("kept", "removed:
# replicate" or "removed: possible swap")
find_duplicate_arrays <- function(rgSet, person, keep, quality, sentrix_id,
                                  predicted_sex, min_cor) {
  p   <- ifelse(keep & !is.na(person) & person != "", person, NA)
  idx <- which(p %in% unique(p[duplicated(p) & !is.na(p)]))
  if (length(idx) == 0) {
    return(data.frame(index = integer(0), person = character(0), n_arrays = integer(0),
                      min_snp_cor = numeric(0), same_sex = logical(0),
                      decision = character(0)))
  }
  snp <- minfi::getSnpBeta(rgSet[, idx])

  groups <- split(seq_along(idx), p[idx])
  do.call(rbind, lapply(groups, function(j) {
    i     <- idx[j]
    r     <- suppressWarnings(cor(snp[, j], use = "pairwise.complete.obs"))
    min_r <- min(r[upper.tri(r)])
    same_sex <- length(unique(predicted_sex[i])) == 1
    swap  <- !isTRUE(min_r >= min_cor) || !same_sex   # NA correlation counts as a swap
    best  <- i[order(-quality[i], -as.numeric(sentrix_id[i]))[1]]
    data.frame(index = i, person = p[i], n_arrays = length(i),
               min_snp_cor = min_r, same_sex = same_sex,
               decision = if (swap) "removed: possible swap"
                          else ifelse(i == best, "kept", "removed: replicate"),
               row.names = NULL)
  }))
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
