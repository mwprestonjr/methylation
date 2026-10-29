# =============================================================================
# Clock estimation helpers
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 29, 2026
# Description: Probe-name collapsing and clock CpG coverage, shared by
#              clocks/1_estimate.R and clocks/5_epicv2_simulation.R
# =============================================================================

# EPICv2 has replicate probes for some CpGs; average them into one row
# (e.g. cg00000029_TC21 -> cg00000029). EPICv1 names are left as they are
collapse_epicv2 <- function(b) {
  cg_id <- sub("_.*$", "", rownames(b))
  keep  <- startsWith(cg_id, "cg") | startsWith(cg_id, "ch.")   # clocks use cg and ch. probes
  b     <- b[keep, , drop = FALSE]
  cg_id <- cg_id[keep]

  dup    <- cg_id %in% cg_id[duplicated(cg_id)]
  single <- b[!dup, , drop = FALSE]
  rownames(single) <- cg_id[!dup]

  sums  <- rowsum(b[dup, , drop = FALSE], cg_id[dup], na.rm = TRUE)
  n     <- rowsum((!is.na(b[dup, , drop = FALSE])) * 1, cg_id[dup])
  multi <- sums / n

  rbind(single, multi)
}

# CpGs per clock (methylclockData coefficient tables; first row is the
# intercept except for Hannum) and DunedinPACE (model probes, and the probes
# PACEProjector uses to normalize; it needs PACE_COV_THRESH of both)
clock_cpg_list <- function() {
  list(
    Horvath     = methylclockData::get_coefHorvath()$CpGmarker[-1],
    Hannum      = methylclockData::get_coefHannum()$CpGmarker,
    Levine      = methylclockData::get_coefLevine()$CpGmarker[-1],
    skinHorvath = methylclockData::get_coefSkin()$CpGmarker[-1],
    PedBE       = methylclockData::get_coefPedBE()$CpGmarker[-1],
    Wu          = methylclockData::get_coefWu()$CpGmarker[-1],
    TL          = methylclockData::get_coefTL()$CpGmarker[-1],
    BLUP        = methylclockData::get_coefBLUP()$CpGmarker[-1],
    EN          = methylclockData::get_coefEN()$CpGmarker[-1],
    DunedinPACE = DunedinPACE::getRequiredProbes()$DunedinPACE,
    DunedinPACE_normalization = DunedinPACE::getRequiredProbes(backgroundList = TRUE)$DunedinPACE
  )
}

# Share of each clock's CpGs among `cpgs` (the CpG IDs on the array).
# pct_on_array is what the coverage thresholds are applied to
clock_coverage_table <- function(clock_cpgs, cpgs) {
  imap_dfr(clock_cpgs, function(clock_set, clock) {
    n_on_array <- sum(clock_set %in% cpgs)
    tibble(clock        = clock,
           n_cpgs       = length(clock_set),
           n_on_array   = n_on_array,
           pct_on_array = 100 * n_on_array / length(clock_set))
  })
}

# Which clocks were estimated: below the coverage threshold a clock is all NA.
# The normalization probes row has no estimate of its own
clocks_estimated <- function(clock_cpgs, dnam_age, pace) {
  c(sapply(setdiff(names(clock_cpgs), c("DunedinPACE", "DunedinPACE_normalization")),
           function(clock) clock %in% names(dnam_age) && !all(is.na(dnam_age[[clock]]))),
    DunedinPACE = !all(is.na(pace$DunedinPACE)))
}
