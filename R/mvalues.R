# =============================================================================
# M value helpers
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Sept 30, 2026
# Description: Helpers for beta and M value matrices, used by qc/2_qc.R and
#              the preprocessing scripts
# =============================================================================

# Sites measured by more than one probe (EPICv2 replicates, e.g.
# cg00002033_TC11 and cg00002033_TC12: probe names sharing the part before
# "_") averaged to one row per site, ignoring NAs. The averaged row keeps the
# name of the site's first probe (sorted), so lookups by probe name (array
# annotation, CPG_ANNOTATION, sesame manifest) still work; rows stay in their
# original order. EPICv1/450k probe names have no suffix, so they're left
# unchanged. Returns list(beta, map): map has one row per probe of a
# replicated site (site, probe, kept = the row name it was averaged into)
collapse_replicate_probes <- function(beta) {
  site      <- sub("_.*$", "", rownames(beta))
  is_rep    <- site %in% site[duplicated(site)]
  rep_probe <- rownames(beta)[is_rep]
  rep_site  <- site[is_rep]
  kept      <- tapply(rep_probe, rep_site, function(p) sort(p)[1])
  map <- data.frame(site = rep_site, probe = rep_probe,
                    kept = unname(kept[rep_site]), row.names = NULL)
  if (!any(is_rep)) return(list(beta = beta, map = map))

  x   <- beta[is_rep, , drop = FALSE]
  avg <- rowsum(x, rep_site, na.rm = TRUE) / rowsum(1 * !is.na(x), rep_site)
  avg[is.nan(avg)] <- NA                       # every probe of the site NA
  beta[match(kept, rownames(beta)), ] <- avg[names(kept), , drop = FALSE]
  list(beta = beta[!is_rep | rownames(beta) %in% kept, , drop = FALSE], map = map)
}

# Betas of exactly 0 or 1 (which Funnorm can produce) give M values of -Inf or
# Inf, which break ComBat, limma and other linear models. Replace them with the
# M value of beta = beta_bound or 1 - beta_bound (0.001 -> about -/+9.97).
# Finite values and NAs are left unchanged
cap_infinite_m <- function(m, beta_bound = 0.001) {
  lim <- log2((1 - beta_bound) / beta_bound)
  n_inf <- sum(is.infinite(m))
  if (n_inf > 0) {
    cat("Capping", n_inf, "infinite M values (beta of exactly 0 or 1) at +/-",
        round(lim, 2), "\n")
    m[m == Inf]  <- lim
    m[m == -Inf] <- -lim
  }
  m
}
