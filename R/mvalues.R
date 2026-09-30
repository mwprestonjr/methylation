# =============================================================================
# M value helpers
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Sept 30, 2026
# Description: Helpers for M value matrices, used by Script 02 and the
#              preprocessing scripts
# =============================================================================

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
