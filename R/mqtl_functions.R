# =============================================================================
# mQTL functions
# Author: GP2 Subtypes and Mechanisms - M.E., M.P.
# Date: Oct 7, 2026
# Description: Steps shared by mqtl/3_mqtl.R, mqtl/4_mqtl_permutations.R
#              and mqtl/tune_meth_pcs.R: loading and preparing the
#              methylation data and covariates, latent methylation PCs, the
#              per-chromosome cis-mQTL scan, and the per-CpG permutations
#              and their significance.
#              Needs config.R, mqtl/config.R, R/mqtl_setup.R,
#              R/array_profiles.R and R/cpg_annotation.R sourced first, and minfi, tidyverse,
#              data.table, matrixStats, MatrixEQTL, GenomicRanges, rtracklayer
#              loaded
# =============================================================================

# Methylation data and known covariates for the samples with genotypes:
#   meth_int    inverse-normal transformed ComBat M values (CpGs x samples)
#   cpg_pos     hg38 CpG positions (MatrixEQTL genepos format)
#   known_covs  covariates: age, sex, cell proportions, genotype PCs,
#               phenotype and chip row indicators
#   sample_map  samples in the final set, in the column order of meth_int
#   geno_dir    folder with the step 2 genotype files
#   probe_regions  hg38 range of SNP positions in each CpG's probe
#               (probe_snp_regions), excluded from the scan
mqtl_prepare <- function() {
  geno_dir <- file.path(MQTL_DIR, "genotypes")

  # --- Inputs ---
  cat("Loading inputs...\n")
  combat_mVals <- readRDS(COMBAT_MVALS)
  targets      <- read.csv(SAMPLE_SHEET_FINAL,
                           colClasses = c(GP2ID       = "character",
                                          clinical_id = "character",
                                          Sentrix_ID  = "character",
                                          Basename    = "character")) %>%
    mutate(meth_id = basename(Basename))
  sample_map   <- read.csv(MQTL_SAMPLE_MAP, colClasses = "character")

  # Genotype PCs; plink2 writes either "#IID" or "#FID IID" as the header
  geno_pcs <- fread(file.path(geno_dir, "geno_pcs.eigenvec"))
  setnames(geno_pcs, sub("^#", "", names(geno_pcs)))
  geno_pcs <- geno_pcs %>%
    select(IID, all_of(paste0("PC", 1:MQTL_N_GENO_PCS))) %>%
    rename_with(~ paste0("geno_", .x), starts_with("PC"))

  # Final sample set: in methylation, sample map, and genotypes after plink QC,
  # with age and sex available (model.matrix would silently drop incomplete rows)
  complete_covs <- targets$meth_id[!is.na(targets$age) & targets$sex %in% c("Female", "Male")]
  sample_map <- sample_map %>%
    filter(GP2ID %in% geno_pcs$IID,           # genotype IIDs are GP2IDs
           meth_id %in% colnames(combat_mVals),
           meth_id %in% complete_covs)
  cat("Samples in final mQTL set:", nrow(sample_map), "\n")

  # Everything below is ordered by sample_map$meth_id
  meth <- combat_mVals[, sample_map$meth_id]
  rm(combat_mVals)
  covs <- sample_map %>%
    select(meth_id) %>%
    left_join(targets, by = "meth_id") %>%
    left_join(geno_pcs, by = c("GP2ID" = "IID"))
  stopifnot(all(covs$meth_id == colnames(meth)))

  # --- CpG positions (hg38) ---
  # From the array's annotation, lifted over for hg19 arrays (cpg_positions_hg38
  # in R/cpg_annotation.R); one array type per data source
  array_type <- unique(targets$Array)
  stopifnot(length(array_type) == 1, array_type %in% names(ARRAY_PROFILES))
  profile <- ARRAY_PROFILES[[array_type]]
  cat("\nArray:", array_type, "- annotation:", profile$anno_pkg, "(", profile$genome, ")\n")

  pos     <- cpg_positions_hg38(rownames(meth), profile, CHAIN_HG19_HG38)
  cpg_pos <- data.frame(geneid = pos$cpg, chr = pos$chr, left = pos$pos, right = pos$pos)
  meth <- meth[cpg_pos$geneid, ]

  # --- Methylation transform ---
  cat("\nInverse-normal transforming M values per CpG...\n")
  int_transform <- function(x) qnorm((rank(x, ties.method = "average") - 0.5) / length(x))
  meth_int <- t(apply(meth, 1, int_transform))
  dimnames(meth_int) <- dimnames(meth)
  rm(meth)

  # --- Known covariates ---
  # Cell proportions without Neu (they sum to ~1). Phenotype is coded as
  # indicator columns against the most common phenotype, so it works for PD
  # vs Control and for more groups (PD, DLB, MSA, ...); it is left out if all
  # samples share one phenotype
  known_covs <- covs %>%
    transmute(age  = as.numeric(age),
              male = as.integer(sex == "Male"),
              CD8T, CD4T, NK, Bcell, Mono,
              across(starts_with("geno_PC")))
  phenotypes <- names(sort(table(covs$phenotype), decreasing = TRUE))
  if (length(phenotypes) > 1) {
    for (ph in phenotypes[-1]) {
      known_covs[[paste0("phenotype_", make.names(ph))]] <- as.integer(covs$phenotype == ph)
    }
  }
  cat("Phenotype reference group:", phenotypes[1], "\n")

  # Position on the chip (row R01-R08 of Sentrix_Position): ComBat corrects for
  # chip or plate but not for the row within a chip, which remains one of the
  # main technical effects (see preprocessing/2_variation_sources.R). Coded as
  # indicator columns against the first row
  chip_row <- substr(covs$Sentrix_Position, 1, 3)
  rows     <- sort(unique(chip_row))
  for (r in rows[-1]) {
    known_covs[[paste0("chip_row_", r)]] <- as.integer(chip_row == r)
  }
  cat("Chip rows:", paste(rows, collapse = ", "), "(reference", rows[1], ")\n")

  # Sites averaged from EPICv2 replicate probes (qc/2_qc.R); missing for QC
  # outputs from before the averaging, which have no averaged sites
  replicate_map <- if (file.exists(REPLICATE_PROBES)) read.csv(REPLICATE_PROBES) else NULL

  list(meth_int = meth_int, cpg_pos = cpg_pos, known_covs = known_covs,
       sample_map = sample_map, geno_dir = geno_dir,
       probe_regions = probe_snp_regions(cpg_pos, profile, replicate_map))
}

# hg38 range of SNP positions [lo, hi] that fall in each CpG's probe, where a
# SNP can change probe binding or the extension base and mimic an mQTL.
# Fixed by the Infinium design; positions counted from the CpG's C towards the
# probe's 5' end:
#   -1      the CpG's G (extension base of Infinium I probes)
#    0      the CpG's C (extension base of Infinium II probes)
#    0-49   Infinium I probe body;  1-50  Infinium II probe body
# The probe strand and CpG position come from the sesame hg38 manifest
# (probe_strand of the Zhou lab InfiniumAnnotation, aligned to hg38): a "+"
# probe lies at lower coordinates than its CpG, a "-" probe at higher ones.
# CpGs missing from the manifest, or on another chromosome there, get the
# range on both sides of their position in cpg_pos. So do sites averaged from
# replicate probes on different strands (replicate_map, from
# REPLICATE_PROBES: site, probe, kept = the row name of the average).
# Returns a data.frame: geneid, lo, hi, in cpg_pos order
probe_snp_regions <- function(cpg_pos, profile, replicate_map = NULL) {
  from <- -1
  to   <- 50
  man  <- sesameData::sesameData_getManifestGRanges(profile$sesame_platform, genome = "hg38")
  m    <- match(cpg_pos$geneid, names(man))
  # GRanges start is the CpG's first base: the C of a "+" strand CpG
  beg    <- start(man)[m]
  strand <- as.character(strand(man))[m]
  known  <- !is.na(m) & strand %in% c("+", "-") &
            sub("^chr", "", as.character(seqnames(man))[m]) == cpg_pos$chr
  if (!is.null(replicate_map) && nrow(replicate_map) > 0) {
    rep_strand <- as.character(strand(man))[match(replicate_map$probe, names(man))]
    n_strands  <- tapply(rep_strand, replicate_map$kept,
                         function(s) length(unique(s[s %in% c("+", "-")])))
    mixed <- cpg_pos$geneid %in% names(n_strands)[n_strands > 1]
    cat("Sites averaged from probes on both strands:", sum(mixed), "\n")
    known <- known & !mixed
  }
  cat("Probe strand from the sesame", profile$sesame_platform, "manifest:", sum(known), "of",
      nrow(cpg_pos), "CpGs; the rest are excluded on both sides\n")

  lo <- ifelse(strand %in% "+", beg - to, beg + 1 + from)
  hi <- ifelse(strand %in% "+", beg - from, beg + 1 + to)
  lo[!known] <- cpg_pos$left[!known] - to - 1
  hi[!known] <- cpg_pos$left[!known] + to + 1
  data.frame(geneid = cpg_pos$geneid, lo = lo, hi = hi)
}

# Latent methylation PCs: the first n_pcs principal components of the
# methylation residuals after the known covariates (top 50,000 most variable
# CpGs). They soak up unmeasured technical/biological variation and boost
# power. Components are nested: the first k of n_pcs components are the same
# as computing k components directly
latent_meth_pcs <- function(meth_int, known_covs, n_pcs) {
  if (n_pcs == 0) return(matrix(numeric(0), nrow = ncol(meth_int), ncol = 0))
  X        <- model.matrix(~ ., data = known_covs)
  resid    <- meth_int - t(X %*% solve(crossprod(X), crossprod(X, t(meth_int))))
  top_var  <- order(rowVars(resid), decreasing = TRUE)[seq_len(min(50000, nrow(resid)))]
  meth_pcs <- prcomp(t(resid[top_var, ]), center = TRUE, scale. = FALSE,
                     rank. = n_pcs)$x
  colnames(meth_pcs) <- paste0("meth_PC", seq_len(ncol(meth_pcs)))
  meth_pcs
}

# Covariates (a data frame, one row per sample) as MatrixEQTL SlicedData
covariates_sliced <- function(all_covs, sample_names) {
  cvrt <- SlicedData$new()
  cvrt$CreateFromMatrix(t(as.matrix(all_covs)))
  colnames(cvrt) <- sample_names
  cvrt
}

# Genotypes of one chromosome (step 2 .traw) for the samples, in the column
# order of meth_int: list(snps = SlicedData, snp_info = data.table), or with
# sliced = FALSE list(geno = matrix, SNPs x samples, snp_info)
load_chr_genotypes <- function(geno_dir, chr, sample_map, sample_names, sliced = TRUE) {
  # .traw sample columns are "<FID>_<IID>"
  traw     <- fread(file.path(geno_dir, paste0("geno_chr", chr, ".traw")))
  geno_ids <- names(traw)[-(1:6)]
  col_idx  <- sapply(sample_map$GP2ID,
                     function(id) which(endsWith(geno_ids, paste0("_", id))))
  stopifnot(is.numeric(col_idx), length(col_idx) == nrow(sample_map))

  snp_info <- traw[, .(snpid = SNP, chr = as.character(CHR), pos = POS,
                       counted = COUNTED, alt = ALT)]
  geno_mat <- as.matrix(traw[, geno_ids[col_idx], with = FALSE])
  dimnames(geno_mat) <- list(snp_info$snpid, sample_names)
  if (!sliced) return(list(geno = geno_mat, snp_info = snp_info))

  snps <- SlicedData$new()
  snps$CreateFromMatrix(geno_mat)
  snps$ResliceCombined(sliceSize = 5000)
  list(snps = snps, snp_info = snp_info)
}

# SNPs in each CpG's cis window and probe, as index ranges into snp_bp (the
# chromosome's SNP positions, sorted). MatrixEQTL tests the SNPs at
# [left - MQTL_CIS_WINDOW, right + MQTL_CIS_WINDOW], bounds included: first
# to last. The SNPs in the probe region (probe_snp_regions; all within the
# window) are probe_first to probe_last. A range is empty when last < first.
# genes and region are the chromosome's rows of cpg_pos and probe_regions.
# Used by run_cis_chr and cis_perm_chr, so both test the same SNPs
cis_snp_ranges <- function(genes, region, snp_bp) {
  data.frame(first       = findInterval(genes$left - MQTL_CIS_WINDOW - 1, snp_bp) + 1,
             last        = findInterval(genes$right + MQTL_CIS_WINDOW, snp_bp),
             probe_first = findInterval(region$lo - 1, snp_bp) + 1,
             probe_last  = findInterval(region$hi, snp_bp))
}

# cis-mQTL scan of one chromosome: tests within +/- MQTL_CIS_WINDOW, keeps
# pairs with p < MQTL_P_CIS_SAVE. Pairs with the SNP in the CpG's probe
# (probe_regions, from probe_snp_regions) are dropped and not counted as
# tests. Returns list(eqtls = data frame, ntests, n_probe_tests = tests
# excluded, n_probe_saved = saved pairs excluded, cpg_tests = one row per
# CpG on the chromosome: cpg, chr, n_snps tested, n_probe_snps excluded,
# probe_lo/probe_hi the excluded range)
run_cis_chr <- function(geno, meth_int, cpg_pos, probe_regions, chr, cvrt) {
  chr_cpgs <- cpg_pos$chr == as.character(chr)
  gene <- SlicedData$new()
  gene$CreateFromMatrix(meth_int[chr_cpgs, , drop = FALSE])
  gene$ResliceCombined(sliceSize = 2000)

  me <- Matrix_eQTL_main(
    snps                  = geno$snps,
    gene                  = gene,
    cvrt                  = cvrt,
    output_file_name      = NULL,          # trans off
    pvOutputThreshold     = 0,
    output_file_name.cis  = NULL,
    pvOutputThreshold.cis = MQTL_P_CIS_SAVE,
    useModel              = modelLINEAR,
    errorCovariance       = numeric(),
    snpspos               = geno$snp_info[, .(snpid, chr, pos)],
    genepos               = cpg_pos[chr_cpgs, ],
    cisDist               = MQTL_CIS_WINDOW,
    pvalue.hist           = FALSE,
    min.pv.by.genesnp     = FALSE,
    noFDRsaveMemory       = FALSE,
    verbose               = FALSE
  )

  eqtls <- me$cis$eqtls %>%
    transmute(cpg = as.character(gene), snp = as.character(snps),
              beta, t_stat = statistic, p = pvalue) %>%
    left_join(as.data.frame(geno$snp_info)[, c("snpid", "pos", "counted", "alt")],
              by = c("snp" = "snpid")) %>%
    mutate(chr = chr)

  # SNPs tested per CpG (cis_snp_ranges), which add up to MatrixEQTL's tests
  genes      <- cpg_pos[chr_cpgs, ]
  region     <- probe_regions[chr_cpgs, ]
  win        <- cis_snp_ranges(genes, region, sort(geno$snp_info$pos))
  n_cis      <- win$last - win$first + 1
  stopifnot(sum(n_cis) == me$cis$ntests)
  n_in_probe <- win$probe_last - win$probe_first + 1
  cpg_tests  <- data.frame(cpg = genes$geneid, chr = chr,
                           n_snps = n_cis - n_in_probe, n_probe_snps = n_in_probe,
                           probe_lo = region$lo, probe_hi = region$hi)

  r <- match(eqtls$cpg, region$geneid)
  in_probe <- eqtls$pos >= region$lo[r] & eqtls$pos <= region$hi[r]
  list(eqtls = eqtls[!in_probe, ], ntests = sum(cpg_tests$n_snps),
       n_probe_tests = sum(n_in_probe), n_probe_saved = sum(in_probe),
       cpg_tests = cpg_tests)
}

# Combine the per-chromosome results: FDR over all tests, SNP-CpG distance,
# and the lead SNP per CpG. Saved rows are exactly the smallest p-values, so
# BH with n = all tests gives the same q-values as adjusting the full set
summarise_cis <- function(cis_results, n_tests, cpg_pos) {
  cis <- bind_rows(cis_results) %>%
    left_join(cpg_pos %>% select(cpg = geneid, cpg_pos = left), by = "cpg") %>%
    mutate(distance = pos - cpg_pos,
           fdr      = p.adjust(p, method = "BH", n = n_tests)) %>%
    arrange(p)
  lead <- cis %>%
    group_by(cpg) %>%
    slice_min(p, n = 1, with_ties = FALSE) %>%
    ungroup()
  list(cis = cis, lead = lead)
}

# The p-value threshold that FDR < fdr corresponds to: the largest p-value
# among the saved pairs that passes. Only pairs with p < MQTL_P_CIS_SAVE are
# saved, so if every saved pair passes, the true threshold lies at or beyond
# MQTL_P_CIS_SAVE: significant pairs were discarded and the counts are
# incomplete. Then this warns and returns NA; loosen MQTL_P_CIS_SAVE in
# mqtl/config.R. cis is the cis table from summarise_cis()
fdr_threshold <- function(cis, fdr = MQTL_FDR, label = "") {
  passing <- cis$fdr < fdr
  if (nrow(cis) > 0 && all(passing)) {
    warning(label, "every saved pair passes FDR < ", fdr, ": the FDR threshold is ",
            "at or beyond MQTL_P_CIS_SAVE (", MQTL_P_CIS_SAVE, "), so the counts are ",
            "incomplete. Loosen MQTL_P_CIS_SAVE in mqtl/config.R", call. = FALSE)
    return(NA_real_)
  }
  if (!any(passing)) return(NA_real_)
  max(cis$p[passing])
}

# --- Permutations (mqtl/4_mqtl_permutations.R) -------------------------------

# Rows of x (features x samples) with the covariates regressed out. Q is an
# orthonormal basis of the covariates incl. the intercept (qr.Q of
# cbind(1, covariates)), so the residuals are x - (x Q) Q'
residualize_rows <- function(x, Q) {
  x - tcrossprod(x %*% Q, Q)
}

# p-value of a correlation r (given as r^2) with df residual degrees of
# freedom: the t-test of the SNP effect, t^2 = df r^2 / (1 - r^2) ~ F(1, df).
# With methylation and genotypes residualized on the covariates, r is the
# partial correlation and df = samples - covariates - 2 (incl. the intercept)
# gives the same p-value as the full linear model (MatrixEQTL). Under the
# null r^2 ~ Beta(1/2, df/2), which gives the same p-value faster than pf
pval_from_r2 <- function(r2, df) {
  pbeta(r2, 0.5, df / 2, lower.tail = FALSE)
}

# Inverse of pval_from_r2: the r^2 that gives p-value p with df
r2_from_pval <- function(p, df) {
  qbeta(p, 0.5, df / 2, lower.tail = FALSE)
}

# Beta approximation of a CpG's permutation p-value, as FastQTL and tensorQTL
# (GTEx): the lead SNP's p-value from each permutation follows roughly a
# Beta(a, b) distribution, which gives p-values far below 1 / MQTL_N_PERM.
# 1. Effective degrees of freedom (true_df): the df at which the permuted
#    p-values give shape a = 1 (method of moments), a small correction that
#    makes the Beta form fit better (FastQTL). It is usually near or below
#    the model's df; well above it means the permuted correlations are too
#    small (e.g. permuted methylation not re-residualized)
# 2. Beta(a, b) fitted by maximum likelihood to the permuted p-values with
#    true_df, starting from the method-of-moments estimates
# 3. p_beta: the Beta CDF at the observed lead p-value with true_df
# r2_lead is the observed lead SNP's r^2, r2_perm the best r^2 of each
# permutation, df the model's residual df. Returns c(true_df, p_true_df,
# shape1, shape2, p_beta)
beta_approx_pval <- function(r2_lead, r2_perm, df) {
  mom <- function(p) {                     # method-of-moments Beta shapes
    m <- mean(p); v <- mean((p - m)^2)
    a <- m * (m * (1 - m) / v - 1)
    c(a, a * (1 / m - 1))
  }
  # Root search from [df / 2, 2 df], widened if needed; if there's no root,
  # the df closest to a = 1
  cost    <- function(d) mom(pval_from_r2(r2_perm, d))[1] - 1
  true_df <- tryCatch(uniroot(cost, c(df / 2, df * 2), extendInt = "yes", tol = 1e-4)$root,
                      error = function(e) optimize(function(d) abs(cost(d)), c(df / 10, df * 10))$minimum)

  p  <- pval_from_r2(r2_perm, true_df)
  s1 <- mean(log(p)); s2 <- mean(log1p(-p))
  # Negative log-likelihood per permutation and its gradient, on the log
  # scale of the shapes so they stay positive
  nll <- function(ls) {
    a <- exp(ls[1]); b <- exp(ls[2])
    lbeta(a, b) - (a - 1) * s1 - (b - 1) * s2
  }
  grad <- function(ls) {
    a <- exp(ls[1]); b <- exp(ls[2])
    dg <- digamma(a + b)
    c(a * (digamma(a) - dg - s1), b * (digamma(b) - dg - s2))
  }
  fit   <- optim(log(pmax(mom(p), 1e-3)), nll, grad, method = "BFGS")
  shape <- exp(fit$par)

  p_true_df <- pval_from_r2(r2_lead, true_df)
  c(true_df = true_df, p_true_df = p_true_df, shape1 = shape[1], shape2 = shape[2],
    p_beta = pbeta(p_true_df, shape[1], shape[2]))
}

# Permutation statistics for the CpGs of one chromosome. Same model as step 3
# (run_cis_chr): the SNPs in each CpG's cis window except those in its probe
# (cis_snp_ranges), missing genotypes set to the SNP's mean (as MatrixEQTL).
# Genotypes and methylation are residualized on the covariates and scaled to
# unit length, so each SNP-CpG correlation is one cross-product. The samples
# are shuffled with the columns of perms (samples x permutations, the same
# for every CpG), as in tensorQTL, and each permutation's best |r| over the
# CpG's SNPs is kept. A shuffled residual is no longer orthogonal to the
# covariates, so, as tensorQTL, it is residualized again: the genotypes are
# orthogonal to the covariates, so only the scaling changes (r divided by
# the length of the shuffled residual's part outside the covariates; without
# this the permuted r^2 is ~covariates/samples too small and the null too
# weak). Neighbouring CpGs (block at a time) share most of
# their SNPs, so they are permuted together in one matrix product.
#   geno      load_chr_genotypes(..., sliced = FALSE)
#   Ys        methylation residuals of the chromosome's CpGs, unit length
#             (CpGs x samples), y_norm their length before scaling
#   genes, region  the chromosome's rows of cpg_pos and probe_regions
#   Q         orthonormal basis of the covariates (residualize_rows)
#   df        residual df of the model (samples - covariates - 2)
# Returns one row per CpG with a SNP tested: the lead SNP (as step 3), its
# effect, t and p-value, n_snps, the empirical permutation p-value p_perm =
# (1 + permutations at least as strong) / (1 + MQTL_N_PERM), and the Beta
# approximation (beta_approx_pval)
cis_perm_chr <- function(geno, Ys, y_norm, genes, region, Q, perms, df,
                         block = MQTL_PERM_BLOCK) {
  n <- ncol(Ys)
  n_perm <- ncol(perms)

  # Genotypes in position order, missing calls set to the SNP's mean,
  # residualized and scaled to unit length (a SNP without variation is left
  # at 0, so it never leads)
  o        <- order(geno$snp_info$pos)
  G        <- geno$geno[o, , drop = FALSE]
  snp_info <- geno$snp_info[o]
  miss     <- which(is.na(G), arr.ind = TRUE)
  if (nrow(miss) > 0) G[miss] <- rowMeans(G, na.rm = TRUE)[miss[, 1]]
  G      <- residualize_rows(G, Q)
  g_norm <- sqrt(rowSums(G^2))
  G      <- G / ifelse(g_norm > 1e-8, g_norm, Inf)

  win <- cis_snp_ranges(genes, region, snp_info$pos)
  win$n_snps <- (win$last - win$first + 1) - pmax(win$probe_last - win$probe_first + 1, 0)

  # Blocks of neighbouring CpGs with a SNP to test
  cpgs   <- which(win$n_snps > 0)
  cpgs   <- cpgs[order(genes$left[cpgs])]
  blocks <- split(cpgs, ceiling(seq_along(cpgs) / block))

  out <- vector("list", length(blocks))
  for (b in seq_along(blocks)) {
    j  <- blocks[[b]]
    u  <- min(win$first[j]):max(win$last[j])          # SNPs of the block
    Gu <- G[u, , drop = FALSE]
    Yb <- Ys[j, , drop = FALSE]
    r_obs <- Gu %*% t(Yb)                              # SNPs x CpGs
    # Permuted methylation: for CpG k, columns (k - 1) * n_perm + 1:n_perm
    Yp     <- matrix(t(Yb)[perms, ], nrow = n)
    r_perm <- Gu %*% Yp                                # SNPs x (CpGs * permutations)
    # Re-residualize the shuffled methylation: Yp has unit length, and Gu is
    # orthogonal to Q, so Gu Yp is unchanged and only the length of Yp's
    # residual part, sqrt(1 - |Q'Yp|^2), is divided out
    r_perm <- r_perm / rep(sqrt(1 - colSums(crossprod(Q, Yp)^2)), each = nrow(r_perm))

    out[[b]] <- do.call(rbind, lapply(seq_along(j), function(k) {
      w    <- j[k]
      rows <- win$first[w]:win$last[w]
      rows <- rows[rows < win$probe_first[w] | rows > win$probe_last[w]] - u[1] + 1
      cols <- (k - 1) * n_perm + seq_len(n_perm)
      best <- pmax(colMaxs(r_perm, rows = rows, cols = cols),
                   -colMins(r_perm, rows = rows, cols = cols))
      lead <- rows[which.max(abs(r_obs[rows, k]))]
      r    <- r_obs[lead, k]
      s    <- u[1] + lead - 1                          # lead SNP, index into snp_info
      ba   <- beta_approx_pval(r^2, best^2, df)
      data.frame(cpg = genes$geneid[w], n_snps = length(rows), snp = snp_info$snpid[s],
                 pos = snp_info$pos[s], counted = snp_info$counted[s], alt = snp_info$alt[s],
                 r = r, beta = r * y_norm[w] / g_norm[s], t_stat = r * sqrt(df / (1 - r^2)),
                 p = pval_from_r2(r^2, df),
                 p_perm = (1 + sum(best >= abs(r))) / (1 + n_perm),
                 true_df = ba[["true_df"]], p_true_df = ba[["p_true_df"]],
                 beta_shape1 = ba[["shape1"]], beta_shape2 = ba[["shape2"]],
                 p_beta = ba[["p_beta"]])
    }))
  }
  bind_rows(out) %>% mutate(chr = genes$chr[1], .after = cpg)
}

# Storey's estimate of the proportion of CpGs without a cis-mQTL (pi0), with
# the defaults of qvalue::pi0est: pi0(lambda) = share of p-values above
# lambda / (1 - lambda) for lambda = 0.05, 0.10, ..., 0.95, smoothed with a
# cubic spline (df = 3), read at the largest lambda and capped at 1. Null
# p-values are uniform, so the flat right side of the histogram gives the
# share of nulls. Falls back to 1 (= Benjamini-Hochberg) with a warning if
# the estimate isn't positive
storey_pi0 <- function(p, lambda = seq(0.05, 0.95, 0.05)) {
  pi0 <- sapply(lambda, function(l) mean(p >= l) / (1 - l))
  pi0 <- min(predict(smooth.spline(lambda, pi0, df = 3), x = max(lambda))$y, 1)
  if (!is.finite(pi0) || pi0 <= 0) {
    warning("Storey pi0 estimate is ", pi0, "; using pi0 = 1 (Benjamini-Hochberg)", call. = FALSE)
    pi0 <- 1
  }
  pi0
}

# CpG-level significance from the permutations (perm: rows from
# cis_perm_chr), as GTEx:
#   q_storey   Storey q-values of p_beta: pi0 (storey_pi0) x the
#              Benjamini-Hochberg q-value (identical to qvalue::qvalue)
#   q_bh       Benjamini-Hochberg q-values of p_beta
#   q          the one used, MQTL_PERM_FDR_METHOD; a CpG has a cis-mQTL if
#              q < fdr
#   p_nominal_threshold  per CpG, the nominal p-value (step 3 scale) a
#              SNP needs to be significant for that CpG: the p_beta cut-off
#              between significant and non-significant CpGs, mapped back
#              through the CpG's Beta distribution (qbeta) and true_df to
#              its r^2, then to the p-value with the model's df. NA if no
#              CpG is significant
# Returns list(perm = perm with these columns, pi0, p_beta_threshold)
perm_significance <- function(perm, df, fdr = MQTL_FDR, method = MQTL_PERM_FDR_METHOD) {
  stopifnot(method %in% c("storey", "bh"))
  pi0  <- storey_pi0(perm$p_beta)
  perm <- perm %>%
    mutate(q_bh     = p.adjust(p_beta, method = "BH"),
           q_storey = pi0 * q_bh,
           q        = if (method == "storey") q_storey else q_bh)

  # Cut-off halfway between the largest significant and smallest
  # non-significant p_beta
  sig <- perm$q < fdr
  thr <- if (!any(sig)) NA_real_
         else if (all(sig)) max(perm$p_beta)
         else (max(perm$p_beta[sig]) + min(perm$p_beta[!sig])) / 2
  perm$p_nominal_threshold <- if (is.na(thr)) NA_real_ else
    pval_from_r2(r2_from_pval(qbeta(thr, perm$beta_shape1, perm$beta_shape2), perm$true_df), df)
  list(perm = perm, pi0 = pi0, p_beta_threshold = thr)
}
