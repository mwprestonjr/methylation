# =============================================================================
# mQTL functions
# Author: GP2 Subtypes and Mechanisms - M.E., M.P.
# Date: Oct 7, 2026
# Description: Steps shared by mqtl/3_mqtl.R and mqtl/tune_meth_pcs.R:
#              loading and preparing the methylation data and covariates,
#              latent methylation PCs, and the per-chromosome cis-mQTL scan.
#              Needs config.R, mqtl/config.R, R/mqtl_setup.R and
#              R/array_profiles.R sourced first, and minfi, tidyverse,
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
  # Probe coordinates from the array's annotation (one array type per data source)
  array_type <- unique(targets$Array)
  stopifnot(length(array_type) == 1, array_type %in% names(ARRAY_PROFILES))
  profile <- ARRAY_PROFILES[[array_type]]
  library(profile$anno_pkg, character.only = TRUE)
  cat("\nArray:", array_type, "- annotation:", profile$anno_pkg, "(", profile$genome, ")\n")

  ann <- getAnnotation(get(profile$anno_pkg))
  ann <- ann[rownames(meth), c("chr", "pos")]
  gr  <- GRanges(ann$chr, IRanges(ann$pos, width = 1), cpg = rownames(ann))

  if (profile$genome == "hg38") {
    gr_hg38 <- gr
  } else if (profile$genome == "hg19") {
    cat("Lifting probe coordinates hg19 -> hg38...\n")
    # import.chain() needs an uncompressed file
    chain_file <- tempfile(fileext = ".chain")
    system2("gunzip", c("-c", CHAIN_HG19_HG38), stdout = chain_file)
    chain   <- import.chain(chain_file)
    gr_hg38 <- liftOver(gr, chain)

    # Keep probes that map to exactly one hg38 position on the same chromosome
    # (compared without the "chr" prefix: the chain file may name chromosomes
    # "19" where the annotation has "chr19")
    one_hit  <- lengths(gr_hg38) == 1
    gr_hg38  <- unlist(gr_hg38[one_hit])
    same_chr <- sub("^chr", "", as.character(seqnames(gr_hg38))) ==
                sub("^chr", "", ann$chr[one_hit])
    gr_hg38  <- gr_hg38[same_chr]
    cat("Probes lifted:", length(gr_hg38), "of", nrow(meth),
        "(dropped", nrow(meth) - length(gr_hg38), ")\n")
  } else {
    stop("No hg38 conversion for genome build ", profile$genome)
  }

  cpg_pos <- data.frame(
    geneid = gr_hg38$cpg,
    chr    = sub("^chr", "", as.character(seqnames(gr_hg38))),
    left   = start(gr_hg38),
    right  = start(gr_hg38)
  )
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

  list(meth_int = meth_int, cpg_pos = cpg_pos, known_covs = known_covs,
       sample_map = sample_map, geno_dir = geno_dir)
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
# order of meth_int: list(snps = SlicedData, snp_info = data.table)
load_chr_genotypes <- function(geno_dir, chr, sample_map, sample_names) {
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

  snps <- SlicedData$new()
  snps$CreateFromMatrix(geno_mat)
  snps$ResliceCombined(sliceSize = 5000)
  list(snps = snps, snp_info = snp_info)
}

# cis-mQTL scan of one chromosome: tests within +/- MQTL_CIS_WINDOW, keeps
# pairs with p < MQTL_P_CIS_SAVE. Returns list(eqtls = data frame, ntests)
run_cis_chr <- function(geno, meth_int, cpg_pos, chr, cvrt) {
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
  list(eqtls = eqtls, ntests = me$cis$ntests)
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
fdr_threshold <- function(cis, fdr = 0.05, label = "") {
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
