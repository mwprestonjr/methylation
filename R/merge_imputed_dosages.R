# =============================================================================
# Merge imputed dosages across imputation batches (mQTL step 2)
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Nov 11, 2026
# Updated: Nov 17, 2026
# Description: Called by mqtl/2_mqtl_genotypes.sh for each chromosome when
#              MQTL_GENOTYPES is "imputed". Reads the per-batch dosage files
#              (genotypes/imputed_part<k>_chr<N>.traw, our samples in each
#              imputation batch, already filtered on R2 and batch MAF), keeps
#              the variants present in every batch, joins the samples, keeps
#              the unrelated ones (KING, genotypes/king.king.cutoff.in.id),
#              applies MQTL_GENO_MISSING and MQTL_MAF in them, and
#              writes genotypes/geno_chr<N>.traw for steps 3-5. plink2 can't
#              merge files with different samples (non-concatenating
#              --pmerge is not implemented), hence this step. Appends the
#              variant counts to genotypes/imputed_variant_counts.tsv
# Usage:       Rscript R/merge_imputed_dosages.R <data source> <ancestry> <chr>
#              run from the repository root (by 2_mqtl_genotypes.sh)
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))

source("config.R")
source("mqtl/config.R")
source("R/mqtl_setup.R")

chr      <- commandArgs(trailingOnly = TRUE)[3]
geno_dir <- file.path(MQTL_DIR, "genotypes")
if (is.na(chr)) stop("Usage: Rscript R/merge_imputed_dosages.R <data source> <ancestry> <chr>")

# --- 1. Read the batches -----------------------------------------------------

# .traw: CHR SNP (C)M POS COUNTED ALT, then one dosage column per sample
# (copies of COUNTED). Variant IDs are chr:pos:ref:alt, so the same ID is the
# same variant and allele coding in every batch
part_files <- Sys.glob(file.path(geno_dir, paste0("imputed_part*_chr", chr, ".traw")))
if (length(part_files) == 0) {
  stop("No imputed_part*_chr", chr, ".traw in ", geno_dir,
       " - run bash mqtl/2_mqtl_genotypes.sh ", DATA_SOURCE, " ", MQTL_ANCESTRY)
}
parts    <- lapply(part_files, fread)
n_part   <- vapply(parts, nrow, integer(1))
info_col <- c("CHR", "SNP", "(C)M", "POS", "COUNTED", "ALT")

# --- 2. Variants in every batch, samples joined ------------------------------

# R2 is estimated per batch, so a borderline variant can pass in one and not
# another; keeping only variants in all batches avoids variants missing for
# whole batches of samples
common <- Reduce(intersect, lapply(parts, `[[`, "SNP"))
merged <- parts[[1]][match(common, SNP)]
for (k in seq_along(parts)[-1]) {
  p <- parts[[k]][match(common, SNP)]
  stopifnot(identical(p$COUNTED, merged$COUNTED), identical(p$ALT, merged$ALT))
  merged <- cbind(merged, p[, setdiff(names(p), info_col), with = FALSE])
}
rm(parts, p); invisible(gc())
setorder(merged, POS)
samples <- setdiff(names(merged), info_col)
stopifnot(!anyDuplicated(samples))

# Unrelated samples only (KING on the array genotypes, step 2). .traw sample
# columns are "<FID>_<IID>" with the VCF ID as IID (the GP2ID, or
# "<GP2ID>_<GP2ID>"), so they're matched on ending in "_<GP2ID>"; GP2IDs can
# contain "_" themselves
king_file <- file.path(geno_dir, "king.king.cutoff.in.id")
if (!file.exists(king_file)) {
  stop(king_file, " not found - run bash mqtl/2_mqtl_genotypes.sh ", DATA_SOURCE, " ", MQTL_ANCESTRY)
}
king      <- fread(king_file)
unrelated <- king[[grep("IID$", names(king))]]
keep      <- vapply(samples, function(s) any(endsWith(s, paste0("_", unrelated))), logical(1))
if (any(!keep)) {
  merged  <- merged[, !samples[!keep], with = FALSE]
  samples <- samples[keep]
}

# --- 3. QC in our samples ----------------------------------------------------

dos      <- as.matrix(merged[, samples, with = FALSE])
missing  <- rowMeans(is.na(dos))
freq     <- rowMeans(dos, na.rm = TRUE) / 2
maf      <- pmin(freq, 1 - freq)
pass     <- missing <= MQTL_GENO_MISSING & !is.na(maf) & maf >= MQTL_MAF
merged   <- merged[pass]

fwrite(merged, file.path(geno_dir, paste0("geno_chr", chr, ".traw")),
       sep = "\t", na = "NA", quote = FALSE)
unlink(part_files)

# --- 4. Counts ---------------------------------------------------------------

counts <- data.table(chr = chr, batches = length(part_files),
                     samples = length(samples),
                     variants_per_batch = paste(n_part, collapse = ","),
                     in_all_batches = length(common),
                     after_qc = nrow(merged))
counts_file <- file.path(geno_dir, "imputed_variant_counts.tsv")
fwrite(counts, counts_file, sep = "\t", append = file.exists(counts_file))
cat("chr", chr, ":", length(samples), "samples;", paste(n_part, collapse = ", "),
    "variants per batch,", length(common), "in all,", nrow(merged),
    "after missingness and MAF >=", MQTL_MAF, "\n")
