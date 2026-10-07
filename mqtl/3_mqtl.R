# =============================================================================
# mQTL Step 3: cis-mQTL mapping
# Author: GP2 Subtypes and Mechanisms - M.E.
# Date: Oct 7, 2026
# Description: Tests SNP-CpG associations within +/- MQTL_CIS_WINDOW using
#              MatrixEQTL, one chromosome at a time. Methylation: ComBat
#              M values, inverse-normal transformed per CpG. Covariates: age,
#              sex, phenotype, cell proportions, genotype PCs and latent
#              methylation PCs. Probe coordinates come from the array's
#              annotation (R/array_profiles.R) and are lifted hg19 -> hg38
#              when needed (EPICv1) to match GP2 genotypes.
# Usage:       Rscript mqtl/3_mqtl.R <data source> <ancestry>   (see R/mqtl_setup.R)
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

library(minfi)                            # getAnnotation
library(tidyverse)
library(data.table)
library(matrixStats)                      # rowVars
library(MatrixEQTL)
library(GenomicRanges)
library(rtracklayer)                      # liftOver

# Load shared configuration and array profiles
source("config.R")
source("R/mqtl_setup.R")
source("R/array_profiles.R")

GENO_DIR <- file.path(MQTL_DIR, "genotypes")

# --- 1. Load inputs ----------------------------------------------------------

cat("Loading inputs...\n")
combat_mVals <- readRDS(file.path(DIR_RESULTS, "combat_mVals.rds"))
targets      <- read.csv(file.path(DIR_RESULTS, "sample_sheet_final.csv"),
                         colClasses = c(GP2ID       = "character",
                                        clinical_id = "character",
                                        Sentrix_ID  = "character",
                                        Basename    = "character")) %>%
  mutate(meth_id = basename(Basename))
sample_map   <- read.csv(MQTL_SAMPLE_MAP, colClasses = "character")

# Genotype PCs; plink2 writes either "#IID" or "#FID IID" as the header
geno_pcs <- fread(file.path(GENO_DIR, "geno_pcs.eigenvec"))
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
covs <- sample_map %>%
  select(meth_id) %>%
  left_join(targets, by = "meth_id") %>%
  left_join(geno_pcs, by = c("GP2ID" = "IID"))
stopifnot(all(covs$meth_id == colnames(meth)))

# --- 2. CpG positions (hg38) -------------------------------------------------

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

# --- 3. Methylation transform and latent factors -----------------------------

cat("\nInverse-normal transforming M values per CpG...\n")
int_transform <- function(x) qnorm((rank(x, ties.method = "average") - 0.5) / length(x))
meth_int <- t(apply(meth, 1, int_transform))
dimnames(meth_int) <- dimnames(meth)

# Known covariates (drop Neu: cell proportions sum to ~1). Phenotype is
# coded as indicator columns against the most common phenotype, so it works
# for PD vs Control and for more groups (PD, DLB, MSA, ...); it is left out if
# all samples share one phenotype
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

# Latent methylation PCs from residuals after known covariates.
# These soak up unmeasured technical/biological variation and boost power;
# MQTL_N_METH_PCS is worth tuning (e.g. 0/5/10/20) by number of cis-mQTLs found.
cat("Computing", MQTL_N_METH_PCS, "latent methylation PCs...\n")
X        <- model.matrix(~ ., data = known_covs)
resid    <- meth_int - t(X %*% solve(crossprod(X), crossprod(X, t(meth_int))))
top_var  <- order(rowVars(resid), decreasing = TRUE)[1:min(50000, nrow(resid))]
meth_pca <- prcomp(t(resid[top_var, ]), center = TRUE, scale. = FALSE,
                   rank. = MQTL_N_METH_PCS)
meth_pcs <- meth_pca$x
colnames(meth_pcs) <- paste0("meth_PC", seq_len(ncol(meth_pcs)))
rm(resid)

all_covs <- cbind(known_covs, meth_pcs)
cat("Covariates in model:", ncol(all_covs), "\n")
print(colnames(all_covs))

# MatrixEQTL wants covariates as rows x samples
cvrt <- SlicedData$new()
cvrt$CreateFromMatrix(t(as.matrix(all_covs)))
colnames(cvrt) <- colnames(meth_int)

# --- 4. cis-mQTL mapping, per chromosome -------------------------------------

cat("\n--- cis-mQTL mapping (window +/-", MQTL_CIS_WINDOW / 1e6, "Mb) ---\n")

cis_results <- list()
n_tests     <- 0

for (chr in 1:22) {
  cat("\nchr", chr, "\n", sep = "")

  # Genotypes: .traw sample columns are "<FID>_<IID>"
  traw     <- fread(file.path(GENO_DIR, paste0("geno_chr", chr, ".traw")))
  geno_ids <- names(traw)[-(1:6)]
  col_idx  <- sapply(sample_map$GP2ID,
                     function(id) which(endsWith(geno_ids, paste0("_", id))))
  stopifnot(is.numeric(col_idx), length(col_idx) == nrow(sample_map))

  snp_info <- traw[, .(snpid = SNP, chr = as.character(CHR), pos = POS,
                       counted = COUNTED, alt = ALT)]
  geno_mat <- as.matrix(traw[, geno_ids[col_idx], with = FALSE])
  dimnames(geno_mat) <- list(snp_info$snpid, colnames(meth_int))
  rm(traw)

  snps <- SlicedData$new()
  snps$CreateFromMatrix(geno_mat)
  snps$ResliceCombined(sliceSize = 5000)
  rm(geno_mat)

  chr_cpgs <- cpg_pos$chr == as.character(chr)
  gene <- SlicedData$new()
  gene$CreateFromMatrix(meth_int[chr_cpgs, , drop = FALSE])
  gene$ResliceCombined(sliceSize = 2000)

  me <- Matrix_eQTL_main(
    snps                  = snps,
    gene                  = gene,
    cvrt                  = cvrt,
    output_file_name      = NULL,          # trans off
    pvOutputThreshold     = 0,
    output_file_name.cis  = NULL,
    pvOutputThreshold.cis = MQTL_P_CIS_SAVE,
    useModel              = modelLINEAR,
    errorCovariance       = numeric(),
    snpspos               = snp_info[, .(snpid, chr, pos)],
    genepos               = cpg_pos[chr_cpgs, ],
    cisDist               = MQTL_CIS_WINDOW,
    pvalue.hist           = FALSE,
    min.pv.by.genesnp     = FALSE,
    noFDRsaveMemory       = FALSE,
    verbose               = FALSE
  )

  n_tests <- n_tests + me$cis$ntests
  cat("  cis tests:", me$cis$ntests,
      "| saved (p <", MQTL_P_CIS_SAVE, "):", nrow(me$cis$eqtls), "\n")

  cis_results[[chr]] <- me$cis$eqtls %>%
    transmute(cpg = as.character(gene), snp = as.character(snps),
              beta, t_stat = statistic, p = pvalue) %>%
    left_join(as.data.frame(snp_info)[, c("snpid", "pos", "counted", "alt")],
              by = c("snp" = "snpid")) %>%
    mutate(chr = chr)

  rm(snps, gene, me); gc()
}

# --- 5. Multiple testing and summary -----------------------------------------

cat("\n--- Combining results ---\n")
cat("Total cis tests:", n_tests, "\n")

# Saved rows are exactly the smallest p-values, so BH with n = all tests
# gives the same q-values as adjusting the full set
cis <- bind_rows(cis_results) %>%
  left_join(cpg_pos %>% select(cpg = geneid, cpg_pos = left), by = "cpg") %>%
  mutate(distance = pos - cpg_pos,
         fdr      = p.adjust(p, method = "BH", n = n_tests)) %>%
  arrange(p)

# Lead SNP per CpG
lead <- cis %>%
  group_by(cpg) %>%
  slice_min(p, n = 1, with_ties = FALSE) %>%
  ungroup()

sig_fdr <- lead %>% filter(fdr < 0.05)
sig_1e8 <- lead %>% filter(p < 1e-8)                  # GoDMC cis threshold
cat("CpGs with a cis-mQTL at FDR < 5%:", nrow(sig_fdr), "\n")
cat("CpGs with a cis-mQTL at p < 1e-8:", nrow(sig_1e8), "\n")

# Distance of lead SNPs from their CpG
png(file.path(MQTL_DIR, "mqtl_01_lead_snp_distance.png"),
    width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = FIG_RES)
print(ggplot(sig_fdr, aes(x = distance / 1e3)) +
  geom_histogram(bins = 100) +
  labs(title = "Lead cis-mQTL SNP position relative to CpG (FDR < 5%)",
       x = "SNP - CpG distance (kb)", y = "CpGs") +
  theme_bw())
dev.off()

# --- 6. Save outputs ---------------------------------------------------------

cat("\n--- Saving outputs ---\n")
fwrite(cis,  file.path(MQTL_DIR, "cis_mqtl_all_saved.tsv.gz"), sep = "\t")
fwrite(lead, file.path(MQTL_DIR, "cis_mqtl_lead_per_cpg.tsv"),  sep = "\t")
write.csv(all_covs %>% mutate(meth_id = colnames(meth_int), .before = 1),
          file.path(MQTL_DIR, "mqtl_covariates.csv"), row.names = FALSE)

mqtl_summary <- data.frame(
  metric = c("Samples", "Ancestry", "CpGs tested", "Cis window (bp)",
             "Genotype PCs", "Methylation PCs", "Total cis tests",
             "CpGs with cis-mQTL (FDR < 0.05)",
             "CpGs with cis-mQTL (p < 1e-8)"),
  value  = c(ncol(meth_int), paste(MQTL_ANCESTRY, collapse = "+"), nrow(meth_int), MQTL_CIS_WINDOW,
             MQTL_N_GENO_PCS, MQTL_N_METH_PCS, n_tests,
             nrow(sig_fdr), nrow(sig_1e8))
)
write.csv(mqtl_summary, file.path(MQTL_DIR, "mqtl_summary.csv"), row.names = FALSE)

cat("\nmQTL step 3 complete!\n")
cat("Outputs saved to:", MQTL_DIR, "\n")
