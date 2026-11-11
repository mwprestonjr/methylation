# =============================================================================
# Configuration - mQTL module (mqtl/ scripts)
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 7, 2026
# Updated: Nov 11, 2026
# Description: Settings used only by the mqtl/ scripts. The ancestry to
#              analyze is a command-line argument, not a setting here;
#              R/mqtl_setup.R sets it and the per-ancestry paths (GENO_PFILE,
#              MQTL_DIR, MQTL_SAMPLE_MAP, MQTL_KEEP). Source after config.R,
#              then source R/mqtl_setup.R
# =============================================================================

# GP2 release bucket, mounted with gcsfuse. A full path, not ~/: the shell
# step (2_mqtl_genotypes.sh) passes these paths to plink2 as they are
GP2_RELEASE_DIR <- "/home/Michael/gp2_release12"

# GP2 NBA array genotypes (plink2 .pgen/.pvar/.psam, hg38), one file set per
# ancestry: folder, and file prefix with {ANCESTRY} where the ancestry label
# goes, e.g. .../nba/gwas/GP2_r12_final_samples_related_removed_AFR. These are
# the directly typed SNPs (~1.6M, not imputed), with related individuals
# removed. They define the samples and the genotype PCs in every run, and are
# the scan genotypes when MQTL_GENOTYPES is "array"
GENO_PFILE_PATH <- file.path(GP2_RELEASE_DIR, "nba/gwas")
GENO_PFILE_NAME <- "GP2_r12_final_samples_related_removed_{ANCESTRY}"
GENO_SOURCE     <- "nba"    # "wgs" or "nba" (NBA array) - column in master key; must match the genotype files

# Genotypes tested in the cis scan:
#   "imputed": GP2 imputed dosages (MQTL_IMPUTED_VCFS). The default: the
#              array has ~1.6M SNPs, too few for colocalization with GWAS
#              loci and for comparison with imputation-based mQTL resources
#              (e.g. GoDMC)
#   "array":   the NBA array calls (GENO_PFILE_*), as in the first analysis
# Imputed runs write to <ancestry>_imputed (R/mqtl_setup.R), so the array
# results in <ancestry> are kept for comparison
MQTL_GENOTYPES <- "imputed"

# Imputed dosages: Minimac4 VCFs (hg38, dosage in the DS field, imputation
# R2 in INFO), one file per chromosome ({CHR}). Most ancestries are in
# imputed_vcf_filtered_r2/ (R2 >= 0.3 already applied). R12 has no EUR there:
# EUR was imputed in four sample batches, whose raw output (no R2 filter) is in
# imputation_out_dir/ (the EUR_merged folder there is empty). Each ancestry is
# a vector of files per chromosome, one per batch; "default" covers the rest
MQTL_IMPUTED_VCFS <- list(
  EUR     = file.path(GP2_RELEASE_DIR, "imputation_out_dir",
                      paste0("EUR_part", 1:4), "chr{CHR}.dose.vcf.gz"),
  default = file.path(GP2_RELEASE_DIR, "imputed_vcf_filtered_r2",
                      "{ANCESTRY}", "chr{CHR}.vcf.gz"))

# Imputed variants kept: imputation R2 >= MQTL_MIN_R2, GP2's own threshold for
# imputed_vcf_filtered_r2, so EUR (filtered here) and the other ancestries
# (filtered by GP2) are treated alike. R2 is from the INFO field, estimated
# by Minimac4 over each imputation batch (~22,000 samples for EUR), so it is
# far more precise than anything estimated in our samples. A variant must
# pass in every batch our samples are in
MQTL_MIN_R2 <- 0.3
# Variants with MAF below this in the imputation batch (INFO MAF) are dropped
# before export, to keep the files small: such a variant reaching MQTL_MAF in
# our samples would be a sampling fluke. MQTL_MAF itself is applied afterwards,
# in our samples
MQTL_IMPUTED_MIN_BATCH_MAF <- 0.01

MQTL_CIS_WINDOW <- 1e6     # +/- 1 Mb around each CpG

# Genotype QC (2_mqtl_genotypes.sh, plink2), applied to the samples in the
# analysis. The array genotypes get all four; imputed dosages get MQTL_MAF and
# MQTL_GENO_MISSING (R/merge_imputed_dosages.R): they have no missing calls
# per sample, and Hardy-Weinberg is a check on genotype calling, done by GP2 on
# the array data before imputation
MQTL_MAF            <- 0.05    # drop variants with minor allele frequency < 5% (--maf)
MQTL_GENO_MISSING   <- 0.05    # drop variants missing in > 5% of samples (--geno)
MQTL_SAMPLE_MISSING <- 0.05    # drop samples missing > 5% of genotypes (--mind)
MQTL_HWE_P          <- 1e-6    # drop variants failing Hardy-Weinberg at p < 1e-6 (--hwe)

# LD pruning for the genotype PCs (--indep-pairwise <window> <r2>)
MQTL_PRUNE_WINDOW <- "500kb"
MQTL_PRUNE_R2     <- 0.1

MQTL_N_GENO_PCS <- 5
MQTL_N_METH_PCS <- 10      # latent methylation PCs; choose with mqtl/tune_meth_pcs.R
MQTL_N_METH_PCS_GRID <- c(0, 5, 10, 15, 20, 30)   # values tested by mqtl/tune_meth_pcs.R

# Multiple testing: a CpG has a cis-mQTL if its lead SNP passes MQTL_FDR
# (Benjamini-Hochberg over all cis tests). MQTL_P_STRICT gives a stricter,
# fixed-threshold set comparable with GoDMC (cis p < 1e-8)
MQTL_FDR        <- 0.05
MQTL_P_STRICT   <- 1e-8
MQTL_P_CIS_SAVE <- 1e-3    # pairs with p below this are saved; must be looser than the FDR threshold (3_mqtl.R checks)

# Permutations (4_mqtl_permutations.R), as FastQTL/tensorQTL (GTEx): each
# CpG's lead SNP is compared with the best SNP after shuffling the samples,
# which accounts for the number of SNPs tested and the LD between them. A CpG
# has a cis-mQTL if the q-value of its Beta-approximated permutation p-value
# is < MQTL_FDR
MQTL_N_PERM          <- 1000       # shuffles, the same for every CpG
MQTL_PERM_SEED       <- 20261027
MQTL_PERM_FDR_METHOD <- "storey"   # "storey" (q-values, as GTEx) or "bh" (Benjamini-Hochberg, more conservative)
MQTL_PERM_BLOCK      <- 20         # neighboring CpGs permuted together; memory ~ SNPs in their windows x MQTL_PERM_BLOCK x MQTL_N_PERM x 8 bytes, per core
# Chromosomes run in parallel. Memory (Psomagen AFR, ~870k CpGs, array
# genotypes): ~8 GB for the main process plus ~4 GB per core, so ~70 GB with
# 15 cores; lower this if other jobs share the VM. Imputed genotypes have ~15x
# more SNPs per window, so expect much more memory per core and a longer run:
# start with fewer cores and watch memory on the first chromosomes
MQTL_PERM_CORES      <- max(1, parallel::detectCores() - 1)

# Plots (5_mqtl_plots.R)
MQTL_PLOT_N_TOP      <- 12      # top mQTLs (distinct loci) shown as genotype boxplots
MQTL_PLOT_N_REGIONAL <- 4       # top loci shown as regional plots
MQTL_PLOT_REGION     <- 250e3   # regional plot window: +/- 250 kb around the CpG
MQTL_PLOT_LOGP_CAP   <- 50      # genome-wide plot: -log10 p capped here so extreme hits don't flatten the rest
MQTL_PLOT_ZOOM       <- 50e3    # zoomed lead SNP distance plot: +/- 50 kb around the CpG
