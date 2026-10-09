# =============================================================================
# Configuration - mQTL module (mqtl/ scripts)
# Author: GP2 Subtypes and Mechanisms - M.E., M.P.
# Date: Oct 7, 2026
# Description: Settings used only by the mqtl/ scripts. The ancestry to
#              analyse is a command-line argument, not a setting here;
#              R/mqtl_setup.R sets it and the per-ancestry paths (GENO_PFILE,
#              MQTL_DIR, MQTL_SAMPLE_MAP, MQTL_KEEP). Source after config.R,
#              then source R/mqtl_setup.R
# =============================================================================

# GP2 master key with genotype availability, GP2 QC and ancestry labels
GP2_MASTER_KEY <- "/mnt/output/metadata/R12_CURRENT_master_key_nba_wgs_20_06_2026.txt"

# GP2 genotypes (plink2 .pgen/.pvar/.psam, hg38), one file set per ancestry:
# folder, and file prefix with {ANCESTRY} where the ancestry label goes, e.g.
# .../gwas/GP2_r12_final_samples_related_removed_AFR
GENO_PFILE_PATH <- "/home/Michael/gp2_release12/gwas"
GENO_PFILE_NAME <- "GP2_r12_final_samples_related_removed_{ANCESTRY}"
GENO_SOURCE     <- "nba"    # "wgs" or "nba" (imputed) - column in master key; must match the genotype files

MQTL_CIS_WINDOW <- 1e6     # +/- 1 Mb around each CpG

# Genotype QC (2_mqtl_genotypes.sh, plink2), applied to the samples in the
# analysis
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

# Plots (4_mqtl_plots.R)
MQTL_PLOT_N_TOP      <- 12      # top mQTLs (distinct loci) shown as genotype boxplots
MQTL_PLOT_N_REGIONAL <- 4       # top loci shown as regional plots
MQTL_PLOT_REGION     <- 250e3   # regional plot window: +/- 250 kb around the CpG
