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

# hg19 -> hg38 liftover for arrays whose annotation is hg19 (EPICv1)
CHAIN_HG19_HG38 <- "/mnt/expansion_working/chain_files/hg19ToHg38.chain.gz"

MQTL_CIS_WINDOW <- 1e6     # +/- 1 Mb around each CpG
MQTL_MAF        <- 0.05
MQTL_N_GENO_PCS <- 5
MQTL_N_METH_PCS <- 10      # latent methylation PCs; choose with mqtl/tune_meth_pcs.R
MQTL_N_METH_PCS_GRID <- c(0, 5, 10, 15, 20, 30)   # values tested by mqtl/tune_meth_pcs.R
MQTL_P_CIS_SAVE <- 1e-5    # nominal p-value threshold to write cis results
