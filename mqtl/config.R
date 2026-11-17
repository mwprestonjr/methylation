# =============================================================================
# Configuration - mQTL module (mqtl/ scripts)
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 7, 2026
# Updated: Nov 12, 2026
# Description: Settings used only by the mqtl/ scripts. The ancestry to
#              analyze is a command-line argument, not a setting here;
#              R/mqtl_setup.R sets it and the per-ancestry paths (GENO_PFILE,
#              MQTL_DIR, MQTL_SAMPLE_MAP, MQTL_KEEP). Source after config.R,
#              then source R/mqtl_setup.R
# =============================================================================

# Genotypes tested, from GP2 R12 (plink2 .pgen/.pvar/.psam, hg38, one file
# set per ancestry):
#   "wgs":   whole-genome sequencing (DRAGEN), the default. Variants are
#            sequenced, so there's no imputation uncertainty, and the common
#            variants needed for colocalization with GWAS loci and for
#            comparison with mQTL resources are all there. 268 of 283 PPMI EUR
#            and all 268 Psomagen AFR methylation samples have WGS
#   "array": the NBA array (~1.6M typed variants, not imputed), as in the
#            first analysis
# WGS runs write to <ancestry>_wgs (R/mqtl_setup.R), so the array results in
# <ancestry> are kept for comparison. GP2's imputed data were not used: the
# EUR imputation output in R12 is raw, unfiltered VCFs for ~22,000 samples per
# batch (5.6 TB), which would take days to read
MQTL_GENOTYPES <- "wgs"

# GP2 release bucket, mounted with gcsfuse. A full path, not ~/: the shell
# step (2_mqtl_genotypes.sh) passes these paths to plink2 as they are
GP2_RELEASE_DIR <- "/home/Michael/gp2_release12"

# File prefixes, with {ANCESTRY} where the ancestry label goes
GENO_PFILE_ARRAY <- file.path(GP2_RELEASE_DIR,
                              "nba/gwas/GP2_r12_final_samples_related_removed_{ANCESTRY}")
GENO_PFILE_WGS   <- file.path(GP2_RELEASE_DIR,
                              "dragen/r12_analyses/plink2_pfiles/{ANCESTRY}/all_chrs")

# WGS samples with a chromosomal aneuploidy (e.g. XXY, trisomy): left out, since
# their genotype dosages on the affected chromosome aren't 0/1/2
WGS_ANEUPLOIDY <- file.path(GP2_RELEASE_DIR, "dragen/r12_wgs_samples_with_aneuploidy.csv")

# GP2 master key columns for the genotype source (genotyped flag, QC prune
# reason <source>_prune_reason, ancestry <source>_label); follows MQTL_GENOTYPES
GENO_SOURCE <- switch(MQTL_GENOTYPES, wgs = "wgs", array = "nba")

MQTL_CIS_WINDOW <- 1e6     # +/- 1 Mb around each CpG

# Genotype QC (2_mqtl_genotypes.sh, plink2), applied to the samples in the
# analysis. WGS calls must also have FILTER = PASS
MQTL_MAF            <- 0.05    # drop variants with minor allele frequency < 5% (--maf)
MQTL_GENO_MISSING   <- 0.05    # drop variants missing in > 5% of samples (--geno)
MQTL_SAMPLE_MISSING <- 0.05    # drop samples missing > 5% of genotypes (--mind)
MQTL_HWE_P          <- 1e-6    # drop variants failing Hardy-Weinberg at p < 1e-6 (--hwe)

# Related individuals among the analyzed samples (2_mqtl_genotypes.sh): plink2
# KING kinship on the QC'd, LD-pruned genotypes; one of each pair above the
# cutoff is left out (--king-cutoff), and steps 3-5 use only the samples kept.
# 0.0884 = second-degree relatives or closer (KING's threshold). Checked within
# our samples rather than with a GP2-wide list: only relatives within the
# analysis matter, and GP2's WGS sample lists (dragen/r12_analyses/sample_lists)
# cover one analysis's samples, not ours. The array files already have GP2's
# relatives removed, so in array runs this is a check
MQTL_KING_CUTOFF <- 0.0884

# LD pruning for the genotype PCs and KING (--indep-pairwise <window> <r2>)
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
# Chromosomes run in parallel, one per core, largest first. Memory: ~8 GB for
# the main process plus, per core, ~4 GB with array genotypes (Psomagen AFR:
# ~70 GB with 15 cores) or up to ~5 GB with WGS (chr1, PPMI EUR). The run
# takes at least as long as chr1 (WGS PPMI EUR: ~5 h), so with WGS more than
# ~10 cores saves little time. Lower this if other jobs share the VM; a run
# that fails resumes from the chromosomes already saved
MQTL_PERM_CORES      <- min(10, max(1, parallel::detectCores() - 1))

# Plots (5_mqtl_plots.R)
MQTL_PLOT_N_TOP      <- 12      # top mQTLs (distinct loci) shown as genotype boxplots
MQTL_PLOT_N_REGIONAL <- 4       # top loci shown as regional plots
MQTL_PLOT_REGION     <- 250e3   # regional plot window: +/- 250 kb around the CpG
MQTL_PLOT_LOGP_CAP   <- 50      # genome-wide plot: -log10 p capped here so extreme hits don't flatten the rest
MQTL_PLOT_ZOOM       <- 50e3    # zoomed lead SNP distance plot: +/- 50 kb around the CpG
