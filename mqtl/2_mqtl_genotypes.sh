#!/usr/bin/env bash
# =============================================================================
# mQTL Step 2: Genotype QC, ancestry PCs and per-chromosome dosage export
# Author: GP2 Subtypes and Mechanisms - M.E.
# Date: Oct 7, 2026
# Description: Subsets GP2 genotypes to the samples in the step 1 keep file,
#              applies variant QC, computes genotype PCs on LD-pruned SNPs,
#              and exports additive dosages (.traw, SNPs x samples) per
#              autosome for MatrixEQTL in step 3
# Usage:       bash mqtl/2_mqtl_genotypes.sh <data source> <ancestry>
#              e.g. psomagen AFR (see R/mqtl_setup.R)
#              run from the repository root
# =============================================================================

set -euo pipefail

DATA_SOURCE=${1:?"usage: bash mqtl/2_mqtl_genotypes.sh <data source> <ancestry>"}
ANCESTRY=${2:?"usage: bash mqtl/2_mqtl_genotypes.sh <data source> <ancestry>"}

# Pull paths/parameters from the R configs so there is one source of truth.
# One R call prints them as shell assignments (VAR='value'), which are then
# evaluated; capture.output hides the configs' own messages
CFG_VARS="GENO_PFILE MQTL_DIR MQTL_KEEP MQTL_MAF MQTL_GENO_MISSING MQTL_SAMPLE_MISSING MQTL_HWE_P MQTL_PRUNE_WINDOW MQTL_PRUNE_R2 MQTL_N_GENO_PCS"
CFG=$(Rscript -e "DATA_SOURCE <- '${DATA_SOURCE}'; MQTL_ANCESTRY <- '${ANCESTRY}'; invisible(capture.output(suppressMessages({source('config.R'); source('mqtl/config.R'); source('R/mqtl_setup.R')}))); for (v in strsplit('${CFG_VARS}', ' ')[[1]]) cat(v, '=', shQuote(as.character(get(v))), '\\n', sep = '')")
eval "${CFG}"
KEEP="${MQTL_KEEP}"
THREADS=${THREADS:-8}

GENO_DIR="${MQTL_DIR}/genotypes"
mkdir -p "${GENO_DIR}"

echo "Genotypes: ${GENO_PFILE}"
if [ ! -f "${GENO_PFILE}.pgen" ]; then
  echo "ERROR: ${GENO_PFILE}.pgen not found - check GENO_PFILE_PATH/GENO_PFILE_NAME in mqtl/config.R"
  exit 1
fi
echo "Keep file: ${KEEP} ($(($(wc -l < "${KEEP}") - 1)) samples)"

# --- 1. Sample subset + variant QC ------------------------------------------
# Thresholds (--geno, --mind, --maf, --hwe) are set in mqtl/config.R. The other
# options are fixed because later steps rely on them:
#   --autosome, --max-alleles 2: genotypes are coded as 0/1/2 copies of one
#     allele, which needs two copies of each chromosome and two alleles
#   --snps-only just-acgt: single-base A/C/G/T variants only, so alleles match
#     simply across datasets
#   --set-all-var-ids: variants are named chr:pos:ref:alt (hg38); the NBA array
#     assays some variants more than once, so --rm-dup keeps one copy of each
#     (removed IDs: geno_qc.rmdup.list)
echo "Running variant QC..."
plink2 --pfile "${GENO_PFILE}" \
  --keep "${KEEP}" \
  --autosome \
  --snps-only just-acgt \
  --max-alleles 2 \
  --geno "${MQTL_GENO_MISSING}" \
  --mind "${MQTL_SAMPLE_MISSING}" \
  --maf "${MQTL_MAF}" \
  --hwe "${MQTL_HWE_P}" \
  --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 50 \
  --rm-dup force-first list \
  --threads "${THREADS}" \
  --make-pgen \
  --out "${GENO_DIR}/geno_qc"

# --- 2. Genotype PCs (population structure covariates) ----------------------
echo "Computing genotype PCs..."
plink2 --pfile "${GENO_DIR}/geno_qc" \
  --indep-pairwise "${MQTL_PRUNE_WINDOW}" "${MQTL_PRUNE_R2}" \
  --threads "${THREADS}" \
  --out "${GENO_DIR}/ld_prune"

plink2 --pfile "${GENO_DIR}/geno_qc" \
  --extract "${GENO_DIR}/ld_prune.prune.in" \
  --pca "${MQTL_N_GENO_PCS}" \
  --threads "${THREADS}" \
  --out "${GENO_DIR}/geno_pcs"

# --- 3. Per-chromosome additive dosages -------------------------------------
# .traw columns: CHR SNP (C)M POS COUNTED ALT <FID_IID ...>
for chr in $(seq 1 22); do
  echo "Exporting chr${chr}..."
  plink2 --pfile "${GENO_DIR}/geno_qc" \
    --chr "${chr}" \
    --export A-transpose \
    --threads "${THREADS}" \
    --out "${GENO_DIR}/geno_chr${chr}"
done

echo "Genotype prep complete. Outputs in ${GENO_DIR}"
echo "Next: Rscript mqtl/3_mqtl.R ${DATA_SOURCE} ${ANCESTRY}"
