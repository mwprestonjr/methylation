# =============================================================================
# Preprocessing 4: CpG annotation
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 8, 2026
# Description: Annotates every probe on the data source's array with its hg38
#              position, CpG island context and GENCODE gene (promoter, gene
#              body or nearest gene; R/cpg_annotation.R), and saves
#              CPG_ANNOTATION for later analyses (e.g. mQTL plots and
#              top-hits tables). Covers all probes on the array, so it doesn't
#              depend on QC; rerun only when the gene settings in config.R change
# Usage:       Rscript preprocessing/4_cpg_annotation.R <data source>   (see config.R)
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

library(minfi)                            # getAnnotation
library(data.table)
library(GenomicRanges)
library(rtracklayer)                      # liftOver

# Load shared configuration, array profiles and CpG annotation helpers
source("config.R")
source("preprocessing/config.R")
source("R/array_profiles.R")
source("R/cpg_annotation.R")

# --- 1. Array probes and hg38 positions --------------------------------------

array_type <- unique(read.csv(SAMPLE_SHEET_QC)$Array)
stopifnot(length(array_type) == 1, array_type %in% names(ARRAY_PROFILES))
profile <- ARRAY_PROFILES[[array_type]]
library(profile$anno_pkg, character.only = TRUE)
ann <- getAnnotation(get(profile$anno_pkg))
cat("Array:", array_type, "-", nrow(ann), "probes\n")

cpg_pos <- cpg_positions_hg38(rownames(ann), profile, CHAIN_HG19_HG38)

# --- 2. Gene annotation ------------------------------------------------------

genes <- load_gencode_genes(GENCODE_GTF, GENE_TYPES)
cpg_annotation <- annotate_cpgs(cpg_pos, genes, PROMOTER_UPSTREAM, PROMOTER_DOWNSTREAM)

# Position and CpG island context (EPICv1 splits shores and shelves into
# north/south; merged here to match EPICv2)
cpg_annotation <- cbind(
  data.table(cpg = cpg_pos$cpg, chr = cpg_pos$chr, pos_hg38 = cpg_pos$pos,
             island = sub("^[NS]_", "", ann[cpg_pos$cpg, "Relation_to_Island"])),
  cpg_annotation[, -"cpg"]
)

cat("\nCpGs by gene region:\n")
print(cpg_annotation[, .N, by = gene_region][order(-N)])
cat("CpGs by island context:\n")
print(cpg_annotation[, .N, by = island][order(-N)])
cat("Median distance of intergenic CpGs to the nearest gene:",
    median(cpg_annotation[gene_region == "intergenic", distance_to_gene]), "bp\n")

# --- 3. Save -----------------------------------------------------------------

fwrite(cpg_annotation, CPG_ANNOTATION)
cat("\nCpG annotation saved to:", CPG_ANNOTATION, "\n")
cat("Genes:", basename(GENCODE_GTF), "| types:", paste(GENE_TYPES, collapse = ", "),
    "| promoter:", PROMOTER_UPSTREAM, "bp up /", PROMOTER_DOWNSTREAM, "bp down of the TSS\n")
