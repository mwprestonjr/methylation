# =============================================================================
# CpG annotation
# Author: GP2 Subtypes and Mechanisms - M.P.
# Date: Oct 8, 2026
# Description: hg38 positions of array probes and their gene annotation from
#              GENCODE (GENCODE_GTF in config.R). Shared by the modules: the
#              annotation table is built once per data source by
#              preprocessing/4_cpg_annotation.R (CPG_ANNOTATION), and the hg38
#              positions are also used by the mQTL scan. Needs
#              R/array_profiles.R sourced, and minfi, data.table,
#              GenomicRanges and rtracklayer loaded
# =============================================================================

# hg38 positions of the given CpGs (probe names) from the array's annotation;
# probes on hg19 arrays (EPICv1) are lifted over with CHAIN_HG19_HG38. Probes
# that don't map to exactly one hg38 position on the same chromosome are
# dropped. Returns a data.frame: cpg, chr (without "chr", e.g. "19"), pos
cpg_positions_hg38 <- function(cpgs, profile, chain_hg19_hg38) {
  library(profile$anno_pkg, character.only = TRUE)
  ann <- getAnnotation(get(profile$anno_pkg))
  ann <- ann[cpgs, c("chr", "pos")]
  gr  <- GRanges(ann$chr, IRanges(ann$pos, width = 1), cpg = rownames(ann))

  if (profile$genome == "hg38") {
    gr_hg38 <- gr
  } else if (profile$genome == "hg19") {
    cat("Lifting probe coordinates hg19 -> hg38...\n")
    # import.chain() needs an uncompressed file
    chain_file <- tempfile(fileext = ".chain")
    system2("gunzip", c("-c", chain_hg19_hg38), stdout = chain_file)
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
    cat("Probes lifted:", length(gr_hg38), "of", length(cpgs),
        "(dropped", length(cpgs) - length(gr_hg38), ")\n")
  } else {
    stop("No hg38 conversion for genome build ", profile$genome)
  }

  data.frame(cpg = gr_hg38$cpg,
             chr = sub("^chr", "", as.character(seqnames(gr_hg38))),
             pos = start(gr_hg38))
}

# Genes from a GENCODE GTF (only the gene records are read), restricted to
# gene_types. Returns a stranded GRanges with gene_id, gene_name, gene_type
load_gencode_genes <- function(gtf, gene_types) {
  cat("Reading genes from", basename(gtf), "...\n")
  g <- fread(cmd = paste("zcat", shQuote(gtf), "| awk -F'\\t' '$3 == \"gene\"'"),
             sep = "\t", header = FALSE, quote = "",
             col.names = c("chr", "source", "feature", "start", "end",
                           "score", "strand", "frame", "attributes"))
  attr_value <- function(key) sub(paste0('.*', key, ' "([^"]*)".*'), "\\1", g$attributes)
  g[, `:=`(gene_id   = attr_value("gene_id"),
           gene_name = attr_value("gene_name"),
           gene_type = attr_value("gene_type"))]
  g <- g[gene_type %in% gene_types]
  cat("Genes kept:", nrow(g), "(", paste(gene_types, collapse = ", "), ")\n")
  GRanges(sub("^chr", "", g$chr), IRanges(g$start, g$end), strand = g$strand,
          gene_id = g$gene_id, gene_name = g$gene_name, gene_type = g$gene_type)
}

# Gene annotation of CpGs (cpg_pos from cpg_positions_hg38()). Each CpG gets
# one gene, chosen in this order:
#   promoter   within promoter_up bp upstream to promoter_down bp downstream of
#              a gene's transcription start (TSS)
#   gene_body  inside a gene
#   intergenic the nearest gene
# Within a category, protein-coding genes come before others, then the gene
# whose TSS is closest. Columns: gene, gene_id, gene_type, gene_region,
# distance_to_gene (0 if inside the gene or its promoter), distance_to_tss
# (CpG - TSS, in the gene's direction: negative = upstream)
annotate_cpgs <- function(cpg_pos, genes, promoter_up, promoter_down) {
  cpgs <- GRanges(cpg_pos$chr, IRanges(cpg_pos$pos, width = 1))
  tss  <- resize(genes, width = 1, fix = "start")
  prom <- promoters(genes, upstream = promoter_up, downstream = promoter_down)

  # Candidate CpG-gene pairs for each category
  hits <- function(ov, region) {
    data.table(i = queryHits(ov), g = subjectHits(ov), region = region)
  }
  near <- distanceToNearest(cpgs, genes, ignore.strand = TRUE)
  cand <- rbindlist(list(
    hits(findOverlaps(cpgs, prom, ignore.strand = TRUE), "promoter"),
    hits(findOverlaps(cpgs, genes, ignore.strand = TRUE), "gene_body"),
    hits(near, "intergenic")
  ))
  cand[, rank_region := match(region, c("promoter", "gene_body", "intergenic"))]
  cand[, non_coding := genes$gene_type[g] != "protein_coding"]   # FALSE (protein-coding) sorts first
  cand[, tss_dist := abs(start(tss)[g] - start(cpgs)[i])]
  setorder(cand, i, rank_region, non_coding, tss_dist)
  best <- cand[!duplicated(i)]

  # Signed distance to the TSS in the gene's direction
  strand_sign <- ifelse(as.character(strand(genes))[best$g] == "-", -1, 1)
  out <- data.table(cpg = cpg_pos$cpg)
  out[best$i, `:=`(
    gene             = genes$gene_name[best$g],
    gene_id          = genes$gene_id[best$g],
    gene_type        = genes$gene_type[best$g],
    gene_region      = best$region,
    distance_to_gene = ifelse(best$region == "intergenic", mcols(near)$distance[match(best$i, queryHits(near))], 0L),
    distance_to_tss  = (start(cpgs)[best$i] - start(tss)[best$g]) * strand_sign
  )]
  out[]
}
