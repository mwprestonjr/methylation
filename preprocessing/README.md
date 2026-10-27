# Preprocessing (`preprocessing/`)

Estimates blood cell proportions, checks sources of variation, applies ComBat batch correction, and annotates every CpG (hg38 position, gene, island context). Reads the `qc` outputs.

Run from the repository root with the `methylation` conda environment active
(see [Running the pipeline](../README.md#running-the-pipeline)):

```bash
make preprocessing SOURCE=<source>
Rscript preprocessing/2_variation_sources.R ppmi_p140 combat   # a single script
```

Settings: [`config.R`](config.R) in this folder; shared paths and settings in
the root [`config.R`](../config.R).

## Scripts

- `1_cell_counts.R` - estimates blood cell proportions.
- `2_variation_sources.R` - tests sources of variation before ComBat to choose `COMBAT_BATCH_VAR`, then reruns after ComBat to check the correction.
- `3_combat.R` - applies ComBat batch correction.
- `4_cpg_annotation.R` - creates `CPG_ANNOTATION`, an annotation for every array probe with hg38 position (liftover for EPICv1), GENCODE gene and gene type, promoter/gene-body/intergenic context (nearest gene for intergenic probes), distances to the gene and its TSS, and CpG island context.
- `config.R` - variation-analysis and ComBat settings.

## Outputs

In `<DIR_OUTPUT>/results/preprocessing/`. The config name is the variable
other scripts use to find the file (root `config.R`, or
`preprocessing/config.R` for files used only within preprocessing).

| Step | File | Config name | Contents | Read by |
|---|---|---|---|---|
| 1 | `cell_proportions.csv` | `CELL_PROPORTIONS` | Estimated blood cell proportions per sample (`CELL_TYPES`) | `2_variation_sources.R`, `3_combat.R` |
| 2 | `variation_sources_<raw\|combat>.rds` | | PCs, variance explained and variable-PC association p-values behind the variation figures | |
| 2 | `variation_svd_pvalues_<raw\|combat>.csv` | | p-values for each sample variable vs each of the top PCs | |
| 2 | `confounder_summary_<raw\|combat>.csv` | | Whether age, sex, batch and chip differ between phenotypes | |
| 3 | `combat_mVals.rds` | `COMBAT_MVALS` | Batch-corrected M values (`COMBAT_BATCH_VAR`, protecting `COMBAT_PROTECT`) | mQTL (`R/mqtl_functions.R`, `4_mqtl_plots.R`) |
| 3 | `sample_sheet_final.csv` | `SAMPLE_SHEET_FINAL` | QC-passed sample sheet + cell proportions + `combat_batch` (the batch ComBat used) | clocks, mQTL |
| 3 | `combat_summary.csv` | | Samples, CpGs, batch variable, number and smallest size of batches, protected variables | |
| 4 | `cpg_annotation.csv` | `CPG_ANNOTATION` | Every probe on the array: hg38 position, island context, GENCODE gene and type, promoter/body/intergenic, distances to gene and TSS | `mqtl/4_mqtl_plots.R` |

## Methods and decisions

**Cell proportions** (`1_cell_counts.R`): CD8T, CD4T, NK, B cells, monocytes
and neutrophils (`CELL_TYPES`), estimated with the IDOL reference
(`FlowSorted.Blood.EPIC`, Salas et al. 2018), from the noob-normalized betas
saved by `qc/2_qc.R` (`BVALS_NOOB`), since the method needs noob rather than
Funnorm values. minfi's `estimateCellCounts2` only supports 450K/EPICv1, so its
IDOL steps are run directly, the same way for both arrays; EPICv2 replicate
probes are averaged per CpG first. In later models, neutrophils are left out because the
proportions sum to about 1.

**EPICv2 coverage of the cell reference**: the IDOL CpGs were selected for
EPICv1, and 30 of the 450 aren't on EPICv2, so Psomagen's cell proportions are
estimated from 420 CpGs (PPMI: all 450). `1_cell_counts.R` prints the number
found. The estimates still work, but they may be somewhat less precise on
EPICv2, and cell proportions aren't directly comparable between the two
arrays.

**Sources of variation** (`2_variation_sources.R`): SVD of the
`VARIATION_TOP_CPGS` (50,000) most variable CpGs; each of the first
`VARIATION_N_PCS` (10) components is tested against each sample variable
(Kruskal-Wallis for categorical variables, linear regression for numeric ones,
as in ChAMP's `champ.SVD`). It also tests whether age, sex, batch and chip
differ between phenotypes: a batch variable confounded with phenotype can't be
corrected without removing phenotype signal. Run before ComBat (`raw`) to
choose the batch variable, and after (`combat`) to check the correction.

**Batch correction** (`3_combat.R`): ComBat (`sva`) on the probe-filtered M
values (one row per site: EPICv2 replicate probes are averaged in QC), with
parametric priors.
- **Batch variable** (`COMBAT_BATCH_VAR`): `Batch` (processing plate) for
  PPMI and `Sentrix_ID` (chip) for Psomagen. PPMI samples are spread over about
  240 chips, often one per chip, too few for ComBat (a batch of one sample
  would have its own signal removed); its 22 plates hold 4-40 samples each.
  Psomagen has about 8 samples per chip, and chip explains more technical
  variation than delivery. How `Batch` is defined: see
  [qc/README.md](../qc/README.md#sample-sheet).
- **Protected variables** (`COMBAT_PROTECT`): phenotype, so ComBat keeps
  phenotype differences. If batch and phenotype overlap strongly, protecting
  phenotype can inflate phenotype differences (Nygaard et al. 2016); check
  `confounder_summary_raw.csv`.
- **Not corrected**: position on the chip (row R01-R08). Later analyses adjust
  for it as a covariate (`chip_row`).
- **Not ComBat-corrected**: the unfiltered betas used by the clocks; the clock
  models adjust for the ComBat batch (`combat_batch` in the final sample
  sheet) instead.

**CpG annotation** (`4_cpg_annotation.R`), for every probe on the array,
independent of QC:
- **Position**: hg38. EPICv1 positions (hg19 annotation) are lifted over with
  the UCSC chain; probes that don't map to exactly one hg38 position on the
  same chromosome are dropped.
- **Gene**: one GENCODE gene per CpG (`GENE_TYPES`: protein-coding and lncRNA),
  chosen in this order: promoter (`PROMOTER_UPSTREAM` 1,500 bp upstream to
  `PROMOTER_DOWNSTREAM` 500 bp downstream of the transcription start), gene
  body, otherwise the nearest gene (intergenic). Within a category,
  protein-coding genes come first, then the gene with the closest start.
- **Island context**: from the Illumina annotation (`Relation_to_Island`);
  EPICv1's north/south shores and shelves are merged to match EPICv2.

## Figures

Saved in `<DIR_OUTPUT>/figures/preprocessing/` for the data source being processed
(`DIR_OUTPUT` in the root `config.R`).

- `cell_proportions.png` - estimated blood-cell proportions across samples.
- `cell_proportions_by_phenotype.png` - estimated blood-cell proportions by phenotype.
- `variation_<raw|combat>_svd_heatmap.png` - associations between sample variables and SVD components, before (`raw`) or after (`combat`) ComBat.
- `variation_<raw|combat>_scree.png` - variance explained by SVD components, before or after ComBat.
- `variation_<raw|combat>_pca_by_<variable>.png` - PC1 vs PC2, colored by each available variable among `Chip`, `Batch`, `phenotype`, `sex`, and `age`; generated only for variables with usable variation in the data.
