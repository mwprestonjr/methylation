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

## Batch variable

`COMBAT_BATCH_VAR` (`config.R`) is `Batch` (processing plate) for PPMI and
`Sentrix_ID` (chip) for Psomagen, chosen with `2_variation_sources.R`. How
`Batch` is defined: see [qc/README.md](../qc/README.md#sample-sheet).

## Figures

Saved in `<DIR_OUTPUT>/figures/preprocessing/` for the data source being processed
(`DIR_OUTPUT` in the root `config.R`).

- `cell_proportions.png` - estimated blood-cell proportions across samples.
- `cell_proportions_by_phenotype.png` - estimated blood-cell proportions by phenotype.
- `variation_<raw|combat>_svd_heatmap.png` - associations between sample variables and SVD components, before (`raw`) or after (`combat`) ComBat.
- `variation_<raw|combat>_scree.png` - variance explained by SVD components, before or after ComBat.
- `variation_<raw|combat>_pca_by_<variable>.png` - PC1 vs PC2, colored by each available variable among `Chip`, `Batch`, `phenotype`, `sex`, and `age`; generated only for variables with usable variation in the data.
