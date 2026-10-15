# QC (`qc/`)

Builds the sample sheet for a data source, runs sample and probe QC, normalizes (Funnorm), and plots the QC results.

Run from the repository root with the `methylation` conda environment active
(see [Running the pipeline](../README.md#running-the-pipeline)):

```bash
make qc SOURCE=<source>
Rscript qc/2_qc.R ppmi_p140   # a single script
```

Settings: [`config.R`](config.R) in this folder; shared paths and settings in
the root [`config.R`](../config.R).

## Scripts

- `1_build_sample_sheet_ppmi.R` - builds the PPMI sample sheet.
- `1_build_sample_sheet_psomagen.R` - builds the Psomagen sample sheet. Both sample-sheet scripts write the common `SAMPLE_SHEET_COLUMNS` from the root `config.R`, so downstream scripts work independently of data source.
- `2_qc.R` - performs sample and probe QC and normalization. It supports any array in `R/array_profiles.R`, with one array type processed per run.
- `3_qc_plots.R` - generates QC figures from the figure data saved by `2_qc.R`.
- `config.R` - QC thresholds and method settings.

## Sample sheet

Both sample sheet scripts write the same columns (`SAMPLE_SHEET_COLUMNS` in the
root `config.R`), so later scripts don't depend on the data source.

- `Batch`: for Psomagen, the delivery (dataset) ID; for PPMI, the top-level
  idat folder in the `Basename` path, i.e. the processing plate, e.g.
  `.../idat/20190604_plate1/SENTRIXID/...` -> `Batch = "20190604_plate1"`

## Figures

Saved in `<DIR_OUTPUT>/figures/qc/` for the data source being processed
(`DIR_OUTPUT` in the root `config.R`).

- `qc_01_median_intensities.png` - methylated vs unmethylated median intensity per sample.
- `qc_02_detection_pvalues.png` - mean detection p-value per sample and the QC threshold.
- `qc_02b_sesame_detection.png` - fraction of detected CpG probes per sample; generated only when SeSAMe QC is enabled.
- `qc_03_sex_prediction.png` - predicted vs reported sex from X/Y intensities.
- `qc_04_density_before_normalization.png` - beta-value densities before normalization, by phenotype.
- `qc_04b_density_by_dataset.png` - beta-value densities before normalization, by dataset.
- `qc_05_density_after_normalization.png` - beta-value densities after normalization, by phenotype.

## Troubleshooting

### `preprocessFunnorm` pthread_create error on Linux/GCP VMs

If normalization in `2_qc.R` fails with a `pthread_create` error, reinstall `preprocessCore` with threading disabled:
```r
BiocManager::install("preprocessCore", 
                     configure.args = "--disable-threading",
                     force = TRUE)
```
