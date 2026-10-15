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

## Outputs

In `<DIR_OUTPUT>/results/qc/`. The config name is the variable other scripts
use to find the file (root `config.R`, or `qc/config.R` for files used only
within qc).

| Step | File | Config name | Contents | Read by |
|---|---|---|---|---|
| 1 | `sample_sheet.csv` | `SAMPLE_SHEET` | One row per array sample: `SAMPLE_SHEET_COLUMNS` (IDs, phenotype, sex, race, age, batch, chip/position, array, idat path) | `2_qc.R` |
| 2 | `qc_sample_metrics.csv` | | Every input sample: median intensities, mean detection p, predicted vs reported sex, density score, and why it was removed | |
| 2 | `qc_sesame_stats.csv` | | SeSAMe per-sample stats (detection, intensity, dye bias); only when `USE_SESAME_QC = TRUE` | |
| 2 | `qc_summary.csv` | | Number of samples and probes removed at each QC step | |
| 2 | `qc_figure_data.rds` | `QC_FIGURE_DATA` | Sample metrics, density curves and thresholds behind the QC figures | `3_qc_plots.R` |
| 2 | `sample_sheet_qc_passed.csv` | `SAMPLE_SHEET_QC` | Sample sheet of the samples that passed QC | preprocessing (all steps) |
| 2 | `mVals.rds` | `MVALS_QC` | M values: QC-passed samples, Funnorm-normalized, probe-filtered, infinite values capped | `preprocessing/2_variation_sources.R`, `3_combat.R` |
| 2 | `bVals.rds` | | Beta values, same samples and probes as `mVals.rds` | not read by the pipeline |
| 2 | `bVals_unfiltered.rds` | `BVALS_UNFILTERED` | Beta values: QC-passed samples, normalized, **no probe filtering** (filters can remove clock CpGs) | `clocks/1_estimate.R` |
| 2 | `mSetSq_qc_passed.rds` | `MSET_QC` | Normalized, filtered `GenomicRatioSet` (minfi object with annotation) | not read by the pipeline |

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
