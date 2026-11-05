# Overview

Describes the samples available in each data source, side by side. These scripts
cover all data sources at once, so they take no arguments, and they aren't part
of `make`.

## Running

From the repository root, on the machine with the GP2 data mounted:

```bash
Rscript overview/1_sample_overview.R
Rscript overview/check_trait_data.R
```

`1_sample_overview.R` needs `make qc` to have finished for every data source,
plus `GP2_MASTER_KEY` and `PPMI_GENETICS`. `check_trait_data.R` needs `make
preprocessing` for every data source, plus `GP2_MASTER_KEY` and `CLINICAL_DATA`.

## Scripts

| Script | Does |
|---|---|
| `1_sample_overview.R` | Sample counts before and after QC, Table 1 and figures |
| `check_trait_data.R` | Samples with a value of each clinical trait and exposure, to decide which have enough data to run an EWAS. Its settings are in the script, and it also reads `ewas/config.R` |

## Outputs

In `/mnt/output/methylation/reports/overview/` (`DIR_OVERVIEW` in
`overview/config.R`), except the trait check:

| File | Config name | Contents |
|---|---|---|
| `overview_counts.csv` | `OVERVIEW_COUNTS` | One row per data source x stage x phenotype x variable x level: `n` and `pct`. `phenotype` is `All` or a phenotype label |
| `overview_table1.csv` | `OVERVIEW_TABLE1` | Summary of the samples after QC: one column per data source x phenotype. `-` means none in that group, or not known for that data source |
| `overview_01_phenotype.png` | - | Samples by phenotype |
| `overview_02_sex.png` | - | Share of each sex by phenotype, with counts |
| `overview_03_ancestry.png` | - | Samples by genetic ancestry |
| `overview_04_age.png` | - | Age at sample collection by phenotype |
| `../ewas_trait_data_check.csv` | `EWAS_TRAIT_CHECK_FILE` | Trait check, counts only |

Nothing here is read by another module.

## Methods and decisions

- **Stages**: `before_qc` is `SAMPLE_SHEET` and `after_qc` is
  `SAMPLE_SHEET_QC`. `SAMPLE_SHEET_FINAL` isn't counted because preprocessing
  removes no samples. `before_qc` isn't everything delivered: each sample sheet
  script has already selected the samples of interest (for PPMI, baseline PD
  and control samples).
- **Samples and participants**: a participant can have more than one array
  before QC (replicates), so both are counted. Percentages are of samples.
- **Phenotype**: labels are shown as they are in the sample sheet rather than
  as case vs control. Psomagen has no controls, and its non-PD samples (DLB,
  MSA, PSP, Other) aren't controls.
- **Ancestry and race**: ancestry is the GP2 genetic label (`GP2_ANCESTRY_COL`
  in `GP2_MASTER_KEY`), joined by GP2ID. Race is self-reported, from the sample
  sheet, and is a coarser measure: for example, it doesn't separate Ashkenazi
  Jewish participants in PPMI. Labels aren't pooled, unlike the EWAS covariate
  (`EWAS_MIN_ANCESTRY_N`), so these counts can differ from that covariate's.
- **Missing values**: a missing phenotype, sex, race or ancestry is counted as
  `Unknown`. Monogenic variant counts appear only for data sources with variant
  data (`PPMI_GENETICS`); participants missing from that file are `No variant
  data`.
- **Batches**: distinct values of the sample sheet's `Batch` column: plates for
  PPMI and deliveries for Psomagen. This isn't the same as `combat_batch` for
  Psomagen, which uses chip.
