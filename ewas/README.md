# Epigenome-wide association studies (`ewas/`)

Tests every CpG for association with diagnosis (case-control), age, and an
exposure (smoking), using limma on the ComBat-corrected M values. Reads the
`preprocessing` outputs and, for smoking, the GP2 extended clinical data. Each data source is run separately.

Run from the repository root with the `methylation` conda environment active
(see [Running the pipeline](../README.md#running-the-pipeline)):

```bash
make ewas SOURCE=<source>
Rscript ewas/2_age.R ppmi_p140   # a single script
```

Settings: [`config.R`](config.R) in this folder; shared paths and settings in
the root [`config.R`](../config.R).

## Inputs

From `preprocessing` (paths in the root `config.R`), so run
`make preprocessing` first, plus the GP2 clinical data for smoking:

| File | Name in code | From | Used for |
|---|---|---|---|
| `results/preprocessing/combat_mVals.rds` | `COMBAT_MVALS` | `preprocessing/3_combat.R` | M values, CpGs x samples (columns are idat basenames) |
| `results/preprocessing/sample_sheet_final.csv` | `SAMPLE_SHEET_FINAL` | `preprocessing/3_combat.R` | phenotype, age, sex, cell proportions, `Sentrix_Position` (chip row), any exposure column |
| `results/preprocessing/cpg_annotation.csv` | `CPG_ANNOTATION` | `preprocessing/4_cpg_annotation.R` | hg38 position, island context and gene, joined to the results by `cpg` |
| `r12_extended_clinical_data.csv` | `CLINICAL_DATA` | GP2 R12, `gp2_release12` bucket (root `config.R`) | `smoking_status` at the blood draw, joined by GP2ID (`3_exposure.R` only) |

## Steps

| Step | What it does | Outputs (`results/ewas/`) |
|---|---|---|
| `1_case_control.R` | `EWAS_GROUPS[2]` vs `EWAS_GROUPS[1]` (PD vs Control), adjusted for `EWAS_COVARIATES`. Skipped if a group has fewer than `EWAS_MIN_GROUP_N` samples | `case_control_*` |
| `2_age.R` | Age at the blood draw (per year), adjusted for the other `EWAS_COVARIATES` and phenotype | `age_*` |
| `3_exposure.R` | Smoking (`EWAS_EXPOSURE_COL`): ever vs never for PPMI, former and current vs never for Psomagen, adjusted for `EWAS_COVARIATES` and phenotype. Any other continuous or categorical sample sheet column can be set instead. Skipped if the exposure isn't available | `smoking_values.csv`, `exposure_<column>_<level>_*` |
| `4_plots.R` | Manhattan, QQ and distribution plots for every model with results; combines the model summaries | `ewas_model_summary.csv`, `figures/ewas/*.png` |

Each model (prefix `case_control`, `age`, `exposure_...`) writes:

| File | Contents |
|---|---|
| `<prefix>_results.csv.gz` | Every CpG, sorted by p: `logFC` (effect on the M value scale), `AveExpr` (mean M value), `t`, `P.Value`, `FDR`, `B`, plus the `CPG_ANNOTATION` columns |
| `<prefix>_top_hits.csv` | The `EWAS_TOP_HITS` (100) CpGs with the smallest p, same columns |
| `<prefix>_summary.csv` | One row: label, coefficient, samples (and per group), covariates used, CpGs tested, CpGs with FDR < `EWAS_FDR`, smallest p, genomic inflation (lambda) |

A model that is skipped deletes its outputs and figures from earlier runs, so
`4_plots.R` and `ewas_model_summary.csv` only show models from the current
settings. `3_exposure.R` always deletes all earlier `exposure_*` outputs, since
the exposure column or its levels may have changed.

## Methods and decisions

- **Model**: for each CpG, `M value ~ variable + covariates`, fitted with
  limma (`lmFit`, then `eBayes` moderated t statistics). Effects (`logFC`) are
  differences in M value (log2 ratio of methylated to unmethylated signal):
  per year for age, per unit for a continuous exposure, vs the reference
  group for categorical variables.
- **Covariates** (`EWAS_COVARIATES`): age, sex, cell proportions (without
  neutrophils, as the proportions sum to ~1) and chip row. The age model
  leaves age out because age is what it tests. Categorical covariates with
  one value among the model's samples are left out (reported in the log).
- **Batch**: the M values are already ComBat-corrected for the batch variable
  (`COMBAT_BATCH_VAR`: plate for PPMI, chip for Psomagen), so batch is not a
  covariate. Adjusting for it again would remove little, and would spend one
  degree of freedom per batch (hundreds of chips for Psomagen). This differs
  from `clocks/`, whose betas aren't ComBat-corrected.
- **Phenotype adjustment**: the age and exposure models also adjust for
  phenotype (`EWAS_ADJUST_PHENOTYPE`), so a PD-control difference that is
  correlated with age or exposure isn't reported as an age or exposure effect.
  Phenotype is left out automatically when all samples have the same one
  (e.g. Psomagen).
- **ComBat protection**: ComBat kept the variation explained by
  `COMBAT_PROTECT` (phenotype), so case-control differences survive batch
  correction. Age and exposures aren't protected: if they differ between
  batches, part of their effect was removed with the batch. Add them to
  `COMBAT_PROTECT` (`preprocessing/config.R`) if this matters.
- **Samples**: samples with a missing value for the variable or any
  covariate are dropped from that model. The case-control model uses only
  samples in `EWAS_GROUPS`. Other phenotypes (e.g. Prodromal) are still in
  the age and exposure models.
- **Smoking** (`3_exposure.R`, `load_smoking()` in `R/ewas_functions.R`):
  the harmonized `smoking_status` of the GP2 extended clinical data
  (`CLINICAL_DATA`), at the first of `EWAS_CLINICAL_VISIT_MONTHS` (baseline,
  then screening) with a known value, joined by GP2ID. In R12 every
  participant with a known status has it at baseline, and it's the same at
  every visit. Two versions are built: `smoking_status` (Never, Former,
  Current) and `smoking_ever` (Never, Ever); "Unknown" or no clinical record
  is missing, and those samples drop out of the model. The values used are
  saved per sample in `smoking_values.csv`.
  - **PPMI** tests `smoking_ever`: PPMI's "Former smoker" means "smoked 100
    or more cigarettes" (PD Risk Factor Questionnaire, `pdrfuq_smoke_100`),
    with no current/former distinction; no PPMI participant is "Current
    smoker". Former and current smokers are therefore mixed in "Ever".
  - **Psomagen** tests `smoking_status` (Former and Current, each vs Never),
    for the studies that record it.
  - Pack-years and smoking years aren't used: few studies record them, and
    in different units (PREDICTPD median ~160 pack-years vs ~10 elsewhere).
    `smoking_ever_smoker_yes` isn't used either: it contradicts
    `smoking_status` for some participants.
- **Exposure model**: for a categorical exposure, all levels are in one
  model; levels with fewer than `EWAS_MIN_GROUP_N` samples are left out (NA),
  and each other level is reported vs the reference in its own result set.
  Empty strings count as missing.
- **Multiple testing**: Benjamini-Hochberg FDR across all CpGs within each
  model (and each exposure level). CpGs with FDR < `EWAS_FDR` are counted as
  significant.
- **Genomic inflation**: lambda is the median chi-square statistic from the
  p-values divided by its expected value under the null (0.456), shown in the
  summary and on the QQ plot. Values well above 1 suggest uncorrected
  confounding (e.g. batch, cell composition) or a highly polygenic signal;
  age EWAS typically shows real, widespread inflation.

## Figures

Saved in `<DIR_OUTPUT>/figures/ewas/`, one set per model:

- `<prefix>_manhattan.png` - -log10(p) along the genome (hg38). Dashed line at
  `EWAS_PLOT_SUGGESTIVE_P` (1e-5), solid line at the FDR < `EWAS_FDR`
  threshold (shown when any CpG passes it).
- `<prefix>_qq.png` - observed vs expected -log10(p), with lambda.
- `<prefix>_distributions.png` - histograms of the effects and of -log10(p).
