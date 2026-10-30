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
`make preprocessing` first, plus the GP2 master key for ancestry, the GP2
clinical data for smoking and, for PPMI, the monogenic PD variants:

| File | Name in code | From | Used for |
|---|---|---|---|
| `results/preprocessing/combat_mVals.rds` | `COMBAT_MVALS` | `preprocessing/3_combat.R` | M values, CpGs x samples (columns are idat basenames) |
| `results/preprocessing/sample_sheet_final.csv` | `SAMPLE_SHEET_FINAL` | `preprocessing/3_combat.R` | phenotype, age, sex, cell proportions, `Sentrix_Position` (chip row), any exposure column |
| `results/preprocessing/cpg_annotation.csv` | `CPG_ANNOTATION` | `preprocessing/4_cpg_annotation.R` | hg38 position, island context and gene, joined to the results by `cpg` |
| GP2 master key | `GP2_MASTER_KEY` | GP2 R12, `gp2_release12` bucket (root `config.R`) | Genetic ancestry (`EWAS_ANCESTRY_COL`), joined by GP2ID (all models) |
| `r12_extended_clinical_data.csv` | `CLINICAL_DATA` | GP2 R12, `gp2_release12` bucket (root `config.R`) | Smoking questionnaire items (`pdrfuq_smoke_*` for PPMI, `merq_cigarettes*` for Nigeria-PD), joined by GP2ID (`3_exposure.R` only) |
| `iu_genetic_consensus_*.csv` | `PPMI_GENETICS` | PPMI (LONI), root `config.R`; `NULL` for Psomagen | Monogenic PD variant per participant (`pd_gene`, `load_pd_gene()` in `R/pd_genetics.R`), joined by `clinical_id` (`1_case_control.R` only) |

## Steps

| Step | What it does | Outputs (`results/ewas/`) |
|---|---|---|
| `1_case_control.R` | `EWAS_GROUPS[2]` vs `EWAS_GROUPS[1]` (PD vs Control), adjusted for `EWAS_COVARIATES`. Where monogenic variants are known (PPMI), also PD vs Control without carriers, and carriers of each gene in `EWAS_CARRIER_GENES` (LRRK2, GBA1) vs non-carrier PD. Each model is skipped if a group has fewer than `EWAS_MIN_GROUP_N` samples (`EWAS_MIN_CARRIER_N`, 10, for the carrier models, so PPMI's 8 GBA1 carriers are skipped) | `case_control_*`, `case_control_no_carriers_*`, `carrier_<gene>_*` |
| `2_age.R` | Age at the blood draw (per year), adjusted for the other `EWAS_COVARIATES` and phenotype | `age_*` |
| `3_exposure.R` | Smoking (`EWAS_EXPOSURE_COL`): ever vs never regular smoking for both data sources, adjusted for `EWAS_COVARIATES` and phenotype. Any other continuous or categorical sample sheet column can be set instead. Skipped if the exposure isn't available | `smoking_values.csv`, `exposure_<column>_<level>_*` |
| `4_plots.R` | Manhattan, QQ and distribution plots for every model with results; combines the model summaries | `ewas_model_summary.csv`, `figures/ewas/*.png` |

Each model (prefix `case_control`, `case_control_no_carriers`, `carrier_<gene>`, `age`, `exposure_...`) writes:

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
  neutrophils, as the proportions sum to ~1), chip row and genetic ancestry.
  The age model leaves age out because age is what it tests. Categorical
  covariates with one value among the model's samples are left out (reported
  in the log).
- **Ancestry** (`add_ancestry()` in `R/ewas_functions.R`): the GP2 NBA
  ancestry label (`EWAS_ANCESTRY_COL`, `nba_label` in `GP2_MASTER_KEY`),
  joined by GP2ID. Methylation differs between ancestries, and in PPMI
  ancestry is unbalanced between groups: 58 of the 65 Ashkenazi Jewish (AJ)
  samples are cases, and most LRRK2 carriers (65 of 88) are AJ or another
  non-EUR group. Samples without a label (not genotyped or not in the
  master key; PPMI 24, Psomagen 8) are "Unknown" rather than dropped, and
  labels with fewer than `EWAS_MIN_ANCESTRY_N` (10) samples are pooled into
  "Other". PPMI's model then has EUR, AJ, Other and Unknown. Psomagen is
  nearly all AFR, so the adjustment does little there. A label is coarse:
  it doesn't capture variation within a group (genetic PCs would, but only
  for genotyped samples).
- **Batch**: the M values are already ComBat-corrected for the batch variable
  (`COMBAT_BATCH_VAR`: plate for PPMI, chip for Psomagen), so batch is not a
  covariate. Adjusting for it again would remove little, and would spend one
  degree of freedom per batch (hundreds of chips for Psomagen). This differs
  from `clocks/`, whose betas aren't ComBat-corrected.
- **Phenotype adjustment**: the age and exposure models also adjust for
  phenotype (`EWAS_ADJUST_PHENOTYPE`), so a PD-control difference that is
  correlated with age or exposure isn't reported as an age or exposure effect.
  Phenotype is left out automatically when all samples have the same one.
  Psomagen keeps it: besides PD it has a few DLB, MSA, PSP and Other samples.
- **ComBat protection**: ComBat kept the variation explained by
  `COMBAT_PROTECT` (phenotype), so case-control differences survive batch
  correction. Age and exposures aren't protected: if they differ between
  batches, part of their effect was removed with the batch. Add them to
  `COMBAT_PROTECT` (`preprocessing/config.R`) if this matters.
- **Samples**: samples with a missing value for the variable or any
  covariate are dropped from that model. The case-control model uses only
  samples in `EWAS_GROUPS`. Other phenotypes (e.g. Prodromal) are still in
  the age and exposure models.
- **Monogenic PD carriers** (`1_case_control.R`, PPMI only): PPMI enrolled
  carriers of known PD variants on purpose (genetic cohort), so about a third
  of the PPMI cases carry one (LRRK2 88, GBA1 8, LRRK2+GBA1 5, PRKN 6, no
  variant 206). The main case-control model includes them. Two further
  model types separate them out, using `pd_gene` from `PPMI_GENETICS`
  (`load_pd_gene()` in `R/pd_genetics.R`):
  - `case_control_no_carriers`: PD vs Control without carriers of any known
    variant, in either group (`pd_gene` other than `EWAS_NONCARRIER`).
    Samples without variant data are kept.
  - `carrier_<gene>`: within PD, carriers of one gene in
    `EWAS_CARRIER_GENES` vs PD without a known variant. Carriers of another
    gene or of two (LRRK2+GBA1) are left out. Phenotype is constant, so it
    isn't a covariate.

  Caveats for the carrier models:
  - **Inherited DNA around the variant**: carriers share the stretch of
    chromosome around their variant (the founder haplotype), and with it any
    nearby SNP that affects methylation. CpGs near the gene can then differ
    for genetic, not disease, reasons. In a test run, the top LRRK2 CpG
    (cg14035730, ~125 kb from LRRK2, SLC2A13) has one of the strongest
    cis-mQTLs in PPMI (`mqtl/`, p ~ 1e-62), and the top GBA1 CpGs (LMNA,
    ~900 kb from GBA1) also have cis-mQTLs. Check hits near the gene against
    the mQTL results before reading them as carrier biology.
  - **Ancestry**: most PPMI LRRK2 carriers (65 of 88) aren't in the GP2 EUR
    group, mostly Ashkenazi Jewish (AJ), where LRRK2 G2019S is common.
    Ancestry is a covariate (above), but with this overlap the LRRK2 effect
    is estimated mostly from the 23 EUR carriers, and differences within
    the AJ group aren't adjusted for.
  - **Group size**: 8 GBA1 carriers is enough to fit
    (`EWAS_MIN_GROUP_N`) but not to estimate effects reliably; treat GBA1
    results as exploratory.
- **Smoking** (`3_exposure.R`, `load_smoking()` in `R/ewas_functions.R`):
  built from each study's own questionnaire items in the GP2 extended
  clinical data (`CLINICAL_DATA`), at the first of
  `EWAS_CLINICAL_VISIT_MONTHS` (baseline, then screening) with an answer,
  joined by GP2ID. Question wording is from the
  [GP2 data dictionary](https://github.com/GP2code/GP2-Data-Dictionary).
  Both data sources test `smoking_ever` (Never vs Ever), where **Ever means
  regular smoking: at least one cigarette a day for 6 months or longer**, the
  definition both questionnaires share:

  | | PPMI | Psomagen (Nigeria-PD) |
  |---|---|---|
  | Questionnaire | PD Risk Factor Questionnaire (PD RFQ-U) | Medical and Environmental Risk Questionnaire (MERQ-PD-B) |
  | Ever | `pdrfuq_smoke_regularly` = Yes: "In your lifetime, have you ever regularly smoked cigarettes, that is, at least one cigarette per day for 6 months or longer?" | `merq_cigarettes` = Yes: "Before you were diagnosed with PD, did you ever smoke cigarettes regularly, that is at least 1 cigarette per day for at least 6 months?", or `merq_cigarettes_current` = Yes |
  | Never | `pdrfuq_smoke_regularly` = No, or `pdrfuq_smoke_100` = No ("In your lifetime, have you smoked 100 or more cigarettes?"; the regular-smoking question is only asked after a yes) | `merq_cigarettes` = No and not `merq_cigarettes_current` = Yes |
  | Current smoking | Not asked | `merq_cigarettes_current`: "Do you smoke cigarettes regularly now?" |
  | Time frame | Lifetime | Before diagnosis (cases), plus now |

  How the definitions still differ:
  - **Current smoking** is only known for Psomagen, so `smoking_status`
    (Never, Former, Current) is only filled in there; PPMI's ever-smokers
    have no status. Psomagen's 5 current smokers are too few to test on their
    own, so both sources test ever vs never.
  - **Time frame**: MERQ asks about regular smoking *before diagnosis*.
    Someone who started after diagnosis answers no, but is counted as Ever
    through the current-smoking question. Someone who started after diagnosis
    and has since stopped would be counted as Never.
  - **Not used**: the harmonized `smoking_status` column. For PPMI it is
    derived from the looser 100-cigarettes question (about 10% of PPMI's
    100-cigarettes smokers never smoked regularly), and the Nigeria-PD MERQ
    answers aren't in it. `smoking_ever_smoker_yes` contradicts it for some
    participants. Pack-years and smoking years: few studies record them, and
    in different units (PREDICTPD median ~160 pack-years vs ~10 elsewhere).

  Answers other than Yes/No ("Don't Know", "Uncertain") count as missing, and
  samples without an answer drop out of the model. The values used, and the
  questionnaire they came from (`smoking_source`), are saved per sample in
  `smoking_values.csv`.
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
