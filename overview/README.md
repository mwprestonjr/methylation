# Overview

Describes the samples available in each data source, side by side, and how
the arrays they were run on compare. These scripts
cover all data sources at once, so they take no arguments, and they aren't part
of `make`.

## Running

From the repository root, on the machine with the GP2 data mounted:

```bash
Rscript overview/sample_overview.R
Rscript overview/array_comparison.R
Rscript overview/check_trait_data.R
```

`sample_overview.R` needs `make qc` to have finished for every data source,
plus `GP2_MASTER_KEY` and `PPMI_GENETICS`. `check_trait_data.R` needs `make
preprocessing` for every data source, plus `GP2_MASTER_KEY` and `CLINICAL_DATA`.
`array_comparison.R` needs `make preprocessing` for every data source, and
uses the mQTL results where `make mqtl` has run; it reads no GP2 inputs, so it
also runs in the dev container. It loads each data source's M values to get
their row names, so it takes about 10 minutes and a few GB of memory.

## Scripts

| Script | Does |
|---|---|
| `sample_overview.R` | Sample counts before and after QC, Table 1 and figures |
| `array_comparison.R` | CpG content of the two arrays (EPICv1, EPICv2): overlap, replicate probes, design and position differences, and overlap after QC and in the mQTL scan |
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
| `array_comparison_summary.csv` | `OVERVIEW_ARRAY_SUMMARY` | One row per metric, by section (Manifest, Replicate probes, Overlap, After QC, mQTL): `value` and a `note` |
| `array_comparison_context.csv` | `OVERVIEW_ARRAY_CONTEXT` | Sites on both arrays or one, by island context and by gene region: `n` and `pct` |
| `array_comparison_sites.csv.gz` | `OVERVIEW_ARRAY_SITES` | One row per CpG site on either array: whether it is on each array, number of probes, Infinium design, 450K, hg38 position and context per array, whether it passed QC per data source and was tested per mQTL run, and its overlap category. A lookup table for combining the data sources |
| `overview_05_array_overlap.png` | - | Sites on both arrays or one: in the manifests, after QC, and tested in the mQTL scan |
| `overview_06_array_context.png` | - | Island context of sites on both arrays and on each array only |
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

**Array comparison** (`array_comparison.R`):
- **Sites**: a site is a probe name without the EPICv2 suffix
  (`cg00002033_TC11` -> `cg00002033`), the same key QC uses to average EPICv2
  replicate probes (`R/mvalues.R`). Sites are matched between arrays by this
  name, as the rest of the pipeline does.
- **Manifests**: probes and sites come from the minfi annotation package of
  each array (`ARRAY_PROFILES` in `R/array_profiles.R`). The SNP (rs) control
  probes aren't in these annotations, so they aren't counted. EPICv2 also has
  `nv` probes, which target known somatic variants rather than CpGs.
- **Other ways to match**: the EPICv2 manifest gives each probe's EPICv1 name
  (`EPICv1_Loci`). For about 1,000 probes it differs from the probe's own name,
  so matching by name misses or mismatches these. A few hundred shared sites
  also have a different hg38 position on the two arrays (`CPG_ANNOTATION`; the
  EPICv1 positions are lifted over). Both are reported; exclude these sites
  before combining results across arrays.
- **Probe design**: a site measured with an Infinium I probe on one array and a
  type II probe on the other can have a different value range, which matters
  when comparing methylation levels, less so for associations within each
  array.
- **After QC and mQTL**: a site passes QC in a data source if it's a row of
  that source's `MVALS_QC`. The mQTL counts use each ancestry's
  `mqtl_cpgs_tested.tsv.gz` (CpGs with at least one cis SNP). Sites in both
  are the CpGs a cross-array or cross-ancestry comparison can use.
- **Context**: island context and gene region come from `CPG_ANNOTATION`; for
  sites on both arrays, from the EPICv2 annotation.
- **Two arrays**: the script expects the data sources to use exactly two
  arrays, one per data source, and stops otherwise.
