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
| 2 | `qc_duplicates.csv` | | Every QC-passed array of a person with more than one: IDs, chip/position, predicted sex, detection metrics, lowest rs-probe correlation between the person's arrays, and whether it was kept, removed as a replicate, or removed as a possible sample swap; empty if there are none | |
| 2 | `qc_sesame_stats.csv` | | SeSAMe per-sample stats (detection, intensity, dye bias); only when `USE_SESAME_QC = TRUE` | |
| 2 | `qc_summary.csv` | | Number of samples and probes removed at each QC step | |
| 2 | `qc_figure_data.rds` | `QC_FIGURE_DATA` | Sample metrics, density curves and thresholds behind the QC figures | `3_qc_plots.R` |
| 2 | `sample_sheet_qc_passed.csv` | `SAMPLE_SHEET_QC` | Sample sheet of the samples that passed QC | preprocessing (all steps) |
| 2 | `mVals.rds` | `MVALS_QC` | M values: QC-passed samples, Funnorm-normalized, probe-filtered, EPICv2 replicate probes averaged per site, infinite values capped | `preprocessing/2_variation_sources.R`, `3_combat.R` |
| 2 | `bVals.rds` | | Beta values, same samples and rows as `mVals.rds`; only if `SAVE_BVALS = TRUE` | not read by the pipeline |
| 2 | `replicate_probes.csv` | `REPLICATE_PROBES` | EPICv2 replicate probes averaged in `mVals.rds`: site, probe, and the row name kept; empty for other arrays or with `COLLAPSE_REPLICATE_PROBES = FALSE` | mQTL (`R/mqtl_functions.R`) |
| 2 | `bVals_noob.rds` | `BVALS_NOOB` | Beta values: QC-passed samples (before density outlier removal), **noob-normalized**, all probes; the IDOL cell reference needs noob, not Funnorm, values | `preprocessing/1_cell_counts.R` |
| 2 | `bVals_unfiltered.rds` | `BVALS_UNFILTERED` | Beta values: QC-passed samples, normalized, **no probe filtering** (filters can remove clock CpGs) | `clocks/1_estimate.R` |
| 2 | `mSetSq_qc_passed.rds` | `MSET_QC` | Normalized, filtered `GenomicRatioSet` (minfi object with annotation); only if `SAVE_MSET_QC = TRUE` | not read by the pipeline |

## Methods and decisions

**Samples**
- **PPMI**: samples come from the Project 140 link list (the only map from
  idat to participant and visit), restricted to one visit
  (`PRIMARY_TIMEPOINT`, baseline) and the enrollment cohorts in
  `COHORTS_OF_INTEREST` (PD and healthy control). Phenotype comes from the PPMI
  enrollment cohort, not from R12 `GP2_phenotype`, which labels
  genetic-cohort unaffected carriers "Control" and later converters "PD".
- **PPMI duplicates**: participants with more than two arrays at the visit are
  removed (likely batch controls or mislabeled visits). Pairs (technical
  replicates) are kept for `2_qc.R` to choose between (see Sample QC).
  `clinical_id` is the PATNO, also for participants not in R12.
- **Psomagen**: all samples in the deliveries listed in `DATASETS`, with
  phenotype, sex and race from the R12 master key (`GP2_MASTER_KEY`) and age
  from the internal-use master key (`GP2_INTERNAL_KEY`), matched on GP2ID.
  Samples appearing in more than one delivery are reported and kept for
  `2_qc.R` to choose between (see Sample QC).
- **Age**: PPMI uses age at the methylation visit (falling back to enrollment
  age); Psomagen uses R12 `age_at_sample_collection`, the age at the GP2 DNA
  draw that was arrayed.
- **Array type** is read from one idat per chip (number of bead types), not
  taken from the delivery files. `2_qc.R` processes one array type per run.

**Sample QC** (`2_qc.R`), in order:
1. **Detection**: with `USE_SESAME_QC = TRUE` (default), samples with fewer
   than `SESAME_MIN_FRAC_DETECTED` (90%) of cg probes detected by SeSAMe
   pOOBAH are removed; otherwise samples with mean minfi detection p >
   `DETECTION_P_THRESHOLD` (0.01).
2. **Sex**: predicted sex (minfi `getSex`, X/Y intensities) is compared with
   reported sex; discordant samples are removed (`REMOVE_SEX_DISCORDANT`).
   Samples with unknown reported sex are never called discordant.
3. **Duplicate arrays**: one array is kept per person, so later modules
   don't count a person twice (in mQTL mapping, two arrays with the same
   genotypes). A person is identified by GP2ID, or `clinical_id` if there is
   no GP2ID. Among a person's arrays that passed steps 1-2:
   - **Sample swap check**: the arrays are compared on the rs SNP probes
     (59-65 per array, minfi `getSnpBeta`) and on predicted sex. If any two
     correlate below `DUPLICATE_MIN_SNP_COR` (0.9; true replicates are
     ~0.95 or higher, different people much lower) or differ in predicted
     sex, the arrays can't all be the same person, and there's no way to tell
     which one is, so all are removed.
   - Otherwise the best array is kept, by `DUPLICATE_METRIC`: `"frac_dt_cg"`
     (default; highest fraction of cg probes detected by SeSAMe pOOBAH, needs
     `USE_SESAME_QC = TRUE`) or `"mean_detP"` (lowest mean minfi detection
     p-value). Ties go to the more recent chip (higher Sentrix ID).
   - Only arrays passing steps 1-2 are compared, so a re-run of a failed
     array is used. Density outliers (step 4) aren't considered.

   The arrays and decisions are in `qc_duplicates.csv`, and the counts in
   `qc_summary.csv`. Removed arrays aren't normalized.
4. **Beta density outliers**: samples with a third peak in the beta
   distribution between 0.15 and 0.75 (e.g. a mixed sample) are scored after
   normalization and flagged above `DENSITY_MID_PEAK_MAX`. They are **flagged,
   not removed** (`REMOVE_DENSITY_OUTLIERS = FALSE`); see
   `qc_sample_metrics.csv`.

**Normalization**: functional normalization (minfi `preprocessFunnorm`) of the
samples passing steps 1-3. Density outliers are part of the normalization even
when removed afterwards.

**Probe filters**, applied after normalization:
- low bead count: fewer than `MIN_BEADS` (3) beads in more than 5% of samples
- cross-reactive: the `maxprobes` EPICv1 list; for EPICv2, mapped through the
  EPICv1 probe IDs, so **probes new on EPICv2 aren't screened**
- failed detection: detection p > `DETECTION_P_THRESHOLD` (0.01) in more than
  `FAILED_SAMPLE_CUTOFF` (10%, the ChAMP default) of samples
- SNP-overlapping: minfi `dropLociWithSnps` with its defaults: any probe
  with a known SNP (allele frequency > 0 in the annotation) at the CpG or the
  single-base extension site. SNPs elsewhere in the probe aren't filtered
- sex chromosomes (X and Y)

**Three versions of the normalized data are saved**:
- probe-filtered M values (`mVals.rds`), Funnorm, used by preprocessing and mQTL
- unfiltered betas (`bVals_unfiltered.rds`), Funnorm, for the clocks, because
  the probe filters can remove CpGs that clocks use
- unfiltered noob betas (`bVals_noob.rds`) for cell deconvolution, which needs
  noob rather than Funnorm values. All probes are kept so a different cell
  reference (other CpGs) doesn't need QC rerun

**M values**: betas of exactly 0 or 1 give infinite M values; these are capped
(`R/mvalues.R`) so ComBat and other models can use them.

**EPICv2 replicate probes**: EPICv2 measures some sites with more than one
probe; the name's suffix differs (e.g. `cg00002033_TC11` and `_TC12`; in
Psomagen about 4,800 sites have 2-6 probes after the probe filters, almost
all identical copies, a few of another Infinium type or strand). With
`COLLAPSE_REPLICATE_PROBES = TRUE` (default), `mVals.rds` (and `bVals.rds`)
have one value per site: the betas of the site's probes that passed the probe
filters are averaged, then converted to M values. The averaged row keeps the
name of the site's first probe (sorted), so lookups by probe name
(annotation, `CPG_ANNOTATION`) still work; `replicate_probes.csv` lists which
probes went into each site. The unfiltered and noob betas and
`mSetSq_qc_passed.rds` keep every probe (cell deconvolution averages
replicates itself). No effect on EPICv1/450k, which have no replicates.

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
