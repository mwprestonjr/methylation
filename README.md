# methylation
methylation data and scripts 

Analysis pipeline for Illumina EPIC methylation arrays (EPICv1 and EPICv2). Each data source is processed separately, with the same QC code; results from different sources can then be combined for downstream analyses.

## Data sources

| `SOURCE` | Data | Array | Sample sheet script |
|---|---|---|---|
| `ppmi_p140` | PPMI Project 140 | EPICv1 | `qc/01_build_sample_sheet_ppmi.R` |
| `psomagen`  | Psomagen deliveries (`DATASETS` in `config.R`) | EPICv2 | `qc/01_build_sample_sheet_psomagen.R` |

Input paths and the output directory for each source are set in `config.R`.

## Running the pipeline

The pipeline is run with `make` from the repository root, with the
`methylation` conda environment active. Each data source is run separately:

```bash
make qc            SOURCE=ppmi_p140                # sample sheet, QC/normalization, QC plots
make preprocessing SOURCE=ppmi_p140                # cell counts, sources of variation, ComBat
make clocks        SOURCE=ppmi_p140                # epigenetic clocks, age acceleration, group comparison, clinical traits
make mqtl          SOURCE=ppmi_p140 ANCESTRY=EUR   # sample map, genotypes, cis-mQTLs
make all           SOURCE=ppmi_p140                # qc, preprocessing, and clocks
make all           SOURCE=ppmi_p140 ANCESTRY=EUR   # qc, preprocessing, clocks, and mqtl
make                                               # usage
```

| Target | Steps, in order |
|---|---|
| `qc` | `qc/01_build_sample_sheet_<source>.R` → `qc/02_qc.R` → `qc/03_qc_plots.R` |
| `preprocessing` | `1_cell_counts.R` → `2_variation_sources.R raw` → `3_combat.R` → `2_variation_sources.R combat` → `4_cpg_annotation.R` |
| `clocks` | `1_estimate.R` → `2_compare.R` → `3_associations.R` → `4_plots.R` |
| `mqtl` | `1_mqtl_sample_map.R` → `2_mqtl_genotypes.sh` (plink2) → `3_mqtl.R` → `4_mqtl_plots.R` |

Each target needs the outputs of the one before it (`preprocessing` reads the
`qc` outputs; `clocks` and `mqtl` read the `preprocessing` outputs), so run
them in this order the first time; afterwards a target can be rerun on its own.

- **`SOURCE`** (required): one of the data sources above.
- **`ANCESTRY`** (`mqtl` only): the GP2 master key ancestry label to analyse,
  e.g. `EUR` for PPMI or `AFR` for Psomagen. The genotype file for that
  ancestry is found from `GENO_PFILE_PATH`/`GENO_PFILE_NAME` in `mqtl/config.R`, and
  results go to `results/mqtl/<ANCESTRY>/` and figures to
  `figures/mqtl/<ANCESTRY>/`. Several labels (`EUR,AJ`) only
  work if one genotype file covers all of them.
- **Logs**: each step writes `logs/<module>/<step>_<SOURCE>_<timestamp>.log`
  (mQTL steps start with the ancestry); all steps of one run share the
  timestamp. `make` itself prints one line per step.
- **Failures**: `make` stops at the first step that fails; see that step's log.
- **Steps never run in parallel**, even with `make -j`, because they are
  memory-heavy and each depends on the previous one.
- **Editing during a run**: the R steps are run through `source()`, which reads
  the whole script before running it, so editing a script mid-run doesn't
  affect the run. The mQTL step 2 shell script is read as it runs, so don't
  edit it while it's running.
- **Long runs**: run in the background so the run survives logging out:
  ```bash
  nohup make all SOURCE=ppmi_p140 > logs/make_ppmi_p140.log 2>&1 < /dev/null &
  ```
- **Latent PCs for mQTL**: `3_mqtl.R` adjusts for `MQTL_N_METH_PCS` latent
  methylation PCs (unmeasured variation), whose best number depends on the
  data. `Rscript mqtl/tune_meth_pcs.R <source> <ancestry>` (after `make mqtl`)
  reruns the cis scan for each value in
  `MQTL_N_METH_PCS_GRID` (`mqtl/config.R`), saves the counts of CpGs with a
  cis-mQTL (`tune_meth_pcs.csv`) and plots them (`mqtl_tune_meth_pcs.png`).
  Choose the value where the curve levels off, set `MQTL_N_METH_PCS`, and
  rerun `3_mqtl.R`. Tune each data source/ancestry separately.
- **mQTL figures** (`4_mqtl_plots.R`, in `figures/mqtl/<ANCESTRY>/`): genotype
  boxplots of the top hits at distinct loci, regional (LocusZoom-style) plots
  of the top loci, genome-wide lead p-values, lead SNP-CpG distance, effect
  size by allele frequency, mQTLs by CpG island context, and
  `mqtl_top_hits.csv`. Numbers of hits shown and the regional window are set in
  `mqtl/config.R` (`MQTL_PLOT_*`). Regional plots recompute every SNP in the
  window with the step 3 model, since step 3 saves only pairs with
  p < `MQTL_P_CIS_SAVE`
- **mQTL genotypes**: the genotype files must be readable from the VM (e.g. the
  `gp2_release12` bucket mounted at `~/gp2_release12` with gcsfuse), and
  `GENO_SOURCE` in `mqtl/config.R` must match them (`nba` for the NBA array files).

Every target reruns all of its steps, even if their outputs already exist; QC
of a few hundred samples takes about 2 hours.

### Running a single script

Scripts can also be run directly from the repository root, with the data
source (and, for some, a second argument) on the command line:

```bash
Rscript qc/02_qc.R ppmi_p140
Rscript preprocessing/2_variation_sources.R ppmi_p140 combat
Rscript mqtl/1_mqtl_sample_map.R psomagen AFR
```

In an interactive R session, set `DATA_SOURCE <- "ppmi_p140"` (and
`MQTL_ANCESTRY <- "AFR"` for mQTL scripts) with the working directory at the
repository root before sourcing a script.

## Layout

```
config.R          shared settings: data sources and their paths, output folders, files passed between modules
Makefile          runs the pipeline (see above)
R/                shared functions, sourced by the scripts
qc/               sample sheets, QC and normalization, QC plots (settings: qc/config.R)
preprocessing/    cell counts, sources of variation, ComBat batch correction (settings: preprocessing/config.R)
clocks/           epigenetic clocks, age acceleration, group comparison, clinical traits (settings: clocks/config.R)
mqtl/             cis-mQTL mapping against GP2 genotypes (settings: mqtl/config.R)
logs/             step logs from make
```

### Configuration

Settings are split between a shared config and one config per module. Every
script sources `config.R` first, then its own module's config:

- **`config.R`** (repository root): anything more than one module uses - the
  data sources (input paths, `DIR_OUTPUT`), figure settings, the files one
  module writes and another reads (e.g. `SAMPLE_SHEET_QC`, `COMBAT_MVALS`,
  `SAMPLE_SHEET_FINAL`, `CPG_ANNOTATION`), the R12 master key
  (`FNAME_METADATA`), `SAMPLE_SHEET_COLUMNS`, `CELL_TYPES`,
  and the hg19->hg38 liftover chain and GENCODE gene settings
- **`qc/config.R`**: QC thresholds, SeSAMe and density outlier settings
- **`preprocessing/config.R`**: sources-of-variation settings, the ComBat batch
  variable (per data source) and protected variables
- **`clocks/config.R`**: clock coverage thresholds, the clocks used per data
  source, the groups compared, the model covariates, and the clinical traits
  tested in cases
- **`mqtl/config.R`**: genotype files and source, genotype QC, cis window,
  model, multiple-testing and plot settings

Rule of thumb: a setting read by another module belongs in `config.R`;
otherwise it goes in its module's config.

Output folders, per data source (`DIR_OUTPUT` in `config.R`), one per module
(`DIR_RESULTS_<MODULE>` and `DIR_FIGURES_<MODULE>` in `config.R`):

```
<DIR_OUTPUT>/                e.g. /mnt/output/methylation/psomagen
├── results/
│   ├── qc/                  sample sheets, normalized mSet, M and beta values, QC metrics
│   ├── preprocessing/       cell proportions, ComBat M values, final sample sheet, CpG annotation, variation sources
│   ├── clocks/              clock estimates, accuracy, coverage, group comparison, clinical trait associations
│   └── mqtl/<ANCESTRY>/     sample map, genotypes/, cis-mQTL results, summary, top hits, tuning table
└── figures/
    ├── qc/
    ├── preprocessing/
    ├── clocks/
    └── mqtl/<ANCESTRY>/
```

### Script inventory

**Pipeline and shared configuration**

- `Makefile` - runs the module workflows in order and writes step logs.
- `config.R` - shared data-source paths, output paths, and settings.
- `requirements.R` - installs the R package dependencies.

**QC (`qc/`)**

- `01_build_sample_sheet_ppmi.R` - builds the PPMI sample sheet.
- `01_build_sample_sheet_psomagen.R` - builds the Psomagen sample sheet. Both sample-sheet scripts write the common `SAMPLE_SHEET_COLUMNS` from `config.R`, so downstream scripts work independently of data source.
- `02_qc.R` - performs sample and probe QC and normalization. It supports any array in `R/array_profiles.R`, with one array type processed per run.
- `03_qc_plots.R` - generates QC figures from the figure data saved by `02_qc.R`.
- `config.R` - QC thresholds and method settings.

**Preprocessing (`preprocessing/`)**

- `1_cell_counts.R` - estimates blood cell proportions.
- `2_variation_sources.R` - tests sources of variation before ComBat to choose `COMBAT_BATCH_VAR`, then reruns after ComBat to check the correction.
- `3_combat.R` - applies ComBat batch correction.
- `4_cpg_annotation.R` - creates `CPG_ANNOTATION`, an annotation for every array probe with hg38 position (liftover for EPICv1), GENCODE gene and gene type, promoter/gene-body/intergenic context (nearest gene for intergenic probes), distances to the gene and its TSS, and CpG island context.
- `config.R` - variation-analysis and ComBat settings.

**Epigenetic clocks (`clocks/`)**

- `1_estimate.R` - estimates clock ages, pace of aging, and age acceleration. Run the clock pipeline separately for each data source.
- `2_compare.R` - compares clock outcomes between diagnosis groups.
- `3_associations.R` - tests clock outcomes against clinical traits in cases.
- `4_plots.R` - generates clock, group-comparison, and trait-association figures.
- `config.R` - clock selection, model, and clinical-trait settings.

**mQTL (`mqtl/`)**

Run steps 1-4 separately for each data source and ancestry.

- `1_mqtl_sample_map.R` - creates the ancestry-specific methylation/genotype sample map.
- `2_mqtl_genotypes.sh` - prepares and QC-filters genotypes with PLINK 2.
- `3_mqtl.R` - performs the cis-mQTL scan.
- `4_mqtl_plots.R` - summarizes and plots mQTL results.
- `tune_meth_pcs.R` - reruns scans across methylation-PC counts to support PC selection.
- `config.R` - genotype, QC, scan, and plot settings.

**Shared R helpers (`R/`)**

- `array_profiles.R` - array-specific annotation packages, genome builds, cross-reactive probes, and array detection from IDAT files.
- `idat_qc.R` - QC helpers for IDAT file checks, SeSAMe QC statistics, sex checks, and low-bead-count probes.
- `density.R` - computes beta-value density curves and the density outlier score.
- `mvalues.R` - caps infinite M values from beta values of exactly 0 or 1.
- `cpg_annotation.R` - provides hg38 positions for array probes (with EPICv1 liftover) and GENCODE gene annotations; `preprocessing/4_cpg_annotation.R` saves `CPG_ANNOTATION` for every array probe, including gene and gene type, promoter/gene-body/intergenic context (nearest gene for intergenic probes), distances to the gene and its TSS, and CpG island context.
- `mqtl_setup.R` - handles the mQTL ancestry argument and configures ancestry-specific genotype inputs and output paths.
- `clock_models.R` - shared clock association models for diagnosis-group comparisons and clinical-trait associations.
- `clock_traits.R` - derives clinical traits at the blood draw for clock association analyses.
- `mqtl_functions.R` - functions shared by `3_mqtl.R` and `tune_meth_pcs.R` for data preparation, latent methylation PCs, and per-chromosome cis scans.

### Figure inventory

Figures are saved under `<DIR_OUTPUT>/figures/` for the data source being
processed (`DIR_OUTPUT` in `config.R`).

**QC (`figures/qc/`, generated by `qc/03_qc_plots.R`)**

- `qc_01_median_intensities.png` - methylated vs unmethylated median intensity per sample.
- `qc_02_detection_pvalues.png` - mean detection p-value per sample and the QC threshold.
- `qc_02b_sesame_detection.png` - fraction of detected CpG probes per sample; generated only when SeSAMe QC is enabled.
- `qc_03_sex_prediction.png` - predicted vs reported sex from X/Y intensities.
- `qc_04_density_before_normalization.png` - beta-value densities before normalization, by phenotype.
- `qc_04b_density_by_dataset.png` - beta-value densities before normalization, by dataset.
- `qc_05_density_after_normalization.png` - beta-value densities after normalization, by phenotype.

**Preprocessing (`figures/preprocessing/`)**

- `cell_proportions.png` - estimated blood-cell proportions across samples.
- `cell_proportions_by_phenotype.png` - estimated blood-cell proportions by phenotype.
- `variation_<raw|combat>_svd_heatmap.png` - associations between sample variables and SVD components, before (`raw`) or after (`combat`) ComBat.
- `variation_<raw|combat>_scree.png` - variance explained by SVD components, before or after ComBat.
- `variation_<raw|combat>_pca_by_<variable>.png` - PC1 vs PC2, colored by each available variable among `Chip`, `Batch`, `phenotype`, `sex`, and `age`; generated only for variables with usable variation in the data.

**Epigenetic clocks (`figures/clocks/`, generated by `clocks/4_plots.R`)**

- `clock_01_vs_chronological_age.png` - clock age vs chronological age, with accuracy and CpG coverage.
- `clock_02_acceleration_by_group.png` - age acceleration and DunedinPACE by phenotype.
- `clock_03_group_effects.png` - adjusted differences between diagnosis groups; generated when a group comparison is available.
- `clock_04_trait_effects.png` - adjusted clinical-trait associations in cases; generated when traits were tested.

**mQTL (`figures/mqtl/<ANCESTRY>/`)**

- `mqtl_01_lead_snp_distance.png` - distance between each significant lead SNP and its CpG; generated by `mqtl/3_mqtl.R`.
- `mqtl_02_top_genotype_boxplots.png` - methylation by genotype for top mQTLs at distinct loci.
- `mqtl_03_regional.png` - regional association plots around top CpGs, with points colored by LD to the lead SNP.
- `mqtl_04_genome.png` - lead association p-values across the genome.
- `mqtl_05_distance_zoom.png` - lead SNP-to-CpG distances within 50 kb.
- `mqtl_06_effect_vs_maf.png` - lead-SNP effect size by minor allele frequency.
- `mqtl_07_island_context.png` - cis-mQTL share by CpG island context.
- `mqtl_tune_meth_pcs.png` - significant CpG counts across tested methylation-PC counts; generated separately by `mqtl/tune_meth_pcs.R`.

## Epigenetic clocks

`make clocks` estimates epigenetic age and pace of aging per sample, compares
them between diagnosis groups, and tests their association with clinical
traits in cases. Each data source is run separately:
the arrays measure different CpGs, so clock ages aren't comparable across them.

| Step | What it does | Outputs (`results/clocks/`) |
|---|---|---|
| `1_estimate.R` | Collapses EPICv2 probe names to 450K/EPIC IDs (replicate probes averaged), records each clock's CpG coverage, estimates the methylclock clocks and DunedinPACE from normalized, unfiltered betas, and computes age acceleration | `clock_age_acceleration.csv` (per sample, with the model covariates), `clock_accuracy.csv`, `clock_cpg_coverage.csv` |
| `2_compare.R` | Compares age acceleration and DunedinPACE between `CLOCK_GROUPS` (first = reference), adjusted for `CLOCK_COVARIATES`; FDR across outcomes | `clock_acceleration_by_diagnosis.csv` |
| `3_associations.R` | In cases (`CLOCK_CASE_GROUP`), tests each trait in `CLOCK_TRAITS` against age acceleration and DunedinPACE, adjusted for `CLOCK_COVARIATES`; FDR across all tests | `clock_case_traits.csv`, `clock_trait_summary.csv`, `clock_trait_associations.csv` |
| `4_plots.R` | Clock vs chronological age (with r, median absolute error, CpG coverage), acceleration by phenotype, adjusted group differences, clinical trait effects | `figures/clocks/clock_*.png` |

- **Age acceleration** is the residual of clock age regressed on chronological
  age, in years. DunedinPACE is already a rate (years of aging per year; 1 =
  average), so it isn't residualized.
- **Comparison model**: `outcome ~ phenotype + CLOCK_COVARIATES` for every
  outcome: age, sex, cell proportions (without neutrophils, as the proportions
  sum to ~1), chip row, and the ComBat batch variable. Clock betas aren't
  ComBat-corrected, so batch is adjusted for in the model. Categorical
  covariates with one value among the compared samples are left out.
- **Skipped comparison**: if a group in `CLOCK_GROUPS` has fewer than
  `CLOCK_MIN_GROUP_N` samples (e.g. Psomagen has no controls), `2_compare.R`
  skips the comparison and `4_plots.R` skips the group effects figure.
- **Clinical traits** (`3_associations.R`): measured at the blood draw, from
  the R12 master key (disease duration, family history), the GP2 extended
  clinical data (`CLINICAL_DATA`; first of `CLOCK_CLINICAL_VISIT_MONTHS` with
  a value, e.g. PPMI MoCA comes from screening), and for PPMI the monogenic PD
  gene (`PPMI_GENETICS`). Each trait is its own model; continuous traits are
  per unit (the figure shows per SD), categorical ones vs `ref`. Traits with
  fewer than `CLOCK_MIN_TRAIT_N` cases are skipped and categorical levels with
  fewer than `CLOCK_MIN_LEVEL_N` are left out; `clock_trait_summary.csv` says
  which and why. Psomagen's clinical data only has Hoehn & Yahr.
- **Multiple testing**: Benjamini-Hochberg FDR, across outcomes in
  `2_compare.R` and across all outcome x trait tests in `3_associations.R`;
  results with FDR < `CLOCK_FDR` are called significant (filled points in the
  figures). Clocks are correlated, so agreement across clocks matters more
  than any single FDR value.
- **Chronological age**: PPMI uses age at the baseline visit (PPMI
  `Age_at_visit`, falling back to enrollment age); Psomagen uses R12
  `age_at_sample_collection`, the age at the GP2 DNA draw that was arrayed.
- **Missing CpGs**: methylclock leaves CpGs that aren't on the array out of the
  clock's weighted sum rather than imputing them, which shifts clock ages
  (values missing in some samples are imputed). A clock is estimated only if
  `CLOCK_COV_THRESH` of its CpGs are on the array (`PACE_COV_THRESH` for
  DunedinPACE); `clock_cpg_coverage.csv` lists the coverage per clock. The
  telomere length clock (TL) isn't estimated on EPICv2 for this reason.
- **Clocks used per array** (`ADULT_CLOCKS`): all clocks are estimated and
  saved, but age acceleration, comparisons and plots use only `ADULT_CLOCKS`.
  For EPICv2 this leaves out Hannum and skinHorvath: EPICv2 lacks CpGs they
  weight heavily (e.g. 5 Hannum CpGs on EPICv1 but not EPICv2 add ~15 years
  to a typical sample), giving clock ages ~31 and ~10 years below
  chronological age. A constant offset wouldn't affect age acceleration, but
  the variation between people carried by the missing CpGs is lost.

## Setup
```bash
conda env create -f environment.yml
conda activate methylation
Rscript requirements.R
```

## Troubleshooting

### preprocessFunnorm pthread_create error on Linux/GCP VMs

If you encounter this error during normalization:

Reinstall `preprocessCore` with threading disabled:
```r
BiocManager::install("preprocessCore", 
                     configure.args = "--disable-threading",
                     force = TRUE)
```

## Batch variable derivation note

`Batch` in the sample sheet is the Psomagen dataset (delivery) ID, or for PPMI
the top-level idat folder name in the Basename path, e.g.
`.../idat/20190604_plate1/SENTRIXID/... -> Batch = "20190604_plate1"`
