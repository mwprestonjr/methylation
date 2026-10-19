# methylation
methylation data and scripts 

Analysis pipeline for Illumina EPIC methylation arrays (EPICv1 and EPICv2). Each data source is processed separately, with the same QC code; results from different sources can then be combined for downstream analyses.

## Data sources

| `SOURCE` | Data | Array | Sample sheet script |
|---|---|---|---|
| `ppmi_p140` | PPMI Project 140 | EPICv1 | `qc/1_build_sample_sheet_ppmi.R` |
| `psomagen`  | Psomagen deliveries (`DATASETS` in `config.R`) | EPICv2 | `qc/1_build_sample_sheet_psomagen.R` |

Input paths and the output directory for each source are set in `config.R`.

### Input files and versions

External files and reference data the pipeline reads. Each path is set by the
variable under "Name in code", in the file under "Set in"; update this table
when a file is replaced.

- **Source**: who provides the file and where it comes from
- **Version**: the release or version of the data (for R packages, the
  package version, and the GitHub commit for packages installed from GitHub)
- "-": not applicable or not recorded

**Methylation data and sample metadata**

| File | Name in code | Set in | Source | Version |
|---|---|---|---|---|
| PPMI Project 140 idat files | `IDAT_DIR` | `config.R` | PPMI (LONI), Project 140 | - |
| PPMI Project 140 link list (SENTRIXID to PATNO) | `LINK_LIST` | `config.R` | PPMI (LONI), Project 140 documentation | - |
| PPMI participant status (cohort, enrollment age) | `PARTICIPANT_STATUS` | `config.R` | PPMI (LONI) study data | - |
| PPMI age at visit | `AGE_AT_VISIT` | `config.R` | PPMI (LONI) study data | - |
| PPMI genetic status (LRRK2, GBA, ... carriers) | `PPMI_GENETICS` | `clocks/config.R` | PPMI (LONI) study data, Indiana University genetic consensus | 2025-10-25 consensus |
| Psomagen deliveries (idat files and QC tables) | `DATASETS`, `DIR_DELIVERY` | `config.R` | Psomagen, GP2 deliveries (`gp2-psomagen` bucket) | Deliveries AB00000952, AB00000963 |

**GP2 data**

| File | Name in code | Set in | Source | Version |
|---|---|---|---|---|
| GP2 master key (IDs, sex, race, phenotype, ages, family history) | `FNAME_METADATA` | `config.R` | GP2, internal-use master key | Release 12 (final) |
| GP2 master key with genotype availability, QC and ancestry labels | `GP2_MASTER_KEY` | `mqtl/config.R` | GP2 | Release 12 |
| GP2 extended clinical data (MDS-UPDRS, MoCA, H&Y, smoking, ...) | `CLINICAL_DATA` | `config.R` | GP2 | Release 12 |
| GP2 genotypes, plink2 files per ancestry, related individuals removed | `GENO_PFILE_PATH`, `GENO_PFILE_NAME`, `GENO_SOURCE` | `mqtl/config.R` | GP2 release bucket (`gp2_release12`) | Release 12, NBA array (imputed), hg38 |

**Reference data and annotation**

| File | Name in code | Set in | Source | Version |
|---|---|---|---|---|
| EPICv1 probe annotation (hg19) | `anno_pkg` | `R/array_profiles.R` | Bioconductor, `IlluminaHumanMethylationEPICanno.ilm10b4.hg19` | ilm10b4 (package 0.6.0) |
| EPICv2 probe annotation (hg38) | `anno_pkg` | `R/array_profiles.R` | Bioconductor, `IlluminaHumanMethylationEPICv2anno.20a1.hg38` | 20a1 (package 1.0.0) |
| Cross-reactive probe list (EPICv1, mapped to EPICv2) | `xreactive` | `R/array_profiles.R` | GitHub, `markgene/maxprobes` | 0.0.2 (commit `c2120db`) |
| Blood cell reference (IDOL) | `CELL_TYPES` | `config.R` | Bioconductor, `FlowSorted.Blood.EPIC` | 2.10.0 |
| hg19 to hg38 liftover chain (EPICv1 positions) | `CHAIN_HG19_HG38` | `config.R` | UCSC, `hg19ToHg38.chain.gz` | - |
| Gene annotation (GRCh38) | `GENCODE_GTF` | `config.R` | GENCODE, `gencode.v50.annotation.gtf.gz` | v50 |
| Epigenetic clock coefficients | `ADULT_CLOCKS` | `clocks/config.R` | Bioconductor, `methylclock` and `methylclockData` | 1.12.0 and 1.14.0 |
| DunedinPACE model | - | `clocks/1_estimate.R` | GitHub, `danbelsky/DunedinPACE` | 0.99.0 (commit `4b56998`) |

## Running the pipeline

The pipeline is run with `make` from the repository root, with the
`methylation` conda environment active. Each data source is run separately:

```bash
make qc            SOURCE=ppmi_p140                # sample sheet, QC/normalization, QC plots
make preprocessing SOURCE=ppmi_p140                # cell counts, sources of variation, ComBat
make clocks        SOURCE=ppmi_p140                # epigenetic clocks, age acceleration, group comparison, clinical traits
make ewas          SOURCE=ppmi_p140                # EWAS of case-control status, age, and an exposure
make mqtl          SOURCE=ppmi_p140 ANCESTRY=EUR   # sample map, genotypes, cis-mQTLs
make all           SOURCE=ppmi_p140                # qc, preprocessing, clocks, and ewas
make all           SOURCE=ppmi_p140 ANCESTRY=EUR   # qc, preprocessing, clocks, ewas, and mqtl
make                                               # usage
```

| Target | Steps, in order |
|---|---|
| `qc` | `qc/1_build_sample_sheet_<source>.R` → `qc/2_qc.R` → `qc/3_qc_plots.R` |
| `preprocessing` | `1_cell_counts.R` → `2_variation_sources.R raw` → `3_combat.R` → `2_variation_sources.R combat` → `4_cpg_annotation.R` |
| `clocks` | `1_estimate.R` → `2_compare.R` → `3_associations.R` → `4_plots.R` |
| `ewas` | `1_case_control.R` → `2_age.R` → `3_exposure.R` → `4_plots.R` |
| `mqtl` | `1_mqtl_sample_map.R` → `2_mqtl_genotypes.sh` (plink2) → `3_mqtl.R` → `4_mqtl_plots.R` |

Each target needs the outputs of the one before it (`preprocessing` reads the
`qc` outputs; `clocks`, `ewas` and `mqtl` read the `preprocessing` outputs), so run
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

Every target reruns all of its steps, even if their outputs already exist; QC
of a few hundred samples takes about 2 hours.

### Running a single script

Scripts can also be run directly from the repository root, with the data
source (and, for some, a second argument) on the command line:

```bash
Rscript qc/2_qc.R ppmi_p140
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
ewas/             epigenome-wide association studies: case-control, age, exposure (settings: ewas/config.R)
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
  (`FNAME_METADATA`) and extended clinical data (`CLINICAL_DATA`), `SAMPLE_SHEET_COLUMNS`, `CELL_TYPES`,
  and the hg19->hg38 liftover chain and GENCODE gene settings
- **`qc/config.R`**: QC thresholds, SeSAMe and density outlier settings
- **`preprocessing/config.R`**: sources-of-variation settings, the ComBat batch
  variable (per data source) and protected variables
- **`clocks/config.R`**: clock coverage thresholds, the clocks used per data
  source, the groups compared, the model covariates, and the clinical traits
  tested in cases
- **`ewas/config.R`**: EWAS covariates, the groups compared, the age column,
  the smoking variables and exposure tested (per data source), FDR and
  top-hit settings
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
│   ├── ewas/                per-CpG results, top hits and summary per model, combined model summary
│   └── mqtl/<ANCESTRY>/     sample map, genotypes/, cis-mQTL results, summary, top hits, tuning table
└── figures/
    ├── qc/
    ├── preprocessing/
    ├── clocks/
    ├── ewas/
    └── mqtl/<ANCESTRY>/
```

### Modules

Each module has its own README with its scripts and figures:

| Module | Does | README |
|---|---|---|
| `qc/` | sample sheets, sample and probe QC, normalization, QC plots | [qc/README.md](qc/README.md) |
| `preprocessing/` | cell proportions, sources of variation, ComBat, CpG annotation | [preprocessing/README.md](preprocessing/README.md) |
| `clocks/` | epigenetic clocks, age acceleration, group comparison, clinical traits | [clocks/README.md](clocks/README.md) |
| `ewas/` | CpG-wise association with case-control status, age and an exposure | [ewas/README.md](ewas/README.md) |
| `mqtl/` | cis-mQTL mapping against GP2 genotypes | [mqtl/README.md](mqtl/README.md) |

### Shared files

**Pipeline and shared configuration**

- `Makefile` - runs the module workflows in order and writes step logs.
- `config.R` - shared data-source paths, output paths, and settings.
- `requirements.R` - installs the R package dependencies.

**Shared R helpers (`R/`)**

- `array_profiles.R` - array-specific annotation packages, genome builds, cross-reactive probes, and array detection from IDAT files.
- `idat_qc.R` - QC helpers for IDAT file checks, SeSAMe QC statistics, sex checks, and low-bead-count probes.
- `density.R` - computes beta-value density curves and the density outlier score.
- `mvalues.R` - caps infinite M values from beta values of exactly 0 or 1.
- `cpg_annotation.R` - provides hg38 positions for array probes (with EPICv1 liftover) and GENCODE gene annotations; `preprocessing/4_cpg_annotation.R` saves `CPG_ANNOTATION` for every array probe, including gene and gene type, promoter/gene-body/intergenic context (nearest gene for intergenic probes), distances to the gene and its TSS, and CpG island context.
- `mqtl_setup.R` - handles the mQTL ancestry argument and configures ancestry-specific genotype inputs and output paths.
- `clock_models.R` - shared clock association models for diagnosis-group comparisons and clinical-trait associations.
- `clock_traits.R` - derives clinical traits at the blood draw for clock association analyses.
- `ewas_functions.R` - loads the `preprocessing` outputs for the EWAS, builds the covariates, fits the limma CpG-wise models, and writes each model's results, top hits and summary.
- `ewas_plot_functions.R` - Manhattan, QQ and effect / p-value distribution plots for `ewas/4_plots.R`.
- `mqtl_functions.R` - functions shared by `3_mqtl.R` and `tune_meth_pcs.R` for data preparation, latent methylation PCs, and per-chromosome cis scans.

## Setup
```bash
conda env create -f environment.yml
conda activate methylation
Rscript requirements.R
```
