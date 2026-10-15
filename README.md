# methylation
methylation data and scripts 

Analysis pipeline for Illumina EPIC methylation arrays (EPICv1 and EPICv2). Each data source is processed separately, with the same QC code; results from different sources can then be combined for downstream analyses.

## Data sources

| `SOURCE` | Data | Array | Sample sheet script |
|---|---|---|---|
| `ppmi_p140` | PPMI Project 140 | EPICv1 | `qc/1_build_sample_sheet_ppmi.R` |
| `psomagen`  | Psomagen deliveries (`DATASETS` in `config.R`) | EPICv2 | `qc/1_build_sample_sheet_psomagen.R` |

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
| `qc` | `qc/1_build_sample_sheet_<source>.R` → `qc/2_qc.R` → `qc/3_qc_plots.R` |
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

### Modules

Each module has its own README with its scripts and figures:

| Module | Does | README |
|---|---|---|
| `qc/` | sample sheets, sample and probe QC, normalization, QC plots | [qc/README.md](qc/README.md) |
| `preprocessing/` | cell proportions, sources of variation, ComBat, CpG annotation | [preprocessing/README.md](preprocessing/README.md) |
| `clocks/` | epigenetic clocks, age acceleration, group comparison, clinical traits | [clocks/README.md](clocks/README.md) |
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
- `mqtl_functions.R` - functions shared by `3_mqtl.R` and `tune_meth_pcs.R` for data preparation, latent methylation PCs, and per-chromosome cis scans.

## Setup
```bash
conda env create -f environment.yml
conda activate methylation
Rscript requirements.R
```
