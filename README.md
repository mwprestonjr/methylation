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
make mqtl          SOURCE=ppmi_p140 ANCESTRY=EUR   # sample map, genotypes, cis-mQTLs
make mqtl-tune     SOURCE=ppmi_p140 ANCESTRY=EUR   # choose the number of latent PCs (after mqtl)
make all           SOURCE=ppmi_p140                # qc and preprocessing
make all           SOURCE=ppmi_p140 ANCESTRY=EUR   # qc, preprocessing, and mqtl
make                                               # usage
```

| Target | Steps, in order |
|---|---|
| `qc` | `qc/01_build_sample_sheet_<source>.R` → `qc/02_qc.R` → `qc/03_qc_plots.R` |
| `preprocessing` | `1_cell_counts.R` → `2_variation_sources.R raw` → `3_combat.R` → `2_variation_sources.R combat` |
| `mqtl` | `1_mqtl_sample_map.R` → `2_mqtl_genotypes.sh` (plink2) → `3_mqtl.R` |
| `mqtl-tune` | `tune_meth_pcs.R` (needs the `mqtl` genotype outputs) |

Each target needs the outputs of the one before it (`preprocessing` reads the
`qc` outputs, `mqtl` reads the `preprocessing` outputs), so run them in this
order the first time; afterwards a target can be rerun on its own.

- **`SOURCE`** (required): one of the data sources above.
- **`ANCESTRY`** (`mqtl` only): the GP2 master key ancestry label to analyse,
  e.g. `EUR` for PPMI or `AFR` for Psomagen. The genotype file for that
  ancestry is found from `GENO_PFILE_PATH`/`GENO_PFILE_NAME` in `mqtl/config.R`, and
  results go to `<results>/mqtl/<ANCESTRY>/`. Several labels (`EUR,AJ`) only
  work if one genotype file covers all of them.
- **Logs**: each step writes `logs/<step>_<SOURCE>_<timestamp>.log`; all steps
  of one run share the timestamp. `make` itself prints one line per step.
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
  data. `make mqtl-tune` reruns the cis scan for each value in
  `MQTL_N_METH_PCS_GRID` (`mqtl/config.R`), saves the counts of CpGs with a
  cis-mQTL (`tune_meth_pcs.csv`) and plots them (`mqtl_tune_meth_pcs.png`).
  Choose the value where the curve levels off, set `MQTL_N_METH_PCS`, and
  rerun `3_mqtl.R`. Tune each data source/ancestry separately.
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
mqtl/             cis-mQTL mapping against GP2 genotypes (settings: mqtl/config.R)
logs/             step logs from make
```

- `R/array_profiles.R` - everything that differs by array type (annotation
  package, genome build, cross-reactive probes, array detection from idat files)
- `R/idat_qc.R` - QC helpers (idat file check, SeSAMe QC stats, sex check, low
  bead count probes)
- `R/density.R` - beta density curves and the density outlier score
- `R/mvalues.R` - caps infinite M values (betas of exactly 0 or 1)
- `R/mqtl_setup.R` - mQTL ancestry argument, genotype file and output paths
- `R/mqtl_functions.R` - mQTL steps shared by `3_mqtl.R` and `tune_meth_pcs.R`
  (data preparation, latent PCs, per-chromosome cis scan)
- `qc/01_build_sample_sheet_*.R` - one per data source; each writes the same
  columns (`SAMPLE_SHEET_COLUMNS` in `config.R`) so later scripts don't depend
  on the source
- `qc/02_qc.R` - works for any array in `R/array_profiles.R`, one array type
  per run; `qc/03_qc_plots.R` draws its figures from the data it saves
- `preprocessing/2_variation_sources.R` - run before ComBat to choose the batch
  variable (`COMBAT_BATCH_VAR`) and after it to check the correction
- `mqtl/` - steps 1-3 of the mQTL analysis, run per data source and ancestry

### Configuration

Settings are split between a shared config and one config per module. Every
script sources `config.R` first, then its own module's config:

- **`config.R`** (repository root): anything more than one module uses - the
  data sources (input paths, `DIR_OUTPUT`), figure settings, the files one
  module writes and another reads (e.g. `SAMPLE_SHEET_QC`, `COMBAT_MVALS`,
  `SAMPLE_SHEET_FINAL`), `SAMPLE_SHEET_COLUMNS` and `CELL_TYPES`
- **`qc/config.R`**: QC thresholds, SeSAMe and density outlier settings, the
  R12 master key used for the sample sheets
- **`preprocessing/config.R`**: sources-of-variation settings, the ComBat batch
  variable (per data source) and protected variables
- **`mqtl/config.R`**: genotype files and source, liftover chain, cis window
  and model settings

Rule of thumb: a setting read by another module belongs in `config.R`;
otherwise it goes in its module's config.

Outputs go to `<DIR_OUTPUT>/results` and `<DIR_OUTPUT>/figures` for each data
source (`DIR_OUTPUT` in `config.R`).

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
