# methylation
methylation data and scripts 

QC pipeline for Illumina EPIC methylation arrays (EPICv1 and EPICv2). Each data
source is processed separately, with the same QC code; results from different
sources can then be combined for downstream analyses.

## Data sources

| `DATA_SOURCE` | Data | Array | Sample sheet script |
|---|---|---|---|
| `ppmi_p140` | PPMI Project 140 | EPICv1 | `scripts/01_build_sample_sheet_ppmi.R` |
| `psomagen`  | Psomagen deliveries (`DATASETS` in `config.R`) | EPICv2 | `scripts/01_build_sample_sheet_psomagen.R` |

Input paths and the output directory for each source are set in `config.R`.

## Pipeline Order

Run scripts from the repository root, in this order:

```bash
Rscript scripts/01_build_sample_sheet_ppmi.R      # or _psomagen.R
Rscript scripts/02_qc.R ppmi_p140                 # QC, normalization, probe filtering
Rscript scripts/03_qc_plots.R ppmi_p140           # QC figures
```

Each script reads input files produced by the previous script. In an
interactive session, set `DATA_SOURCE <- "ppmi_p140"` (working directory = repo
root) before sourcing a script.

## Layout

- `config.R` - shared settings: QC thresholds, figure settings, output file
  names, and the paths for each data source
- `R/array_profiles.R` - everything that differs by array type (annotation
  package, cross-reactive probes, array detection from idat files)
- `scripts/01_build_sample_sheet_*.R` - one per data source; each writes the
  same columns (`SAMPLE_SHEET_COLUMNS` in `config.R`) so later scripts don't
  depend on the source
- `scripts/02_qc.R` - works for any array in `R/array_profiles.R`, one array
  type per run
- `scripts/03_qc_plots.R` - figures from the data saved by Script 02

Outputs go to `<DIR_OUTPUT>/results` and `<DIR_OUTPUT>/figures`.

## Setup
```bash
conda env create -f environment.yml
conda activate methylation
Rscript scripts/requirements.R
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
