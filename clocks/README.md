# Epigenetic clocks (`clocks/`)

Estimates epigenetic age and pace of aging per sample, compares them between
diagnosis groups, and tests their association with clinical traits in cases.
Reads the `qc` and `preprocessing` outputs. Each data source is run
separately: the arrays measure different CpGs, so clock ages aren't comparable
across them.

Run from the repository root with the `methylation` conda environment active
(see [Running the pipeline](../README.md#running-the-pipeline)):

```bash
make clocks SOURCE=<source>
Rscript clocks/3_associations.R ppmi_p140   # a single script
```

Settings: [`config.R`](config.R) in this folder; shared paths and settings in
the root [`config.R`](../config.R).

## Steps

| Step | What it does | Outputs (`results/clocks/`) |
|---|---|---|
| `1_estimate.R` | Collapses EPICv2 probe names to 450K/EPIC IDs (replicate probes averaged), records each clock's CpG coverage, estimates the methylclock clocks and DunedinPACE from normalized, unfiltered betas, and computes age acceleration | `clock_age_acceleration.csv` (per sample, with the model covariates), `clock_accuracy.csv`, `clock_cpg_coverage.csv` |
| `2_compare.R` | Compares age acceleration and DunedinPACE between `CLOCK_GROUPS` (first = reference), adjusted for `CLOCK_COVARIATES`; FDR across outcomes. Repeated without cell proportions | `clock_acceleration_by_diagnosis.csv`, `clock_acceleration_by_diagnosis_no_cells.csv` |
| `3_associations.R` | In cases (`CLOCK_CASE_GROUP`), tests each trait in `CLOCK_TRAITS` against age acceleration and DunedinPACE, adjusted for `CLOCK_COVARIATES`; FDR across all tests. Repeated without cell proportions | `clock_case_traits.csv`, `clock_trait_summary.csv`, `clock_trait_associations.csv`, `clock_trait_associations_no_cells.csv` |
| `4_plots.R` | Clock vs chronological age (with r, median absolute error, CpG coverage), acceleration by phenotype, adjusted group differences, clinical trait effects | `figures/clocks/clock_*.png` |
| `5_epicv2_simulation.R` | Optional, EPICv1 source only (`make clocks-epicv2-sim SOURCE=ppmi_p140`, after `make clocks` for both sources): re-estimates the clocks from only the CpGs also on EPICv2 and compares accuracy with the full-array estimates, with whole-year age, and with `EPICV2_SIM_SOURCE`; agreement with the full-array estimates; group comparison refit | `clock_epicv2_sim_age_acceleration.csv`, `clock_epicv2_sim_cpg_coverage.csv`, `clock_epicv2_sim_accuracy.csv`, `clock_epicv2_sim_agreement.csv`, `clock_epicv2_sim_by_diagnosis.csv`; `figures/clocks/clock_05_epicv2_simulation.png` |

## Methods and decisions

- **Age acceleration** is the residual of clock age regressed on chronological
  age, in years. DunedinPACE is already a rate (years of aging per year; 1 =
  average), so it isn't residualized.
- **Comparison model**: `outcome ~ phenotype + CLOCK_COVARIATES` for every
  outcome: age, sex, cell proportions (without neutrophils, as the proportions
  sum to ~1), chip row, and the ComBat batch variable. Clock betas aren't
  ComBat-corrected, so batch is adjusted for in the model. Categorical
  covariates with one value among the compared samples are left out.
- **Sensitivity analysis without cell proportions**: `2_compare.R` and
  `3_associations.R` repeat every model with `CLOCK_COVARIATES_NO_CELLS`
  (the covariates without the cell proportions), saved to the `_no_cells`
  files, with FDR within each model. Part of what the extrinsic clocks
  (Hannum, Levine) measure is immune cell composition, which the main
  model adjusts away; an effect seen only without cell adjustment points to
  cell composition rather than intrinsic aging. The figures show the main
  model only.
- **Skipped comparison**: if a group in `CLOCK_GROUPS` has fewer than
  `CLOCK_MIN_GROUP_N` samples (e.g. Psomagen has no controls), `2_compare.R`
  skips the comparison and `4_plots.R` skips the group effects figure.
- **Clinical traits** (`3_associations.R`): measured at the blood draw, from
  the R12 master key (`GP2_MASTER_KEY`; disease duration, family history), the GP2 extended
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

## Figures

Saved in `<DIR_OUTPUT>/figures/clocks/` for the data source being processed
(`DIR_OUTPUT` in the root `config.R`).

- `clock_01_vs_chronological_age.png` - clock age vs chronological age, with accuracy and CpG coverage.
- `clock_02_acceleration_by_group.png` - age acceleration and DunedinPACE by phenotype.
- `clock_03_group_effects.png` - adjusted differences between diagnosis groups; generated when a group comparison is available.
- `clock_04_trait_effects.png` - adjusted clinical-trait associations in cases; generated when traits were tested.
- `clock_05_epicv2_simulation.png` - EPICv1 source only, from `5_epicv2_simulation.R`: age acceleration and DunedinPACE from EPICv2 CpGs only vs all CpGs.
