# mQTL (`mqtl/`)

Maps cis-mQTLs (SNP-methylation associations within +/- 1 Mb) against GP2 genotypes, one data source and ancestry at a time. Reads the `preprocessing` outputs.
The scan uses GP2's imputed genotypes by default (`MQTL_GENOTYPES` in
`mqtl/config.R`), or the NBA array genotypes, as in the first analysis (see
Genotypes below).

Run from the repository root with the `methylation` conda environment active
(see [Running the pipeline](../README.md#running-the-pipeline)):

```bash
make mqtl SOURCE=<source> ANCESTRY=<label>
Rscript mqtl/3_mqtl.R psomagen AFR   # a single script
```

Step 4 runs chromosomes in parallel (`MQTL_PERM_CORES`), each worker with
single-threaded BLAS (`RhpcBLASctl`), so the workers don't compete for the
CPUs. It needs about 8 GB for the main process plus about 4 GB per core (~70
GB with 15 cores) and took 30-40 minutes for each of PPMI and Psomagen; lower
`MQTL_PERM_CORES` if other jobs share the VM. Those figures are for array genotypes;
imputed genotypes have about 15 times more SNPs per cis window, so steps 3 and
4 take much longer and step 4 needs more memory per core. Start with fewer
cores and watch memory on the first chromosomes.

Settings: [`config.R`](config.R) in this folder; shared paths and settings in
the root [`config.R`](../config.R).

## Before running

The genotype files must be readable from the VM (e.g. the `gp2_release12` bucket
mounted at `~/gp2_release12` with gcsfuse; `GP2_RELEASE_DIR` in
`mqtl/config.R`), and `GENO_SOURCE` must match the array files (`nba` for the
NBA array). Imputed runs also need the ancestry's imputed VCFs
(`MQTL_IMPUTED_VCFS`); `make check SOURCE=<source> ANCESTRY=<label>` checks all
of them.

## Scripts

Run steps 1-5 separately for each data source and ancestry.

- `1_mqtl_sample_map.R` - creates the ancestry-specific methylation/genotype sample map.
- `2_mqtl_genotypes.sh` - prepares and QC-filters genotypes with PLINK 2; for imputed genotypes, merges the imputation batches with `R/merge_imputed_dosages.R`.
- `3_mqtl.R` - performs the cis-mQTL scan.
- `4_mqtl_permutations.R` - CpG-level significance by permutation (see Significance).
- `5_mqtl_plots.R` - summarizes and plots mQTL results.
- `tune_meth_pcs.R` - reruns scans across methylation-PC counts to support PC selection.
- `config.R` - genotype, QC, scan, permutation and plot settings.

## Outputs

In `<DIR_OUTPUT>/results/mqtl/<ANCESTRY>_imputed/` for imputed genotypes, or
`<DIR_OUTPUT>/results/mqtl/<ANCESTRY>/` for array genotypes (`MQTL_DIR`, set
in `R/mqtl_setup.R`), so each ancestry and genotype set has its own set.

| Step | File | Contents | Read by |
|---|---|---|---|
| 1 | `mqtl_sample_map.csv` (`MQTL_SAMPLE_MAP`) | Methylation samples with genotypes for this ancestry: GP2ID, clinical ID, methylation sample ID, phenotype, ancestry (imputed: also imputation batch and VCF sample ID) | steps 3, 4, 5, `tune_meth_pcs.R` |
| 1 | `mqtl_keep.txt` (`MQTL_KEEP`) | plink2 `--keep` file (FID, IID) of those samples in the array genotypes | step 2 |
| 1 | `mqtl_keep_imputed_part<k>.txt` (`MQTL_KEEP_IMPUTED`) | Imputed only: plink2 `--keep` file (VCF sample IDs) per imputation batch holding any of the samples | step 2 |
| 2 | `genotypes/geno_qc.{pgen,pvar,psam}` | Genotypes after variant and sample QC (autosomal biallelic SNPs; `MQTL_MAF`, `MQTL_GENO_MISSING`, `MQTL_SAMPLE_MISSING`, `MQTL_HWE_P`), variant IDs `chr:pos:ref:alt` (hg38) | step 2 |
| 2 | `genotypes/geno_qc.rmdup.list` | Duplicate variant IDs removed | |
| 2 | `genotypes/ld_prune.prune.{in,out}` | LD-pruned variant list used for the genotype PCs | step 2 |
| 2 | `genotypes/geno_pcs.{eigenvec,eigenval}` | `MQTL_N_GENO_PCS` genotype PCs (ancestry covariates) | step 3, `tune_meth_pcs.R` |
| 2 | `genotypes/geno_chr<N>.traw` | Genotypes per chromosome, SNPs x samples, as copies of the `COUNTED` allele: dosages (0-2) for imputed genotypes, 0/1/2 hard calls for array genotypes | steps 3, 4, 5, `tune_meth_pcs.R` |
| 2 | `genotypes/imputed_variant_counts.tsv` | Imputed only: per chromosome, samples, variants per imputation batch after the R2 and batch MAF filters, variants in all batches, and after missingness and MAF in our samples | |
| 2 | `genotypes/*.log` | plink2 logs | |
| 3 | `cis_mqtl_all_saved.tsv.gz` | Every SNP-CpG pair with p < `MQTL_P_CIS_SAVE`: effect, t, p, positions, distance, FDR | steps 4, 5 |
| 3 | `cis_mqtl_lead_per_cpg.tsv` | Lead (strongest) SNP per CpG, same columns; one row per CpG with any saved pair | steps 4, 5 |
| 3 | `mqtl_cpgs_tested.tsv.gz` | One row per CpG with an hg38 position: SNPs tested (`n_snps`; 0 = no cis SNP, so not tested), SNPs excluded in the probe (`n_probe_snps`) and the excluded range (`probe_lo`, `probe_hi`) | steps 4, 5 |
| 3 | `mqtl_covariates.csv` | Covariates used in the model, per sample in model order (known covariates, genotype PCs, latent methylation PCs) | steps 4, 5 |
| 3 | `mqtl_summary.csv` | Samples, CpGs tested (at least one cis SNP), window, PCs, number of tests, CpGs with a cis-mQTL (FDR and strict thresholds), FDR p-value cut-off | |
| 4 | `cis_mqtl_perm_per_cpg.tsv.gz` | One row per tested CpG: lead SNP (as step 3: effect, t, p, position, distance), SNPs tested, direct permutation p-value (`p_perm`), effective df (`true_df`) and Beta shapes, Beta-approximated p-value (`p_beta`), its Storey and BH q-values (`q_storey`, `q_bh`; `q` is the one used, `MQTL_PERM_FDR_METHOD`) and the CpG's nominal p-value threshold (`p_nominal_threshold`) | step 5 |
| 4 | `cis_mqtl_perm_signif_pairs.tsv.gz` | Saved step 3 pairs of significant CpGs with p below the CpG's nominal threshold: every significant SNP, not only the lead | |
| 4 | `mqtl_perm_summary.csv` | Permutations, seed, FDR method, Storey pi0, CpGs with a cis-mQTL (Storey, BH and step 3's pair-level FDR), p_beta threshold, significant pairs | |
| 5 | `mqtl_top_hits.csv` | Top mQTLs at distinct loci with gene, island context and MAF (and permutation q-value and nominal threshold) | |
| tune | `tune_meth_pcs.csv` | CpGs with a cis-mQTL for each number of latent methylation PCs in `MQTL_N_METH_PCS_GRID` | |

## Methods and decisions

**Samples** (`1_mqtl_sample_map.R`): methylation samples are linked to GP2
genotypes by GP2ID through the GP2 master key (`GP2_MASTER_KEY`), keeping
samples genotyped on `GENO_SOURCE` (NBA array) that passed GP2 genotype QC and
have the requested ancestry label. The genotype files have related
individuals removed. Each person has one methylation array: `qc/2_qc.R`
keeps the best of duplicate arrays (see the QC README), and
`1_mqtl_sample_map.R` stops if a GP2ID still appears twice, since both
arrays would get the same genotypes and count as independent samples.

**Genotypes** (`MQTL_GENOTYPES`): the scan uses GP2's imputed dosages by
default. The NBA array has about 1.6M directly typed variants, so most GWAS
lead variants and their LD neighbors aren't on it; colocalization with PD risk
loci and comparison with imputation-based resources (GoDMC) need the imputed
variants. The array genotypes are still read in every run: they define the
samples (the array files have relatives removed; the imputed VCFs don't) and
give the genotype PCs, for which directly typed, LD-pruned SNPs are the
standard input. The first analysis used the array genotypes; its results are
kept in `<ANCESTRY>/` for comparison with `<ANCESTRY>_imputed/`.

The imputed data are Minimac4 output (hg38; dosage in the `DS` field,
imputation R2 in INFO). Where they are in GP2 R12 (`MQTL_IMPUTED_VCFS`):
- most ancestries, including AFR: `imputed_vcf_filtered_r2/<ANCESTRY>/`,
  already filtered by GP2 to R2 >= 0.3 (`bcftools filter -i R2>=0.3`, per
  the VCF header)
- EUR: not in that folder. EUR was imputed in four sample batches of about
  22,300 people each, whose raw (unfiltered) output is in
  `imputation_out_dir/EUR_part1` to `EUR_part4`; the `EUR_merged` folder there
  is empty. PPMI's EUR samples are spread over all four batches.

**Imputed variant QC** (`2_mqtl_genotypes.sh` and
`R/merge_imputed_dosages.R`), per chromosome and imputation batch, on the
analyzed samples:
- imputation R2 >= `MQTL_MIN_R2` (0.3), GP2's threshold for its filtered
  files, so EUR (filtered here) and the other ancestries (filtered by GP2) are
  treated alike. R2 comes from INFO, estimated by Minimac4 over the whole
  batch, far more precisely than it could be in our samples. A variant must
  pass in every batch our samples are in: R2 is estimated per batch, and a
  variant missing from one batch would be missing for all its samples
- MAF in the imputation batch (INFO) >= `MQTL_IMPUTED_MIN_BATCH_MAF` (1%),
  only to keep the files small
- the same SNP filters and IDs as the array: biallelic A/C/G/T SNPs, IDs
  `chr:pos:ref:alt`, duplicates reduced to one
- after joining the batches: missingness <= `MQTL_GENO_MISSING` and MAF >=
  `MQTL_MAF` (5%) in our samples
- no Hardy-Weinberg or sample missingness filter: dosages have no missing
  calls, and Hardy-Weinberg checks genotype calling, which GP2 did on the
  array data before imputation

The batches are merged in R because plink2 can't merge files with different
samples (non-concatenating `--pmerge` is not implemented). Reading each
batch's VCF is the slow part of step 2.

**Array genotype QC** (`2_mqtl_genotypes.sh`, PLINK 2), on the analyzed
samples; for the genotype PCs in every run, and the scan genotypes in array
runs:
- autosomal, biallelic, single-base (A/C/G/T) SNPs only
- minor allele frequency >= `MQTL_MAF` (5%), variant missingness <=
  `MQTL_GENO_MISSING` (5%), sample missingness <= `MQTL_SAMPLE_MISSING` (5%),
  Hardy-Weinberg p >= `MQTL_HWE_P` (1e-6)
- variants renamed `chr:pos:ref:alt` (hg38); duplicates (assayed more than once
  on the NBA array) are reduced to one
- genotypes exported as hard calls: 0/1/2 copies of the counted allele (array
  runs)

**Methylation**: ComBat-corrected M values, inverse-normal transformed per CpG
(each CpG's values ranked and mapped to a standard normal). This limits the
influence of outliers and puts every CpG on the same scale, so effects are in
standard deviations of methylation per allele. CpG positions are hg38 (EPICv1
lifted over, as in the CpG annotation).

**Model** (`3_mqtl.R`, MatrixEQTL, linear): methylation ~ SNP dosage +
covariates, for every SNP within `MQTL_CIS_WINDOW` (1 Mb) of each CpG, one
chromosome at a time. Covariates:
- age, sex, cell proportions (without neutrophils), phenotype (left out if all
  samples share one), chip row
- `MQTL_N_GENO_PCS` (5) genotype PCs from LD-pruned SNPs (`MQTL_PRUNE_WINDOW`
  500 kb, r2 < `MQTL_PRUNE_R2` 0.1), for ancestry structure
- `MQTL_N_METH_PCS` (10) latent methylation PCs: PCs of the methylation
  residuals after the other covariates (50,000 most variable CpGs), capturing
  unmeasured technical and biological variation; tuned with `tune_meth_pcs.R`
  (below)
- batch isn't a covariate, since the M values are already ComBat-corrected

SNP-CpG pairs with the SNP in the probe are excluded: a SNP there can change
probe binding or the extension base and mimic an mQTL. The excluded positions
(`probe_snp_regions()` in `R/mqtl_functions.R`) are the CpG itself and the 50 bases on the probe's
side of it, which covers the Infinium I and II probe bodies and their
extension bases. The probe's side comes from the probe strand in the sesame
hg38 manifest (`sesameData`, Zhou lab InfiniumAnnotation); CpGs missing from
it (a handful), and EPICv2 sites averaged from replicate probes on both
strands (about 30), are excluded on both sides. Excluded pairs aren't counted as
tests for the FDR, and the counts are in the log and `mqtl_summary.csv`.

The probe side was checked against the unfiltered results: significant pairs
were about twice as frequent 5-45 bp into the probe side as at the same
distance on the other side, falling to background by ~50 bp. The first ~4 bp
(the probe's 3' end) were depleted instead, presumably because probes are
designed to avoid common SNPs there.

Samples without age or sex are left out of the model. Missing genotype calls
are set to the SNP's mean in these samples (MatrixEQTL's handling).

**Significance** (`4_mqtl_permutations.R`): a CpG has a cis-mQTL if its
permutation q-value is < `MQTL_FDR` (0.05). This is the approach of FastQTL
and tensorQTL, used by GTEx:
- **One test per CpG.** Testing every SNP in a CpG's window and keeping the
  best inflates its p-value, more so with more SNPs, and the SNPs aren't
  independent tests because of LD. So each CpG's lead SNP is compared with
  the lead SNP after shuffling the samples, `MQTL_N_PERM` (1,000) times, which
  accounts for both. The shuffles (seed `MQTL_PERM_SEED`) are the same for
  every CpG. Methylation and genotypes are first residualized on the step 3
  covariates (read from `mqtl_covariates.csv`), so each shuffle tests the same
  model as step 3, with the same SNPs (probe SNPs excluded) and missing
  genotypes set to the mean. A shuffled residual is no longer orthogonal to
  the covariates, so it is residualized again, as tensorQTL does; without
  this the permuted correlations are too small (by about covariates/samples,
  ~13% in r² here), the null too weak and the p-values anti-conservative (in
  a pure-null simulation, 10% of CpGs at p_beta < 0.05). The step warns if
  the median `true_df` is well above the model's df, the visible sign of
  this. The step checks that its lead p-values and SNP
  counts match step 3's.
- **Beta approximation.** 1,000 shuffles can't give a p-value below 1/1001, so,
  as FastQTL, a Beta distribution is fitted to each CpG's 1,000 shuffled lead
  p-values and evaluated at the observed one (`p_beta`). As FastQTL, the
  degrees of freedom of those p-values are first adjusted (`true_df`) so the
  Beta's first shape is 1, a small correction to the Beta form.
  The direct permutation p-value (`p_perm`) is saved too, and `p_beta` should
  agree with it where 1,000 shuffles can resolve it
  (`mqtl_perm_02_beta_vs_empirical.png`).
- **FDR across CpGs.** Storey q-values of `p_beta` (`MQTL_PERM_FDR_METHOD =
  "storey"`, as GTEx): Benjamini-Hochberg multiplied by pi0, the estimated
  share of CpGs without a cis-mQTL, from the flat part of the `p_beta`
  histogram (`mqtl_perm_01_pvalues.png`). BH assumes pi0 = 1, which is
  conservative when many CpGs have mQTLs. pi0 is estimated as in
  `qvalue::pi0est` with its defaults (`storey_pi0()` in
  `R/mqtl_functions.R`), so q-values match the qvalue package without
  depending on it. Both are saved (`q_storey`, `q_bh`); set
  `MQTL_PERM_FDR_METHOD = "bh"` for the conservative call.
- **Significant SNPs per CpG.** The `p_beta` cut-off between significant and
  non-significant CpGs is mapped back through each CpG's Beta distribution to a
  nominal p-value threshold (`p_nominal_threshold`), as GTEx's
  `pval_nominal_threshold`. Every SNP below it is significant for that CpG
  (`cis_mqtl_perm_signif_pairs.tsv.gz`). Only step 3's saved pairs (p <
  `MQTL_P_CIS_SAVE`) can be listed. A CpG with few SNPs in its window has a
  looser threshold (with one SNP, about the `p_beta` cut-off itself), so a
  small number of significant CpGs (e.g. 32 in Psomagen AFR, 130 in PPMI EUR)
  have thresholds above `MQTL_P_CIS_SAVE` and incomplete pair lists; their
  per-CpG results are complete. The step reports the count, and warns only if
  it is over 1% of significant CpGs.
- **Step 3's pair-level FDR** (Benjamini-Hochberg over all SNP-CpG tests) is
  still computed, used by `tune_meth_pcs.R` (much faster, and fine for
  comparing PC counts) and as a fallback by `5_mqtl_plots.R` if step 4 hasn't
  been run. A stricter fixed set uses p < `MQTL_P_STRICT` (1e-8), comparable
  with GoDMC. Only pairs with p < `MQTL_P_CIS_SAVE` (1e-3) are saved; the
  pair-level FDR is still over all tests, so its saved q-values are exact, and
  `3_mqtl.R` warns if its threshold lies beyond the saved range.
- **Possible improvements**:
  - *Adaptive permutations* (FastQTL `--permute 1000 10000`): stop early for
    CpGs that are clearly null, and permute more where the direct p-value is
    near the significance cut-off. With the Beta approximation, a fixed 1,000
    is standard and enough, so this would mainly save time.
  - *Conditional (independent) signals* (tensorQTL `cis_independent`): a CpG
    can have several independent mQTLs; only the lead signal is tested here.
    Forward-backward stepwise regression on the lead SNPs, with the
    permutation threshold, would find the others.

**Known caveats**
- **CpGs with a SNP at the CpG or extension base aren't tested**: QC removes
  these probes (minfi `dropLociWithSnps`), since the array would read the
  genotype rather than methylation. So allele-specific methylation caused by
  alleles that destroy a CpG can't be studied with these data. The probes
  are chosen from a fixed list of common SNPs (dbSNP, included in the array
  annotation: build 137 for EPICv1, 151 for EPICv2), not from the study's
  genotypes, so the same probes are removed for every ancestry.
- **Probe-SNP artifacts**: SNP-CpG pairs with the SNP inside the probe are
  excluded (see Model), but a SNP outside the probe in LD with one inside it
  can still pick up the artifact. Only variants in the QC'd genotypes can be
  excluded: indels, SNPs with MAF < `MQTL_MAF` and, in array runs, variants
  not on the array stay in the probe unexcluded, and a common SNP in LD with
  one of them can show the artifact. Imputed genotypes cover far more of the
  SNPs in probes. Before this exclusion, lead SNPs within
  50 bp of their CpG were 3.5% of significant CpGs in Psomagen AFR and 1.7%
  in PPMI EUR.
- **Autosomes only**: QC removes X and Y probes, and step 2 keeps autosomal
  SNPs.
- **EPICv2 replicate probes are averaged per site** in QC (see the QC
  README), so each site is tested and counted once, under the name of its
  first probe (e.g. `cg00002033_TC11` for the average of `_TC11` and
  `_TC12`). Before averaging, the copies gave consistent mQTLs: of the 317
  replicated sites with two or more probes at p < 1e-8 in Psomagen AFR, 88%
  had the same lead SNP and 96% the same effect direction.
- **Cross-reactive probes new on EPICv2 aren't screened** (see the QC
  README): a probe that also binds elsewhere in the genome measures a mix of
  sites, so an mQTL can belong to a site other than the annotated CpG.
- **Lead SNP is not necessarily causal**: it tags a block of correlated SNPs
  (see the regional plots); identifying the causal variant needs
  fine-mapping.
- **Imputation**: imputed dosages carry imputation uncertainty into the
  model (a poorly imputed variant has dosages shrunk toward its mean, which
  weakens its association rather than creating one), down to R2 0.3. The
  imputation batches of EUR were imputed separately; if batch membership
  differed by phenotype or plate, batch-specific dosage differences could act
  like a batch effect, so check that batch (`imputed_part` in the sample map)
  isn't associated with phenotype or the methylation PCs.
- **Psomagen samples are PD cases** only (plus a few other parkinsonisms), so
  its mQTLs are estimated in cases.

## Figures

Saved in `<DIR_OUTPUT>/figures/mqtl/<ANCESTRY>[_imputed]/` for the data source being processed
(`DIR_OUTPUT` in the root `config.R`).

- `mqtl_01_lead_snp_distance.png` - distance between each significant lead SNP (step 3's pair-level FDR) and its CpG; generated by `mqtl/3_mqtl.R`.
- `mqtl_perm_01_pvalues.png` - histogram of the permutation p-values (`p_beta`) with Storey's pi0; generated by `mqtl/4_mqtl_permutations.R`.
- `mqtl_perm_02_beta_vs_empirical.png` - Beta-approximated against direct permutation p-values; generated by `mqtl/4_mqtl_permutations.R`.
- `mqtl_02_top_genotype_boxplots.png` - methylation by genotype for top mQTLs at distinct loci.
- `mqtl_03_regional.png` - regional association plots around top CpGs, with points colored by LD to the lead SNP and the CpG's nominal threshold.
- `mqtl_04_genome.png` - lead association p-values across the genome.
- `mqtl_05_distance_zoom.png` - lead SNP-to-CpG distances within 50 kb.
- `mqtl_06_effect_vs_maf.png` - lead-SNP effect size by minor allele frequency.
- `mqtl_07_island_context.png` - cis-mQTL share by CpG island context.
- `mqtl_tune_meth_pcs.png` - significant CpG counts across tested methylation-PC counts; generated separately by `mqtl/tune_meth_pcs.R`.

The number of hits shown and the regional window are set in `mqtl/config.R`
(`MQTL_PLOT_*`). Regional plots recompute every SNP in the window with the
step 3 model, since step 3 saves only pairs with p < `MQTL_P_CIS_SAVE`: the
same covariates, missing genotypes set to the mean, and SNPs in the probe left
out. `5_mqtl_plots.R` warns if a recomputed lead p-value doesn't match step 3.
The island context plot's denominator is the CpGs step 3 tested
(`mqtl_cpgs_tested.tsv.gz`). Figures 02-07 use the step 4 permutation
significance, or step 3's pair-level FDR (with a warning) if step 4 hasn't
been run.

## Choosing the number of latent methylation PCs

`3_mqtl.R` adjusts for `MQTL_N_METH_PCS` latent methylation PCs (unmeasured
variation), whose best number depends on the data. After `make mqtl`, run

```bash
Rscript mqtl/tune_meth_pcs.R <source> <ancestry>
```

It reruns the cis scan for each value in `MQTL_N_METH_PCS_GRID`
(`mqtl/config.R`), saves the number of CpGs with a cis-mQTL for each
(`tune_meth_pcs.csv`) and plots them (`mqtl_tune_meth_pcs.png`). Choose the
value where the curve levels off, set `MQTL_N_METH_PCS`, and rerun steps 3-5
(`3_mqtl.R`, `4_mqtl_permutations.R`, `5_mqtl_plots.R`). The counts use step
3's pair-level FDR, which is much faster than permutations and fine for
comparing PC counts.
Tune each data source and ancestry separately.
