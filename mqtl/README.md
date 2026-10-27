# mQTL (`mqtl/`)

Maps cis-mQTLs (SNP-methylation associations within +/- 1 Mb) against GP2 genotypes, one data source and ancestry at a time. Reads the `preprocessing` outputs.

Run from the repository root with the `methylation` conda environment active
(see [Running the pipeline](../README.md#running-the-pipeline)):

```bash
make mqtl SOURCE=<source> ANCESTRY=<label>
Rscript mqtl/3_mqtl.R psomagen AFR   # a single script
```

Settings: [`config.R`](config.R) in this folder; shared paths and settings in
the root [`config.R`](../config.R).

## Before running

The genotype files must be readable from the VM (e.g. the `gp2_release12` bucket
mounted at `~/gp2_release12` with gcsfuse), and `GENO_SOURCE` in `mqtl/config.R`
must match them (`nba` for the NBA array files).

## Scripts

Run steps 1-4 separately for each data source and ancestry.

- `1_mqtl_sample_map.R` - creates the ancestry-specific methylation/genotype sample map.
- `2_mqtl_genotypes.sh` - prepares and QC-filters genotypes with PLINK 2.
- `3_mqtl.R` - performs the cis-mQTL scan.
- `4_mqtl_plots.R` - summarizes and plots mQTL results.
- `tune_meth_pcs.R` - reruns scans across methylation-PC counts to support PC selection.
- `config.R` - genotype, QC, scan, and plot settings.

## Outputs

In `<DIR_OUTPUT>/results/mqtl/<ANCESTRY>/` (`MQTL_DIR`, set in
`R/mqtl_setup.R`), so each ancestry has its own set.

| Step | File | Contents | Read by |
|---|---|---|---|
| 1 | `mqtl_sample_map.csv` (`MQTL_SAMPLE_MAP`) | Methylation samples with genotypes for this ancestry: GP2ID, clinical ID, methylation sample ID, phenotype, ancestry | steps 3, 4, `tune_meth_pcs.R` |
| 1 | `mqtl_keep.txt` (`MQTL_KEEP`) | plink2 `--keep` file (FID, IID) of those samples | step 2 |
| 2 | `genotypes/geno_qc.{pgen,pvar,psam}` | Genotypes after variant and sample QC (autosomal biallelic SNPs; `MQTL_MAF`, `MQTL_GENO_MISSING`, `MQTL_SAMPLE_MISSING`, `MQTL_HWE_P`), variant IDs `chr:pos:ref:alt` (hg38) | step 2 |
| 2 | `genotypes/geno_qc.rmdup.list` | Duplicate variant IDs removed | |
| 2 | `genotypes/ld_prune.prune.{in,out}` | LD-pruned variant list used for the genotype PCs | step 2 |
| 2 | `genotypes/geno_pcs.{eigenvec,eigenval}` | `MQTL_N_GENO_PCS` genotype PCs (ancestry covariates) | step 3, `tune_meth_pcs.R` |
| 2 | `genotypes/geno_chr<N>.traw` | Genotypes per chromosome, SNPs x samples, as 0/1/2 copies of the `COUNTED` allele (hard calls) | steps 3, 4, `tune_meth_pcs.R` |
| 2 | `genotypes/*.log` | plink2 logs | |
| 3 | `cis_mqtl_all_saved.tsv.gz` | Every SNP-CpG pair with p < `MQTL_P_CIS_SAVE`: effect, t, p, positions, distance, FDR | step 4 |
| 3 | `cis_mqtl_lead_per_cpg.tsv` | Lead (strongest) SNP per CpG, same columns; one row per CpG with any saved pair | step 4 |
| 3 | `mqtl_cpgs_tested.tsv.gz` | One row per CpG with an hg38 position: SNPs tested (`n_snps`; 0 = no cis SNP, so not tested), SNPs excluded in the probe (`n_probe_snps`) and the excluded range (`probe_lo`, `probe_hi`) | step 4 |
| 3 | `mqtl_covariates.csv` | Covariates used in the model, per sample in model order (known covariates, genotype PCs, latent methylation PCs) | step 4 |
| 3 | `mqtl_summary.csv` | Samples, CpGs tested (at least one cis SNP), window, PCs, number of tests, CpGs with a cis-mQTL (FDR and strict thresholds), FDR p-value cut-off | |
| 4 | `mqtl_top_hits.csv` | Top mQTLs at distinct loci with gene, island context and MAF | |
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

**Genotype QC** (`2_mqtl_genotypes.sh`, PLINK 2), on the analysed samples:
- autosomal, biallelic, single-base (A/C/G/T) SNPs only
- minor allele frequency >= `MQTL_MAF` (5%), variant missingness <=
  `MQTL_GENO_MISSING` (5%), sample missingness <= `MQTL_SAMPLE_MISSING` (5%),
  Hardy-Weinberg p >= `MQTL_HWE_P` (1e-6)
- variants renamed `chr:pos:ref:alt` (hg38); duplicates (assayed more than once
  on the NBA array) are reduced to one
- genotypes exported as hard calls: 0/1/2 copies of the counted allele

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

**Significance**:
- A CpG has a cis-mQTL if its lead SNP has Benjamini-Hochberg FDR <
  `MQTL_FDR` (0.05), with FDR computed over all SNP-CpG tests. A stricter set
  uses p < `MQTL_P_STRICT` (1e-8), comparable with GoDMC.
- Only pairs with p < `MQTL_P_CIS_SAVE` (1e-3) are saved. FDR is still computed
  over all tests, so the saved q-values are exact; `3_mqtl.R` warns if the FDR
  threshold lies beyond the saved range (counts would then be incomplete).
- **Limitation**: FDR over SNP-CpG pairs treats SNPs near a CpG as independent,
  although they're correlated through LD. The standard alternative,
  permutations per CpG (as in FastQTL, QTLtools or tensorQTL), isn't
  implemented yet.

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
  excluded: indels, SNPs with MAF < `MQTL_MAF` and variants not on the array
  stay in the probe unexcluded, and a common SNP in LD with one of them can
  show the artifact. Before this exclusion, lead SNPs within
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
- **Hard-call genotypes**: imputation uncertainty isn't carried into the model
  or the LD estimates.
- **Psomagen samples are PD cases** only (plus a few other parkinsonisms), so
  its mQTLs are estimated in cases.

## Figures

Saved in `<DIR_OUTPUT>/figures/mqtl/<ANCESTRY>/` for the data source being processed
(`DIR_OUTPUT` in the root `config.R`).

- `mqtl_01_lead_snp_distance.png` - distance between each significant lead SNP and its CpG; generated by `mqtl/3_mqtl.R`.
- `mqtl_02_top_genotype_boxplots.png` - methylation by genotype for top mQTLs at distinct loci.
- `mqtl_03_regional.png` - regional association plots around top CpGs, with points colored by LD to the lead SNP.
- `mqtl_04_genome.png` - lead association p-values across the genome.
- `mqtl_05_distance_zoom.png` - lead SNP-to-CpG distances within 50 kb.
- `mqtl_06_effect_vs_maf.png` - lead-SNP effect size by minor allele frequency.
- `mqtl_07_island_context.png` - cis-mQTL share by CpG island context.
- `mqtl_tune_meth_pcs.png` - significant CpG counts across tested methylation-PC counts; generated separately by `mqtl/tune_meth_pcs.R`.

The number of hits shown and the regional window are set in `mqtl/config.R`
(`MQTL_PLOT_*`). Regional plots recompute every SNP in the window with the
step 3 model, since step 3 saves only pairs with p < `MQTL_P_CIS_SAVE`: the
same covariates, missing genotypes set to the mean, and SNPs in the probe left
out. `4_mqtl_plots.R` warns if a recomputed lead p-value doesn't match step 3.
The island context plot's denominator is the CpGs step 3 tested
(`mqtl_cpgs_tested.tsv.gz`).

## Choosing the number of latent methylation PCs

`3_mqtl.R` adjusts for `MQTL_N_METH_PCS` latent methylation PCs (unmeasured
variation), whose best number depends on the data. After `make mqtl`, run

```bash
Rscript mqtl/tune_meth_pcs.R <source> <ancestry>
```

It reruns the cis scan for each value in `MQTL_N_METH_PCS_GRID`
(`mqtl/config.R`), saves the number of CpGs with a cis-mQTL for each
(`tune_meth_pcs.csv`) and plots them (`mqtl_tune_meth_pcs.png`). Choose the
value where the curve levels off, set `MQTL_N_METH_PCS`, and rerun `3_mqtl.R`.
Tune each data source and ancestry separately.
