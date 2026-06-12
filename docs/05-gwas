---
layout: default
title: GWAS
---

[← Back to Home]({{ '/' | relative_url }})

# GWAS

To build polygenic scores we first need effect size estimates, which come from genome-wide association studies (GWAS) run in each UK Biobank GWAS panel. We run GWAS with [regenie](https://rgcgithub.github.io/regenie/), under a factorial set of analysis choices so that we can later compare how each choice affects stratification bias:

- **Correction** (`covars`): with principal components (`pcs`) or without (`nopcs`).
- **Model** (`gtype`): linear regression (`LR`) or a linear mixed model (`LMM`).

Crossing these gives four GWAS configurations per phenotype per panel, for each of the 17 phenotypes. The resulting effect sizes feed directly into the [polygenic score association tests](06-pgs-association-tests).

---

## Overview

The GWAS step has four parts:

1. **Format the phenotype and covariate files** for regenie.
2. **Fit the whole-genome model (regenie step 1)** to produce the leave-one-chromosome-out predictions.
3. **Run the association tests (regenie step 2)** under each model (`LR`/`LMM`) and correction (`pcs`/`nopcs`).
4. **Format the summary statistics** into a consistent SNP/effect table.

---

## Step 1: Format phenotypes and covariates

Phenotypes are assembled into a single regenie-format table (FID, IID, then one column per trait) by [`format_regenie_phenotype.R`](https://github.com/jgblanc/quantifying-susceptibility-PGS.github.io/blob/master/scripts/run_gwas/format_regenie_phenotype.R):

```
rule format_phenotypes:
    input:
        "output/phenotypes/participant_metadata_all.csv"
    output:
        "output/phenotypes/ukb_phenotypes_QT.txt"
    shell:
        """
        Rscript code/run_GWAS/format_regenie_phenotype.R {input} {output}
        """
```

Covariates come in two flavors. The **no-PC** version keeps only age, sex, and genotyping batch ([`format_regenie_covars_no_PCs.R`](https://github.com/jgblanc/quantifying-susceptibility-PGS.github.io/blob/master/scripts/run_gwas/format_regenie_covars_no_PCs.R)):

```r
colnames(df) <- c("FID", "IID", "Age", "Sex", "Batch")
df <- df %>% mutate(Batch = replace(Batch, Batch < 0, 0)) %>%
             mutate(Batch = replace(Batch, Batch > 0, 1))
```

The **PC** version additionally joins the 40 common and 40 rare principal components from the [PCA step](03-pca-correction), prefixing the rare PCs with `R` to keep the column names distinct ([`format_regenie_covars_PCs.R`](https://github.com/jgblanc/quantifying-susceptibility-PGS.github.io/blob/master/scripts/run_gwas/format_regenie_covars_PCs.R)):

```r
# Rename rare PCs R1..R40 so they don't collide with common PC1..PC40
names <- colnames(dfRare)[3:42]
colnames(dfRare)[3:42] <- paste0("R", names)
dfPCs <- cbind(dfCommon, dfRare[, 3:42])
dfOut <- inner_join(df, dfPCs)
```

```
rule format_covariates_PCs:
    input:
        info = "output/phenotypes/participant_covariates_all.csv",
        pcs_common = "output/PCA/commonPCA_{gwas}.eigenvec",
        pcs_rare = "output/PCA/rarePCA_{gwas}.eigenvec"
    output:
        "output/phenotypes/ukb_covariates_{gwas}_pcs.txt"
    shell:
        """
        Rscript code/run_GWAS/format_regenie_covars_PCs.R {input.info} {input.pcs_common} {input.pcs_rare} {output}
        """
```

---

## Step 2: Whole-genome model (regenie step 1)

Regenie's first step fits a whole-genome ridge-regression model on the QC-passed, directly-genotyped calls (`ukb_cal_allChrs` restricted to `qc_pass.snplist`). This model accounts for relatedness and polygenic background in step 2, and produces the per-phenotype prediction list (`ukb_step1_QT_pred.list`) that the step-2 rules take as input via `--pred`. Step 1 is run separately for each GWAS panel and correction setting:

```
rule run_step_one:
    input:
        bed = "data/genotypes/ukb_cal_allChrs.bed",
        bim = "data/genotypes/ukb_cal_allChrs.bim",
        fam = "data/genotypes/ukb_cal_allChrs.fam",
        snplist = "data/genotypes/qc_pass.snplist",
        keep = "data/ids/gwas_ids/{gwas}.txt",
        phenotype = "data/regenie_files/ukb_phenotypes_QT.txt",
        covar = "data/regenie_files/ukb_covariates_{gwas}_{covars}.txt"
    output:
        "data/step1-results/{gwas}/{covars}/ukb_step1_QT_pred.list"
    params:
        bed_prefix = "data/genotypes/ukb_cal_allChrs",
        out_prefix = "data/step1-results/{gwas}/{covars}/ukb_step1_QT",
        lowmem_prefix = "regenie_tmp_preds_{gwas}_{covars}"
    threads: 16
    resources:
        mem_mb=38000,
        time="24:00:00"
    shell:
        """
        regenie --step 1 \
                --bed {params.bed_prefix} \
                --extract {input.snplist} \
                --keep {input.keep} \
                --phenoFile {input.phenotype} \
                --covarFile {input.covar} \
                --bsize 1000 \
                --lowmem \
                --lowmem-prefix {params.lowmem_prefix} \
                --out {params.out_prefix}
        """
```

---

## Step 3: Association tests (regenie step 2)

Step 2 runs the per-variant association tests. The two models differ only in how they use the step-1 predictions: the **linear-regression** model passes `--ignore-pred` (no random effect), while the **linear mixed model** uses `--pred` to include the whole-genome model as a random effect.

```
rule run_step_two_LR:
    input:
        pgen = ".../ukb_imp_chr{chr}_v3.pgen",
        covar = "data/regenie_files/ukb_covariates_{gwas}_{covars}.txt",
        phenotype = "data/regenie_files/ukb_phenotypes_QT.txt",
        pred_list = "data/step1-results/{gwas}/{covars}/ukb_step1_QT_pred.list"
    output:
        "data/step2-results/{gwas}/{covars}/LR/raw/chr{chr}_Standing_Height.regenie",
        # ... one output per phenotype ...
    shell:
        """
        regenie --step 2 \
                --pgen {params.plink_prefix} \
                --ref-first \
                --ignore-pred \
                --phenoFile {input.phenotype} \
                --covarFile {input.covar} \
                --pThresh 0.01 \
                --pred {input.pred_list} \
                --bsize 400 \
                --threads 8 \
                --out {params.regenie_prefix}
        """
```

The LMM rule (`run_step_two_LMM`) is identical except that it **omits `--ignore-pred`**, so the whole-genome predictions enter as a random effect:

```
        regenie --step 2 \
                --pgen {params.plink_prefix} \
                --ref-first \
                --phenoFile {input.phenotype} \
                --covarFile {input.covar} \
                --pThresh 0.01 \
                --pred {input.pred_list} \
                --bsize 400 \
                --threads 8 \
                --out {params.regenie_prefix}
```

Both are run for every chromosome, phenotype, GWAS panel, and `covars` setting (`pcs`/`nopcs`), writing raw per-chromosome regenie output to `data/step2-results/{gwas}/{covars}/{LR,LMM}/raw/`.

---

## Step 4: Format summary statistics

Finally, the raw regenie output is harmonized into a consistent SNP/effect table — selecting the effect (`BETA`) and significance (`LOG10P`) columns and standardizing the header — by [`format_regenie_SS.R`](https://github.com/jgblanc/quantifying-susceptibility-PGS.github.io/blob/master/scripts/run_gwas/format_regenie_SS.R):

```r
dfSS <- dfSS %>% select("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "BETA", "LOG10P")
colnames(dfSS) <- c("#CHROM", "POS", "rsID", "REF", "ALT", "BETA", "LOG10_P")
dfOut <- inner_join(dfSNPs, dfSS) %>%
  select("#CHROM", "POS", "ID", "REF", "ALT", "BETA", "LOG10_P")
```

These harmonized summary statistics are the input to the clumping-and-thresholding step on the [next page](06-pgs-association-tests).

---

## Scripts to copy into the website repo

To make the script links on this page resolve, copy the following files from the [original code repository](https://github.com/jgblanc/strat2) into the website repo, preserving the subfolder layout:

| Copy from `strat2` | To website repo |
|--------------------|-----------------|
| `code/run_gwas/format_regenie_phenotype.R` | `scripts/run_gwas/format_regenie_phenotype.R` |
| `code/run_gwas/format_regenie_covars_no_PCs.R` | `scripts/run_gwas/format_regenie_covars_no_PCs.R` |
| `code/run_gwas/format_regenie_covars_PCs.R` | `scripts/run_gwas/format_regenie_covars_PCs.R` |
| `code/run_gwas/format_regenie_SS.R` | `scripts/run_gwas/format_regenie_SS.R` |

> **Note on the directory name.** The scripts live in `code/run_gwas/` (lowercase) in the repository, but the `snakefile_main` rules call them as `code/run_GWAS/...` (capitalized). On a case-sensitive filesystem (e.g. the Linux HPC cluster) this mismatch will break the rule — rename one or the other so they agree. I've used the lowercase `run_gwas/` path in the copy table above to match what's actually in the repo.

> **Note on regenie step 1.** The step-1 rule shown above is adapted from `run_step_one` in `snakefile_DNANexus`, reformatted as a plain Snakemake rule (the DNANexus `dx run`/`dx upload` wrapping and file-staging removed). It also depends on two upstream prep steps from that snakefile — merging the per-chromosome genotype calls into `ukb_cal_allChrs` and writing `qc_pass.snplist` (the `--maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15` QC pass) — which you'd similarly need to reformat if you want them documented here.
