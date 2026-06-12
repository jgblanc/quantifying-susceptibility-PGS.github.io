---
layout: default
title: PCA Correction
---

[← Back to Home]({{ '/' | relative_url }})

# PCA Correction

The [uncorrected susceptibility](02-uncorrected-susceptibility) $H$ measures how much GWAS-panel genetic variation aligns with a target ancestry axis *before* any attempt to control for population structure. The standard tool for that control is **Principal Component Analysis (PCA)**: we estimate the top axes of genetic variation in the GWAS panel and include them as fixed-effect covariates, which is geometrically equivalent to projecting the phenotype onto the subspace orthogonal to those axes.

This page documents how we compute the principal components used for correction. We estimate two sets of PCs that capture structure at different allele-frequency scales:

- **Common variant PCs** — computed from the LD-pruned, high-quality common SNP set, capturing broad/deep population structure.
- **Rare variant PCs** — computed from LD-pruned variants with $\text{MAF} < 0.01$, capturing finer-scale, more recent structure that common variants miss.

These PCs feed directly into the estimation of [residual susceptibility H' and PCA efficacy](04-residual-susceptibility) and into the [GWAS](05-gwas) as covariates.

---

## Common variant PCA

Common-variant PCs are computed on the merged, LD-pruned high-quality SNP set (the 147,604 UK Biobank PC SNPs; see [Uncorrected Susceptibility](02-uncorrected-susceptibility)). We extract the top 40 PCs using plink2's randomized (`approx`) algorithm, which is accurate and tractable at biobank scale:

```
rule calc_PCA:
    input:
        IDs="data/ids/gwas_ids/{gwas}.txt"
    output:
        "output/calculate_PCA/{gwas}/commonPCA.eigenvec",
        "output/calculate_PCA/{gwas}/commonPCA.eigenval"
    params:
        plink_prefix="/scratch/jgblanc/ukbb/plink2-files/ALL/ukb_pcSNPs_ALL_v3",
        out_prefix="output/calculate_PCA/{gwas}/commonPCA"
    threads: 16
    resources:
        mem_mb=60000,
        time="36:00:00"
    shell:
        """
        plink2 --pfile {params.plink_prefix} \
        --keep {input.IDs} \
        --threads 16 \
        --pca 40 approx \
        --memory 60000 \
        --out {params.out_prefix}
        """
```

PCs are estimated **separately within each GWAS panel** (`{gwas}` = `WBS`, `ALL`, `e5`, `e10`, `e50`, `e200`), so that each panel is corrected using axes estimated from its own genetic composition.

---

## Rare variant PCA

Rare-variant PCs require a dedicated SNP set. We extract biallelic SNPs with $\text{MAF} < 0.01$ directly from the imputed bgen files, applying the same quality filters as the common set and LD-pruning with `--indep-pairwise 1000 80 0.1`:

```
rule extract_rare_variants:
    input:
        bgen=".../ukb_imp_chr{chr}_v3.bgen",
        sample=".../ukb22828_c22_b0_v3_s487192.sample",
        ids="data/ids/gwas_ids/{gwas}.txt"
    output:
        ".../RareVariants/{gwas}/rare_variants_{chr}.prune.in",
        ".../RareVariants/{gwas}/rare_variants_{chr}.pgen"
    shell:
        """
        plink2 --bgen {input.bgen} ref-first \
        --sample {input.sample} \
        --mind 0.1 --geno 0.1 --max-maf 0.01 \
        --rm-dup exclude-all --snps-only --max-alleles 2 \
        --keep {input.ids} \
        --indep-pairwise 1000 80 0.1 \
        --make-pgen --set-all-var-ids @:# \
        --threads 16 --memory 38000 \
        --out {params.prefix}
        """
```

The per-chromosome pruned variants are then merged into a single file:

```
rule concat_chr:
    input:
        pgen=expand(".../RareVariants/{gwas}/rare_variants_{chr}.pgen", chr=CHR),
        snps=expand(".../RareVariants/{gwas}/rare_variants_{chr}.prune.in", chr=CHR)
    output:
        ".../RareVariants/{gwas}/pcSNPs_ALL.pgen"
    shell:
        """
        cat {input.snps} > snplist_{wildcards.gwas}.txt
        plink2 --pmerge-list tmp_chr_list_{wildcards.gwas}.txt \
        --extract snplist_{wildcards.gwas}.txt \
        --make-pgen \
        --out {params.prefix_out}
        """
```

Finally, we subsample to 1,000,000 variants (`--thin-count`) and extract the top 40 rare-variant PCs:

```
rule calc_PCA_rare:
    input:
        IDs="data/ids/gwas_ids/{gwas}.txt"
    output:
        "output/calculate_PCA/{gwas}/rarePCA.eigenvec",
        "output/calculate_PCA/{gwas}/rarePCA.eigenval"
    params:
        plink_prefix=".../RareVariants/{gwas}/pcSNPs_ALL",
        out_prefix="output/calculate_PCA/{gwas}/rarePCA"
    shell:
        """
        plink2 --pfile {params.plink_prefix} \
        --keep {input.IDs} \
        --thin-count 1000000 \
        --pca 40 approx \
        --memory 400000 \
        --out {params.out_prefix}
        """
```
---

## Even/odd chromosome PCs (for cross-validation)

To estimate how well the PCs capture the target axis without overfitting (the efficacy statistic $V_K$ on the [next page](04-residual-susceptibility)), we also compute PCs on **disjoint halves of the genome** — odd chromosomes and even chromosomes — for both common and rare variants. PCs estimated on one half are used to predict the target axis projected from the other half, giving an unbiased estimate of captured variance.

```
rule calc_odd_PCA:
    input:
        IDs="data/ids/gwas_ids/{gwas}.txt"
    output:
        "output/calculate_PCA/{gwas}/odd_PCA.eigenvec"
    params:
        plink_prefix="/scratch/jgblanc/ukbb/plink2-files/ALL/ukb_pcSNPs_ALL_v3",
        out_prefix="output/calculate_PCA/{gwas}/odd_PCA"
    shell:
        """
        plink2 --pfile {params.plink_prefix} \
        --keep {input.IDs} \
        --chr 1,3,5,7,9,11,13,15,17,19,21 \
        --pca 40 approx \
        --memory 60000 \
        --out {params.out_prefix}
        """
```

The `calc_even_PCA`, `calc_odd_PCA_rare`, and `calc_even_PCA_rare` rules are identical except for the chromosome list (`--chr 2,4,...,22`) and, for the rare-variant versions, the rare-variant plink prefix plus `--thin-count 1000000`.

As a further diagnostic, we also estimate PCs from *cumulative* chromosome subsets (chromosome 1, then 1–2, then 1–3, …) to see how quickly captured variance saturates with genomic breadth. This is handled by [`calc_PCA_chr.R`](https://github.com/jgblanc/quantifying-susceptibility-PGS.github.io/blob/master/scripts/calculate_PCA/calc_PCA_chr.R), which loops plink2 over growing chromosome ranges.

---

## Scripts to copy into the website repo

Most PCA steps call plink2 directly inside the Snakemake rules, so only the cumulative-chromosome diagnostic needs an R script:

| Copy from `strat2` | To website repo |
|--------------------|-----------------|
| `code/calculate_PCA/calc_PCA_chr.R` | `scripts/calculate_PCA/calc_PCA_chr.R` |

> **Note.** The rare-variant preparation rules (`extract_rare_variants`, `concat_chr`, the rare PCA rule) live in `snakefile_UKBB`, while the common-variant and even/odd PCA rules live in `snakefile_main2`. The rule names above match those snakefiles; the rule I've labeled `calc_PCA_rare` corresponds to `calc_even_PCA_rare`/`calc_odd_PCA_rare` for the half-genome splits — adjust to a full-genome rare PCA rule if you compute one for the correction covariates directly.
