---
layout: default
title: Genotype Contrasts
---

# Genotype Contrasts

Genotype contrasts capture the systematic differences in allele frequencies between groups or along a continuous axis of ancestry. For each prediction panel, we compute a contrast vector $\hat{r}$ that summarizes how allele frequencies vary along the target axis. This vector is then used in downstream steps to estimate the susceptibility of polygenic scores to stratification bias.

Formally, the contrast vector is defined as:

$$\hat{r} = \frac{1}{N} X^{\top} t$$

where $X$ is the $N \times L$ matrix of prediction panel genotypes and $t$ is the standardized test vector positioning each individual along the ancestry axis of interest.

---

## HGDP1kGP Contrasts

### Data Preparation

## HGDP1kGP Contrasts

### Data Preparation

The HGDP1kGP data were distributed in plink2 format aligned to the hg19 genome build. Before computing genotype contrasts, we harmonized the HGDP1kGP data to match the UK Biobank reference/alternate allele coding, and identified the set of overlapping SNPs between each prediction panel and GWAS panel.

#### Step 1: Recode SNPs to match UK Biobank allele coding

Because the HGDP1kGP and UK Biobank data were processed independently, the reference and alternate alleles at some SNPs may be flipped between the two datasets. We first computed allele frequencies in the HGDP1kGP data for each chromosome:

```
rule HGDP_freq:
    input:
        pgen="/gpfs/data/berg-lab/data/HGDP1KG/plink2-files-hg19/gnomad.genomes.v3.1.2.hgdp_tgp.chr{chr}.pgen"
    output:
        freq="data/HGDP1KG/variants-raw/gnomad.genomes.v3.1.2.hgdp_tgp.chr{chr}.afreq"
    shell:
        """
        plink2 --pfile {params.prefix_in} \
        --freq \
        --threads 8 \
        --memory 38000 \
        --out {params.prefix_out}
        """
```

We then compared allele frequencies between the HGDP1kGP and UK Biobank to identify SNPs that needed to be flipped using [`flip_snps.R`](https://github.com/jgblanc/quantifying-susceptibility-PGS.github.io/blob/master/scripts/HGDP1KG/flip_snps.R):

```
rule recode_list:
    input:
        freq_hgdp="data/HGDP1KG/variants-raw/gnomad.genomes.v3.1.2.hgdp_tgp.chr{chr}.afreq",
        freq_ukbb="data/ukbb/variants/ukb_imp_chr{chr}_v3.afreq"
    output:
        "data/HGDP1KG/variants-fliped/flipped_snps_{chr}.txt"
    shell:
        """
        Rscript code/HGDP1KG/flip_snps.R {input.freq_ukbb} {input.freq_hgdp} {output}
        """
```

The identified SNPs were then recoded so that the HGDP1kGP data uses the same reference/alternate allele orientation as the UK Biobank:

```
rule HGDP_recode:
    input:
        snp_list="data/HGDP1KG/variants-fliped/flipped_snps_{chr}.txt"
    output:
        pgen="data/HGDP1KG/plink2-files/gnomad.genomes.v3.1.2.hgdp_tgp.chr{chr}.pgen"
    shell:
        """
        plink2 --pfile {params.prefix_in} \
        --extract {input.snp_list} \
        --ref-allele force {input.snp_list} \
        --make-pgen \
        --out {params.prefix_out}
        """
```

#### Step 2: Identify overlapping SNPs

For each combination of prediction panel and GWAS panel, we identified the set of SNPs with greater than 1% minor allele frequency in both panels. We first computed allele frequencies within each prediction panel subdataset:

```
rule get_AF_Test:
    input:
        test="data/HGDP1KG/ids/test_ids/{subdataset}.txt"
    output:
        test="data/HGDP1KG/variantFreq/{subdataset}_{chr}.afreq"
    shell:
        """
        plink2 --pfile {params.prefix_in} \
        --freq \
        --keep {input.test} \
        --threads 8 \
        --memory 38000 \
        --out {params.prefix_out_test}
        """
```

We then identified the overlapping SNPs between each prediction panel and GWAS panel using [`overlapping_snps.R`](https://github.com/jgblanc/quantifying-susceptibility-PGS.github.io/blob/master/scripts/HGDP1KG/overlapping_snp.R):

```
rule get_overlapping_snps:
    input:
        freq_test="data/HGDP1KG/variantFreq/{subdataset}_{chr}.afreq",
        freq_gwas="data/ukbb/variantFreq/{gwas}_{chr}.afreq"
    output:
        "data/HGDP1KG/variants/{gwas}/{subdataset}/overlappingSNPs_chr{chr}.txt"
    shell:
        """
        Rscript code/HGDP1KG/overlapping_snps.R {input.freq_gwas} {input.freq_test} {output}
        """
```

### Computing Contrasts

---

## Country of Birth Contrasts

### Data Preparation

### Computing Contrasts
