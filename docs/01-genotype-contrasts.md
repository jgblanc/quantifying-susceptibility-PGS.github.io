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


For each prediction panel, we first construct a standardized test vector $t$ that positions each individual along the ancestry axis of interest, then compute the allele frequency contrast $\hat{r} = \frac{1}{N} X^{\top} t$ across all overlapping SNPs. Contrasts are computed separately for each chromosome and stored as `.rvec` files.

---

#### All Continental Pairwise Contrasts

Test vectors for all $\binom{5}{2} = 10$ pairwise contrasts are constructed simultaneously using [`get_Tvec_all.R`](https://github.com/jgblanc/quantifying-susceptibility-PGS.github.io/blob/master/scripts/HGDP1KG/get_Tvec_all.R), which codes each individual as $+1$ or $-1$ depending on their population label, then mean-centers and standardizes the vector to unit variance:

```
rule make_test_vectors_all:
    input:
        "data/HGDP1KG/ids/test_ids/all.txt"
    output:
        p1="data/HGDP1KG/TestVecs/eas-nfe.txt",
        p2="data/HGDP1KG/TestVecs/eas-sas.txt",
        p3="data/HGDP1KG/TestVecs/eas-afr.txt",
        p4="data/HGDP1KG/TestVecs/eas-amr.txt",
        p5="data/HGDP1KG/TestVecs/nfe-sas.txt",
        p6="data/HGDP1KG/TestVecs/nfe-afr.txt",
        p7="data/HGDP1KG/TestVecs/nfe-amr.txt",
        p8="data/HGDP1KG/TestVecs/sas-afr.txt",
        p9="data/HGDP1KG/TestVecs/sas-amr.txt",
        p10="data/HGDP1KG/TestVecs/afr-amr.txt"
    shell:
        """
        Rscript code/HGDP1KG/get_Tvec_all.R {output.p1} ... {output.p10} {input}
        """
```

The allele frequency contrasts are then computed for each chromosome using [`get_contrasts_all.R`](https://github.com/jgblanc/quantifying-susceptibility-PGS.github.io/blob/master/scripts/HGDP1KG/get_contrasts_all.R):

```
rule compute_r_all:
    input:
        p1="data/HGDP1KG/TestVecs/eas-nfe.txt",
        ...
        p10="data/HGDP1KG/TestVecs/afr-amr.txt",
        SNPs="data/HGDP1KG/variants/{gwas}/all/overlappingSNPs_chr{chr}.txt"
    output:
        p1="data/HGDP1KG/r/{gwas}/eas-nfe_chr{chr}.rvec",
        ...
        p10="data/HGDP1KG/r/{gwas}/afr-amr_chr{chr}.rvec"
    params:
        prefix_tp="data/HGDP1KG/plink2-files/gnomad.genomes.v3.1.2.hgdp_tgp.chr{chr}"
    shell:
        """
        Rscript code/HGDP1KG/get_contrasts_all.R {input.SNPs} {params.prefix_tp} \
        {input.p1} ... {input.p10} {output.p1} ... {output.p10}
        """
```

The overlapping SNP list is then copied to each pairwise contrast directory:

```
rule rename_overlap_snps_all:
    input:
        "data/HGDP1KG/variants/{gwas}/all/overlappingSNPs_chr{chr}.txt"
    output:
        p1="data/HGDP1KG/variants/{gwas}/eas-nfe/overlappingSNPs_chr{chr}.txt",
        ...
        p10="data/HGDP1KG/variants/{gwas}/afr-amr/overlappingSNPs_chr{chr}.txt"
    shell:
        """
        cp {input} {output.p1}
        ...
        cp {input} {output.p10}
        """
```

---

#### Sardinia vs. Mainland Europe

The Sardinia contrast test vector is constructed using [`get_Tvec_sdi.R`](https://github.com/jgblanc/quantifying-susceptibility-PGS.github.io/blob/master/scripts/HGDP1KG/get_Tvec_sdi.R), which codes Sardinian individuals as $+1$ and all other non-Finnish Europeans as $-1$ (before standardization):

```
rule make_test_vectors_sdi:
    input:
        "data/HGDP1KG/ids/test_ids/nfe.txt"
    output:
        "data/HGDP1KG/TestVecs/sdi-eur.txt"
    shell:
        """
        Rscript code/HGDP1KG/get_Tvec_sdi.R {input} {output}
        """
```

The contrast is then computed using [`get_contrasts_sdi.R`](https://github.com/jgblanc/quantifying-susceptibility-PGS.github.io/blob/master/scripts/HGDP1KG/get_contrasts_sdi.R):

```
rule compute_r_sdi:
    input:
        p1="data/HGDP1KG/TestVecs/sdi-eur.txt",
        SNPs="data/HGDP1KG/variants/{gwas}/all/overlappingSNPs_chr{chr}.txt"
    output:
        p1="data/HGDP1KG/r/{gwas}/sdi-eur_chr{chr}.rvec"
    params:
        prefix_tp="data/HGDP1KG/plink2-files/gnomad.genomes.v3.1.2.hgdp_tgp.chr{chr}"
    shell:
        """
        Rscript code/HGDP1KG/get_contrasts_sdi.R {input.SNPs} {params.prefix_tp} \
        {input.p1} {output.p1}
        """
```

---

#### Latitude and Longitude in Eurasia

Test vectors for Eurasian latitude and longitude are constructed using [`get_Tvec_cord.R`](https://github.com/jgblanc/quantifying-susceptibility-PGS.github.io/blob/master/scripts/HGDP1KG/get_Tvec_cord.R), which mean-centers and standardizes the provided geographic coordinates to unit variance:

```
rule make_test_vectors_Eurasia:
    input:
        testids="data/HGDP1KG/ids/test_ids/eurasia.txt"
    output:
        Lat="data/HGDP1KG/TestVecs/eurasia-lat.txt",
        Long="data/HGDP1KG/TestVecs/eurasia-long.txt"
    shell:
        """
        Rscript code/HGDP1KG/get_Tvec_cord.R {input.testids} {output.Lat} {output.Long}
        """
```

Contrasts are computed using [`get_contrasts_cord.R`](https://github.com/jgblanc/quantifying-susceptibility-PGS.github.io/blob/master/scripts/HGDP1KG/get_contrasts_cord.R):

```
rule compute_r_eurasia:
    input:
        p1="data/HGDP1KG/TestVecs/eurasia-lat.txt",
        p2="data/HGDP1KG/TestVecs/eurasia-long.txt",
        SNPs="data/HGDP1KG/variants/{gwas}/eurasia/overlappingSNPs_chr{chr}.txt"
    output:
        p1="data/HGDP1KG/r/{gwas}/eurasia-lat_chr{chr}.rvec",
        p2="data/HGDP1KG/r/{gwas}/eurasia-long_chr{chr}.rvec"
    params:
        prefix_tp="data/HGDP1KG/plink2-files/gnomad.genomes.v3.1.2.hgdp_tgp.chr{chr}"
    shell:
        """
        Rscript code/HGDP1KG/get_contrasts_cord.R {input.SNPs} {params.prefix_tp} \
        {input.p1} {output.p1} {input.p2} {output.p2}
        """
```

---

#### Latitude and Longitude in Europe

The same scripts are used for the non-Finnish European latitude and longitude contrasts, using the CEU-excluded sample:

```
rule get_TestVector_Eur:
    input:
        testids="data/HGDP1KG/ids/test_ids/nfe_no_ceu.txt"
    output:
        Lat="data/HGDP1KG/TestVecs/eur-lat.txt",
        Long="data/HGDP1KG/TestVecs/eur-long.txt"
    shell:
        """
        Rscript code/HGDP1KG/get_Tvec_cord.R {input.testids} {output.Lat} {output.Long}
        """

rule compute_r_eur:
    input:
        p1="data/HGDP1KG/TestVecs/eur-lat.txt",
        p2="data/HGDP1KG/TestVecs/eur-long.txt",
        SNPs="data/HGDP1KG/variants/{gwas}/nfe_no_ceu/overlappingSNPs_chr{chr}.txt"
    output:
        p1="data/HGDP1KG/r/{gwas}/eur-lat_chr{chr}.rvec",
        p2="data/HGDP1KG/r/{gwas}/eur-long_chr{chr}.rvec"
    params:
        prefix_tp="data/HGDP1KG/plink2-files/gnomad.genomes.v3.1.2.hgdp_tgp.chr{chr}"
    shell:
        """
        Rscript code/HGDP1KG/get_contrasts_cord.R {input.SNPs} {params.prefix_tp} \
        {input.p1} {output.p1} {input.p2} {output.p2}
        """
```
---

## Country of Birth Contrasts

### Data Preparation

### Computing Contrasts
