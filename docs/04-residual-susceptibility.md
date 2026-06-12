---
layout: default
title: Residual Susceptibility
---

<a href="{{ '/' | relative_url }}" class="btn">← Back to Home</a>

# Residual Susceptibility

After computing the [principal components](03-pca-correction), we can ask the two questions that determine whether PCA has actually protected a polygenic score along a given target axis:

1. **How much susceptibility remains?** The residual susceptibility $H'$ is the proportion of GWAS-panel genetic variance still aligned with the target axis *after* projecting out the PCs. It is the absolute measure of residual stratification risk.
2. **How effectively did the PCs capture the axis?** The efficacy statistic $V_K$ is the fraction of the target-axis variance captured by the top $K$ PCs.

The residual susceptibility is defined exactly like $H$, but using the PC-corrected genotypes $G' = (I - UU^\top)\,G$:

$$H' = \text{Cor}^2(G', r) = \sigma_{f'}^2 \cdot \sigma_{r}^2$$

where $f' = (I - UU^\top)\,f$ is the residual target axis. The efficacy statistic is then

$$V_K = 1 - \frac{H'}{H} = 1 - \frac{\sigma_{f'}^2}{\sigma_f^2}.$$

---

## Overview of the estimator

We estimate $H'$ and $V_K$ in three pieces:

1. **Residualize the genotypes against the PCs.** Each SNP's dosages are regressed on the principal components and replaced by the residuals, yielding $G'$. To avoid overfitting, SNPs on each chromosome are residualized against PCs estimated from the *opposite* half of the genome (odd vs. even chromosomes).
2. **Estimate $H'$.** We project the residualized genotypes onto the contrast $\hat{r}$, form $\hat{H}' = \hat{\sigma}_{f'}^2 \cdot \hat{\sigma}_{r}^2$, and test it against the same noise floor $1/L$ using a block jackknife.
3. **Estimate efficacy $V_K$.** Using an even/odd cross-validation scheme, we compute the noise-corrected fraction of target-axis variance captured by the top $K$ common, rare, or combined PCs.

---

## Step 1: Residualize genotypes against the PCs

We first extract per-chromosome dosages for the PC SNP set, then residualize each SNP against the principal components. The dosage extraction:

```
rule dosage_chr_snps:
    input:
        IDs="data/ids/gwas_ids/{gwas}.txt"
    output:
        "/scratch/jgblanc/ukbb/dosages/{gwas}/chr{chr}_dosage.traw"
    params:
        plink_prefix="/scratch/jgblanc/ukbb/plink2-files/ALL/ukb_pcSNPs_ALL_v3",
        out_prefix="/scratch/jgblanc/ukbb/dosages/{gwas}/chr{chr}_dosage"
    shell:
        """
        plink2 --pfile {params.plink_prefix} \
        --keep {input.IDs} \
        --chr {wildcards.chr} \
        --recode A-transpose \
        --out {params.out_prefix}
        """
```

The residualization is handled by [`residualize.R`](https://github.com/jgblanc/quantifying-susceptibility-PGS.github.io/blob/master/scripts/calculate_HPrime/residualize.R), which uses **odd-chromosome PCs to residualize even-chromosome SNPs and vice versa**, so that the PCs used for correction are never estimated from the SNPs being corrected:

```
rule residulize_snps:
    input:
        SNPs="/scratch/jgblanc/ukbb/dosages/{gwas}/chr{chr}_dosage.traw",
        oddPCs="output/calculate_PCA/{gwas}/odd_PCA.eigenvec",
        evenPCs="output/calculate_PCA/{gwas}/even_PCA.eigenvec",
        oddPCs_rare="output/calculate_PCA/{gwas}/rare_odd_PCA.eigenvec",
        evenPCs_rare="output/calculate_PCA/{gwas}/rare_even_PCA.eigenvec"
    output:
        traw="/scratch/jgblanc/ukbb/dosages/{gwas}/chr{chr}_residual.traw"
    shell:
        """
        Rscript code/calculate_HPrime/residualize.R {input.SNPs} {params.pca_prefix} \
        {wildcards.chr} {output.traw}
        rm {input.SNPs}
        """
```

Internally, common and rare PCs are concatenated into a single covariate matrix (plus an intercept), and each SNP's dosages are residualized via a precomputed QR decomposition for speed:

```r
# Choose the opposite half-genome PCs for this chromosome
if (chr_num %in% c(1,3,5,...,21)) {
  pca_file_path_common <- paste0(pc_prefix, "even_PCA.eigenvec")
  pca_file_path_rare   <- paste0(pc_prefix, "rare_even_PCA.eigenvec")
} else {
  pca_file_path_common <- paste0(pc_prefix, "odd_PCA.eigenvec")
  pca_file_path_rare   <- paste0(pc_prefix, "rare_odd_PCA.eigenvec")
}

# Combine common + rare PCs into one covariate matrix (with intercept)
covars_with_intercept <- cbind(1, covars)
qr_covars <- qr(covars_with_intercept)

# Residualize each SNP's dosages and standardize
fitted_vals <- qr.fitted(qr_covars, dosages)
resids <- (dosages - fitted_vals) / sd(dosages)
```

---

## Step 2: Estimating $\hat{H}'$

With residualized genotypes in hand, $\hat{H}'$ is computed exactly as $\hat{H}$ was on the [previous page](02-uncorrected-susceptibility): project $G'$ onto the contrast block-by-block, form $\hat{H}' = \hat{\sigma}_{f'}^2 \cdot \hat{\sigma}_{r}^2$, and obtain a block-jackknife standard error.

```
rule calcPrime_H_Jacknife:
    input:
        r=expand("data/{dataset}/r/{gwas}/{contrasts}_chr{chr}.rvec", chr=CHR),
        IDs="data/ids/gwas_ids/{gwas}.txt",
        snps="data/ukbb/ukbb_pc_snps_block.txt",
        resids=expand("/scratch/jgblanc/ukbb/dosages/{gwas}/chr{chr}_residual.traw", chr=CHR)
    output:
        H = "output/calculate_HPrime/{dataset}/{gwas}/H_BJ_{contrasts}.txt"
    shell:
        """
        Rscript code/calculate_HPrime/calc_H.R {params.resid_prefix} {params.r_prefix} \
        {input.snps} {input.IDs} {output.H}
        """
```

[`calculate_HPrime/calc_H.R`](https://github.com/jgblanc/quantifying-susceptibility-PGS.github.io/blob/master/scripts/calculate_HPrime/calc_H.R) reads the residualized `.traw` files block-by-block, projects them onto the standardized contrast, and accumulates the residual target axis:

```r
# Standardize the contrast
dfALL$r <- dfALL$r / sd(dfALL$r)

# For each block, project residualized genotypes onto r
matBlock <- t(as.matrix(dfBlock[, 4:ncol(dfBlock)]))
matBlock <- scale(matBlock)
matf <- matBlock %*% as.matrix(dfR_block$r)
dfFGr_mat[, counter] <- apply(matf, 1, sum)

# Accumulate and scale
FGr <- apply(dfFGr_mat, 1, sum) * (1 / sqrt(L - 1))
H   <- (1 / (M * (L - 1))) * (t(FGr) %*% FGr)
```

As with $\hat{H}$, significance is assessed by a one-sided test against the detection limit $1/L$, using the block-jackknife variance:

```r
pvalNorm <- pnorm(H, mean = 1/L, sd = sqrt(varH), lower.tail = FALSE)
```

Because $\hat{H}'$ is typically driven down close to the detection limit, a non-significant result indicates only that residual structure is *below the noise floor* — not that it is exactly zero. Results across contrasts are concatenated into `plots/overlap_stats_prime/H_{dataset}.txt`.

<!-- ============================================================
FIGURE PLACEHOLDER — could not be retrieved from Overleaf (login-gated).
Suggested: paper Fig. 4 — H' after correction across the six panels,
showing that correction flattens the disparity between diverse and
homogeneous panels (uncorrected H shown as semi-transparent points,
dashed line = 1/L detection limit).
Export to assets/images/ and fill in the filename + caption.

![Residual susceptibility H' across GWAS panels](../assets/images/FIGURE_FILENAME.png)
============================================================ -->

---

## Step 3: PCA efficacy ($V_K$)

To quantify how much of the target axis the PCs capture — without the result being inflated by estimation noise in $\hat{f}$ — we use an **even/odd cross-validation** scheme. PCs estimated on odd chromosomes are used to predict the target axis projected from even-chromosome SNPs (and vice versa). This is implemented in [`R2_EO.R`](https://github.com/jgblanc/quantifying-susceptibility-PGS.github.io/blob/master/scripts/calculate_FGr/R2_EO.R) (common or rare PCs) and [`R2_EO_both.R`](https://github.com/jgblanc/quantifying-susceptibility-PGS.github.io/blob/master/scripts/calculate_FGr/R2_EO_both.R) (common + rare PCs together):

```
rule calc_R2_even:
    input:
        PCs="output/calculate_PCA/{gwas}/odd_PCA.eigenvec",
        FGr="output/calculate_FGr/{dataset}/{gwas}/FGrMat_{contrasts}.txt",
        SNPs="output/calculate_FGr/{dataset}/{gwas}/SNPNum_{contrasts}.txt"
    output:
        R2="output/calculate_FGr/{dataset}/{gwas}/R2_EO/{contrasts}.txt"
    shell:
        """
        Rscript code/calculate_FGr/R2_EO.R {input.PCs} {input.FGr} {input.SNPs} {output.R2}
        """
```

The script first estimates the **noise-corrected signal** in $\hat{f}$ — the fraction of its variance due to real ancestry variation rather than estimation noise ($R^2_f$ in the manuscript) — using a leave-one-block-out jackknife to estimate the error variance:

```r
# Leave-one-block-out estimate of the estimation-error variance in fhat
for (i in 1:numBlocks) {
  fhat_i <- calc_fhat(dfMat[, -i], dfR_not_i$r)
  jckFGr[, i] <- ((L - mi) / mi) * (fhat - fhat_i)^2
}
error  <- mean(apply(jckFGr, 1, mean)) / var(fhat)
signal <- 1 - error
```

It then regresses $\hat{f}$ on the top $K$ cross-validated PCs and reports the cumulative $R^2$, and the ratio $R^2 / \text{signal}$ — the noise-corrected estimate of the captured fraction $\hat{V}_K$:

```r
for (i in seq_len(ncol(PC_nums))) {
  mod   <- lm(fhat ~ PC_nums[, 1:i])
  R2    <- summary(mod)$r.squared
  Ratio <- R2 / signal          # noise-corrected efficacy V_K
  dfOut[i, ] <- c(i, H, varH, pvalNorm, signal, w, R2, Ratio)
}
```

The same rule is run with rare PCs (`calc_R2_even_rare`, also via `R2_EO.R`) and with common + rare PCs jointly (`calc_R2_even_common_rare`, via `R2_EO_both.R`), letting us compare how much additional structure rare variants capture beyond common ones.

<!-- ============================================================
FIGURE PLACEHOLDER — could not be retrieved from Overleaf (login-gated).
Suggested: paper Fig. 3 — V_K (fraction of target-axis variance captured)
for common, rare, and combined PCs, and the cumulative V_K across PC index.
Export to assets/images/ and fill in the filename + caption.

![Fraction of target-axis variance captured by common and rare PCs](../assets/images/FIGURE_FILENAME.png)
============================================================ -->

---

## Scripts to copy into the website repo

To make the script links on this page resolve, copy the following files from the [original code repository](https://github.com/jgblanc/strat2) into the website repo, preserving the subfolder layout:

| Copy from `strat2` | To website repo |
|--------------------|-----------------|
| `code/calculate_HPrime/residualize.R` | `scripts/calculate_HPrime/residualize.R` |
| `code/calculate_HPrime/calc_H.R` | `scripts/calculate_HPrime/calc_H.R` |
| `code/calculate_FGr/R2_EO.R` | `scripts/calculate_FGr/R2_EO.R` |
| `code/calculate_FGr/R2_EO_both.R` | `scripts/calculate_FGr/R2_EO_both.R` |

> **Note.** `calculate_HPrime/calc_H.R` shares its filename with the `H` estimator on page 02 (`calculate_H/calc_H.R`) but is a distinct script that reads PC-residualized genotypes — keep them in separate `scripts/calculate_H/` and `scripts/calculate_HPrime/` folders so the links don't collide. The `calcPrime_FGrMat`/`calcPrime_FGr` rules in `snakefile_main2` also reference `calc_FGr_prime.R` and `calc_fhat_sd.R` for producing a residualized $\hat{f}'$ for visualization; add those to `scripts/calculate_FGr/` if you want to document that path as well.

> **Reminder.** As with pages 02–03, inline math underscores (e.g. `$\sigma_f^2$`) need the same kramdown fix you applied earlier; the display equations (`$$ … $$` on their own lines) render fine as-is.
