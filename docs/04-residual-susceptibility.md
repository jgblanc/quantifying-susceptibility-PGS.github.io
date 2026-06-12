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

The key efficiency in the pipeline is that we never re-project residualized genotypes from scratch. Instead, we residualize the **per-block projection matrix** computed for $H$ on the [previous page](02-uncorrected-susceptibility): because projection onto the contrast is linear, residualizing each block's score against the PCs is equivalent to projecting PC-corrected genotypes. We then reuse the *same* $\hat{H}$ estimator on the residualized matrix.

1. **Residualize the FGr matrix against the PCs.** Each block column of the projection matrix is regressed on the principal components and replaced by the residuals, yielding the corrected matrix used for $f'$. To avoid overfitting, blocks on each chromosome are residualized against PCs estimated from the *opposite* half of the genome (odd vs. even chromosomes).
2. **Estimate $H'$.** Run the same `calc_H.R` estimator on the residualized matrix to form $\hat{H}' = \hat{\sigma}_{f'}^2 \cdot \hat{\sigma}_{r}^2$, tested against the noise floor $1/L$ with a block jackknife.
3. **Estimate efficacy $V_K$.** Using an even/odd cross-validation scheme, compute the noise-corrected fraction of target-axis variance captured by the top $K$ common, rare, or combined PCs.

---

## Step 1: Residualize the FGr matrix against the PCs

[`calc_GR_prime.R`](https://github.com/jgblanc/quantifying-susceptibility-PGS.github.io/blob/master/scripts/calculate_FGr/calc_GR_prime.R) takes the per-block projection matrix `FGrMat` (from the [uncorrected step](02-uncorrected-susceptibility)) and residualizes each block column against the principal components, using **odd-chromosome PCs for even-chromosome blocks and vice versa** so that the PCs are never estimated from the blocks they correct:

```
rule residulize_GR:
    input:
        FGr = "output/calculate_FGr/{dataset}/{gwas}/FGrMat_{contrasts}.txt",
        oddPCs="output/calculate_PCA/{gwas}/odd_PCA.eigenvec",
        evenPCs="output/calculate_PCA/{gwas}/even_PCA.eigenvec",
        oddPCs_rare="output/calculate_PCA/{gwas}/rare_odd_PCA.eigenvec",
        evenPCs_rare="output/calculate_PCA/{gwas}/rare_even_PCA.eigenvec",
        snp_num = "output/calculate_FGr/{dataset}/{gwas}/SNPNum_{contrasts}.txt"
    output:
        GRPrime="output/calculate_FGr/{dataset}/{gwas}/FGrMatPrime_{contrasts}.txt"
    params:
        pca_prefix="output/calculate_PCA/{gwas}/"
    shell:
        """
        Rscript code/calculate_FGr/calc_GR_prime.R {input.FGr} {params.pca_prefix} \
        {output.GRPrime} {input.snp_num}
        """
```

Internally, common and rare PCs are joined into one covariate matrix, and each block column is replaced by the residuals of its regression on those PCs:

```r
for (i in seq_along(dfMat)) {
  chr_num <- as.numeric(colnames(dfMat)[i])

  # Use the opposite half-genome PCs for this block's chromosome
  if (chr_num %% 2 == 1) {                      # odd chromosome
    pca_file_path_common <- paste0(pc_prefix, "even_PCA.eigenvec")
    pca_file_path_rare   <- paste0(pc_prefix, "rare_even_PCA.eigenvec")
  } else {                                      # even chromosome
    pca_file_path_common <- paste0(pc_prefix, "odd_PCA.eigenvec")
    pca_file_path_rare   <- paste0(pc_prefix, "rare_odd_PCA.eigenvec")
  }

  dfPCs <- inner_join(dfPCs_common, dfPCs_rare, by = c("FID","IID"))
  covars_df <- dfPCs %>% select(starts_with("PC"))

  y <- as.numeric(dfMat[[i]])
  dfMat_resids[, i] <- resid(lm(y ~ ., data = covars_df))
}
```

The output, `FGrMatPrime`, is the PC-corrected analogue of the projection matrix.

---

## Step 2: Estimating $\hat{H}'$

Because $H'$ has the same form as $H$, it is computed with the **same script** ([`calc_H.R`](https://github.com/jgblanc/quantifying-susceptibility-PGS.github.io/blob/master/scripts/calculate_H/calc_H.R)) — the only change is that the input is the residualized matrix `FGrMatPrime` rather than `FGrMat`:

```
rule calc_HPrime_Jacknife:
    input:
        r=expand("data/{dataset}/r/{gwas}/{contrasts}_chr{chr}.rvec", chr=CHR),
        FGr = "output/calculate_FGr/{dataset}/{gwas}/FGrMatPrime_{contrasts}.txt",
        snp_num = "output/calculate_FGr/{dataset}/{gwas}/SNPNum_{contrasts}.txt",
        Tvec = "data/{dataset}/TestVecs/{contrasts}.txt",
        pc_snps = "data/ukbb/ukbb_pc_snps_block.txt"
    output:
        H = "output/calculate_H/{dataset}/{gwas}/HPrime_{contrasts}.txt"
    shell:
        """
        Rscript code/calculate_H/calc_H.R {params.r_prefix} {input.FGr} {input.snp_num} \
        {input.Tvec} {output.H} {input.pc_snps}
        """
```

The estimator forms $\hat{H}' = \hat{\sigma}_{f'}^2 \cdot \hat{\sigma}_{r}^2$ and applies the same one-sided block-jackknife test against the detection limit $1/L$ described on the [uncorrected page](02-uncorrected-susceptibility):

```r
pvalNorm <- pnorm(H, mean = 1/L, sd = sqrt(varH), lower.tail = FALSE)
```

Because $\hat{H}'$ is typically driven down close to the detection limit, a non-significant result indicates only that residual structure is *below the noise floor* — not that it is exactly zero. Results across contrasts are concatenated into `plots/overlap_stats/HPrime_{dataset}.txt`.

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

To quantify how much of the target axis the PCs capture — without the result being inflated by estimation noise in $\hat{f}$ — we use an **even/odd cross-validation** scheme. PCs estimated on odd chromosomes are used to predict the target axis projected from even-chromosome SNPs. This is implemented in [`R2_EO.R`](https://github.com/jgblanc/quantifying-susceptibility-PGS.github.io/blob/master/scripts/calculate_FGr/R2_EO.R) (common or rare PCs) and [`R2_EO_both.R`](https://github.com/jgblanc/quantifying-susceptibility-PGS.github.io/blob/master/scripts/calculate_FGr/R2_EO_both.R) (common + rare PCs together):

```
rule calc_R2_even:
    input:
        PCs="output/calculate_PCA/{gwas}/odd_PCA.eigenvec",
        FGr="output/calculate_FGr/{dataset}/{gwas}/FGrMat_{contrasts}.txt",
        SNPs="output/calculate_FGr/{dataset}/{gwas}/SNPNum_{contrasts}.txt",
        r=expand("data/{dataset}/r/{gwas}/{contrasts}_chr{chr}.rvec", chr=CHR),
        pc_snps = "data/ukbb/ukbb_pc_snps_block.txt"
    output:
        R2="output/calculate_FGr/{dataset}/{gwas}/R2_EO/{contrasts}.txt"
    shell:
        """
        Rscript code/calculate_FGr/R2_EO.R {input.PCs} {input.FGr} {input.SNPs} \
        {output.R2} {params.r_prefix} {input.pc_snps}
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

The same step is run with rare PCs (`calc_R2_even_rare`, also via `R2_EO.R`) and with common + rare PCs jointly (`calc_R2_even_common_rare`, via `R2_EO_both.R`), letting us compare how much additional structure rare variants capture beyond common ones.

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
| `code/calculate_FGr/calc_GR_prime.R` | `scripts/calculate_FGr/calc_GR_prime.R` |
| `code/calculate_H/calc_H.R` | `scripts/calculate_H/calc_H.R` |
| `code/calculate_FGr/R2_EO.R` | `scripts/calculate_FGr/R2_EO.R` |
| `code/calculate_FGr/R2_EO_both.R` | `scripts/calculate_FGr/R2_EO_both.R` |

> **Note.** $H$ and $H'$ use the **same** estimator, `code/calculate_H/calc_H.R` — the only difference is whether the input projection matrix is `FGrMat` (uncorrected) or `FGrMatPrime` (PC-residualized). So the file you copy for page 02 is the same one this page links to; you don't need a separate copy.

> **Reminder.** As elsewhere, inline math underscores (e.g. `$\sigma_{f'}^2$`) need the same kramdown fix you applied earlier; the display equations (`$$ … $$` on their own lines) render fine as-is.

