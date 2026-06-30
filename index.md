---
layout: default
title: Quantifying the Susceptibility of Polygenic Scores to Ancestry Stratification
---

# Quantifying the Susceptibility of Polygenic Scores to Ancestry Stratification

**Jennifer Blanc, Walid Mawass, Jeremy J. Berg**  
Department of Human Genetics, University of Chicago

[📄 Read the Preprint](https://www.biorxiv.org/content/10.64898/2025.12.04.692430v1) | [💻 View the Original Code](https://github.com/jgblanc/strat2)

This website was created with the aide of Claude Opus 4.8 (Anthropic).

---

## Overview

Polygenic scores (PGS) aggregate the effects of thousands of genetic variants to predict complex traits. However, they remain vulnerable to spurious correlations arising from environmental variation that covaries with population structure — a phenomenon known as **ancestry stratification**.

This site documents the computational workflow accompanying Blanc et al. (2025), in which we develop a framework to quantify how susceptible a polygenic score is to stratification bias along specific, researcher-defined ancestry gradients. We introduce two key statistics:

- **H** — the proportion of genetic variance in a GWAS panel explained by an external ancestry gradient (uncorrected susceptibility)
- **H′** — the residual susceptibility remaining after correction via Principal Component Analysis (PCA)

Applying this framework to the UK Biobank, we find that robust PCA correction in diverse cohorts reduces residual susceptibility to levels comparable to — or lower than — those achieved by restricting analysis to homogeneous subsets. We then perform polygenic score association tests for 17 phenotypes to demonstrate the practical utility of our diagnostic, recovering well-known and novel signal divergence. 

---

## Pipeline Overview

The analysis is implemented as a Snakemake workflow. The pages below document each major step:

| Step | Description |
|------|-------------|
| [Datasets](docs/07-datasets) | HGDP1kGP, UK Biobank panels, and prediction sets |
| [1. Genotype Contrasts](docs/01-genotype-contrasts) | Computing allele frequency contrasts ($\hat{r}$) |
| [2. Uncorrected Susceptibility](docs/02-uncorrected-susceptibility) | Estimating H and the target axis f |
| [3. PCA Correction](docs/03-pca-correction) | Common and rare variant principal components |
| [4. Residual Susceptibility](docs/04-residual-susceptibility) | Estimating H′ and PCA efficacy $V_K$ |
| [5. GWAS](docs/05-gwas) | Running GWAS with regenie |
| [6. PGS Association Tests](docs/06-pgs-association-tests) | Clumping, thresholding, and the $\hat{q}$ statistic |

---

The original Snakefiles used for the analysis on the UChicago Center for Research Informatics [High Performance Computing Cluster](https://cri.uchicago.edu/hpc/) are located in [snakefiles_raw](snakefiles_raw)  

## Citation

Blanc J, Mawass W, Berg JJ. *Quantifying the susceptibility of polygenic scores to ancestry stratification.* bioRxiv 2025. [https://doi.org/10.64898/2025.12.04.692430](https://doi.org/10.64898/2025.12.04.692430)


