---
layout: default
title: Datasets
---

# Datasets

This page describes the datasets used to construct prediction panels, genotype contrasts, and GWAS panels.

---

## HGDP and 1000 Genomes Project (HGDP1kGP)

We used high-coverage whole genome sequencing data from the harmonized Human Genome Diversity Project (HGDP) and 1000 Genomes Project (1kGP), comprising 4,094 individuals. Individuals were grouped according to project meta population codes. We selected groups with at least 200 samples, resulting in five broad ancestry groups:

| Group | Code | Sample Size |
|-------|------|-------------|
| East Asia | eas | 825 |
| Africa | afr | 1,003 |
| South Asia | sas | 790 |
| Non-Finnish Europe | nfe | 689 |
| America | amr | 552 |

We constructed three sets of genotype contrasts from this dataset:

- **Pairwise continental contrasts** — all $\binom{5}{2} = 10$ pairwise combinations of the five ancestry groups above
- **Continuous geographic gradients** — latitude and longitude within Eurasia, latitude and longitude within non-Finnish Europe, and Sardinia vs. mainland Europe
- **Sardinia vs. mainland Europe** — a contrast of particular interest given prior work on height PGS divergence

---

## UK Biobank (UKB)

### British Isles Country of Birth (CoB) Prediction Panel

To investigate fine-scale axes of variation, we selected UK Biobank individuals whose country of birth was listed as one of five British Isles countries and whose self-identified ethnic background was recorded as White. This gave $\binom{5}{2} = 10$ pairwise country-of-birth contrasts.

| Country | Code | Sample Size |
|---------|------|-------------|
| England | ENG | 142,215 |
| Northern Ireland | NI | 1,120 |
| Republic of Ireland | RoI | 1,802 |
| Scotland | SCT | 14,591 |
| Wales | WAL | 8,125 |

### GWAS Panels

We constructed six GWAS panels from the UK Biobank representing a gradient from homogeneous to maximally diverse:

| Panel | Description | Sample Size |
|-------|-------------|-------------|
| WBS | White British ancestry subset | ~408,625 |
| 5ε | Random sample within 5ε of PC centroid | 100,000 |
| 10ε | Random sample within 10ε of PC centroid | 100,000 |
| 50ε | Random sample within 50ε of PC centroid | 100,000 |
| 200ε | Random sample within 200ε of PC centroid | 100,000 |
| ALL | Full UK Biobank | ~486,612 |

The intermediate panels (5ε–200ε) were constructed using a continuous sampling scheme based on Euclidean distance from the centroid of the full biobank in PC space (where $\varepsilon = 10^{-4}$). This allows us to evaluate the trade-off between sample diversity and stratification control across a spectrum of study designs.

---

## Phenotypes

We conducted GWAS and polygenic score association tests for 17 phenotypes from the UK Biobank:

| Phenotype | UKB Code |
|-----------|----------|
| Standing Height | 50 |
| Body Weight | 21002 |
| Systolic Blood Pressure | 4080 |
| Cholesterol | 30690 |
| HDL | 30760 |
| LDL | 30780 |
| Triglycerides | 30870 |
| Glucose | 30740 |
| HbA1c | 30750 |
| Platelet Count | 30080 |
| Mean Corpuscular Volume (MCV) | 30040 |
| Mean Corpuscular Haemoglobin Concentration (MCHC) | 30060 |
| Alkaline Phosphatase | 30610 |
| Aspartate Aminotransferase | 30650 |
| Basophil Percentage | 30220 |
| Eosinophil Percentage | 30210 |
| Total Protein | 30860 |
