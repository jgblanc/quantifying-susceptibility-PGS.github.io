---
layout: default
title: Datasets
---

[← Back to Home]({{ '/' | relative_url }})

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

- **A) Pairwise continental contrasts** — all $\binom{5}{2} = 10$ pairwise combinations of the five ancestry groups above
- **B) Continuous geographic gradients in Eurasia** — latitude and longitude within Eurasia
- **C) Continuous geographic gradients in Europe** — latitude and longitude within non-Finnish Europe
- **D) Sardinia vs. mainland Europe** — a contrast of particular interest given prior work on height PGS divergence

```r
# Read in metadata
df <- fread("../data/HGDP1KG/gnomad_meta_updated.tsv", sep = "\t")

# Select relevant columns
dfFilter <- df %>% dplyr::select("project_meta.project_subpop", "latitude", "longitude",
                                  "population", "region", "project_meta.project_pop")
colnames(dfFilter)[1] <- "IID"
```

![Map of HGDP1kGP samples colored by population](../assets/images/HGDPMaps.png)

---

### A. All Continental (N = 3,859)

We included all five major population groups to test all $\binom{5}{2} = 10$ pairwise continental contrasts.

| Population | Code | Sample Size |
|------------|------|-------------|
| East Asia | eas | 825 |
| Africa | afr | 1,003 |
| South Asia | sas | 790 |
| Non-Finnish Europe | nfe | 689 |
| America | amr | 552 |

```r
# Select all five continental populations
dfAll <- dfFilter %>%
  filter(project_meta.project_pop %in% c("nfe", "afr", "amr", "eas", "sas")) %>%
  select("IID", "project_meta.project_pop", "population", "latitude", "longitude")
colnames(dfAll) <- c("IID", "pop", "subpop", "lat", "long")
dfAll <- dfAll %>% mutate(FID = 0) %>% select(FID, everything())
```
--- 

### B. Eurasia (N = 2,224)

We selected individuals from four broad Eurasian ancestry groups to test for polygenic selection along latitude and longitude, following [Berg et al.](https://elifesciences.org/articles/39725). CEU samples were excluded because their recorded coordinates correspond to the sample collection site in Utah rather than their ancestral geography.

| Population | Code | Sample Size |
|------------|------|-------------|
| East Asia | eas | 825 |
| Finland | fin | 99 |
| Non-Finnish Europe | nfe | 510 |
| South Asia | sas | 790 |

```r
# Select Eurasian samples, removing CEU
dfEurAsia <- dfFilter %>%
  filter(project_meta.project_pop %in% c("fin", "nfe", "eas", "sas")) %>%
  filter(population != "CEU")
colnames(dfEurAsia) <- c("IID", "pop", "subpop", "lat", "long")
dfEurAsia <- dfEurAsia %>% mutate(FID = 0) %>% select(FID, everything())

fwrite(dfEurAsia, "../data/HGDP1KG/ids/test_ids/eurasia.txt", sep = "\t")
```

---

### C. Non-Finnish Europe (N = 510)

We selected only non-Finnish European individuals to test for polygenic selection along latitude and longitude within Europe, again following [Berg et al.](https://elifesciences.org/articles/39725). CEU samples were again excluded.

| Population | Code | Sample Size |
|------------|------|-------------|
| Non-Finnish Europe | nfe | 510 |

```r
# Select non-Finnish European samples, removing CEU
dfEur <- dfFilter %>%
  filter(project_meta.project_pop %in% c("nfe")) %>%
  filter(population != "CEU")
dfEur <- dfEur %>% dplyr::select("IID", "project_meta.project_pop", "population", "latitude", "longitude")
colnames(dfEur) <- c("IID", "pop", "subpop", "lat", "long")
dfEur <- dfEur %>% mutate(FID = 0) %>% select(FID, everything())
```

---

### D. Sardinia vs. Mainland Europe (N = 689)

Following [Chen et al.](https://www.sciencedirect.com/science/article/pii/S0002929720301610), we tested for polygenic selection acting on height in Sardinia compared to the rest of Europe. Here CEU samples are included, and individuals are coded as Sardinian (`sdi`) or mainland European (`eur`).

| Group | Code | Sample Size |
|-------|------|-------------|
| Sardinian | sdi | 27 |
| Mainland Europe | eur | 662 |

```r
# Select NFE samples (including CEU) and code Sardinians vs mainland
dfEur <- dfFilter %>% filter(project_meta.project_pop %in% c("nfe"))
dfEur <- dfEur %>% select("IID", "project_meta.project_pop", "population", "latitude", "longitude")
colnames(dfEur) <- c("IID", "pop", "subpop", "lat", "long")
dfEur <- dfEur %>%
  mutate(SDI = case_when(subpop == "Sardinian" ~ "sdi", TRUE ~ "eur")) %>%
  dplyr::select(IID, SDI, everything())
dfEur <- dfEur %>% mutate(FID = 0) %>% select(FID, everything())
```

---

## UK Biobank (UKB)

### British Isles Country of Birth (CoB) Prediction Panel

To investigate fine-scale axes of variation, we selected UK Biobank individuals whose country of birth was listed as one of five British Isles countries and whose self-identified ethnic background was recorded as White. This gave $\binom{5}{2} = 10$ pairwise country-of-birth contrasts.

```r
# Load country of birth data and assign labels
dfCoBUK <- fread("../data/InUKBB/data/CountryOfBirthUK_1647.txt")
dfCoBUK <- dfCoBUK %>% filter(CountryOfBirthUK_1647 > 0) %>%
  mutate(CoB = case_when(CountryOfBirthUK_1647 == 1 ~ "England",
                         CountryOfBirthUK_1647 == 2 ~ "Wales",
                         CountryOfBirthUK_1647 == 3 ~ "Scotland",
                         CountryOfBirthUK_1647 == 4 ~ "NorthernIreland",
                         CountryOfBirthUK_1647 == 5 ~ "RepublicOfIreland",
                         CountryOfBirthUK_1647 == 6 ~ "Elsewhere"))

# Join with PCA data and restrict to self-identified European ancestry
dfCoBUK <- inner_join(dfCoBUK, df) %>% filter(continental == "Europe")
```

Individuals coded as "Elsewhere" were excluded from the prediction panel. The final sample sizes for each country of birth are:

| Country | Code | Sample Size |
|---------|------|-------------|
| England | ENG | 142,215 |
| Northern Ireland | NI | 1,120 |
| Republic of Ireland | RoI | 1,802 |
| Scotland | SCT | 14,591 |
| Wales | WAL | 8,125 |

---

## UK Biobank GWAS Panels

We constructed six GWAS panels from the UK Biobank representing a gradient from homogeneous to maximally diverse. The intermediate panels ($5\varepsilon$ through $200\varepsilon$) were constructed using a continuous sampling scheme based on Euclidean distance from the centroid of the full biobank in PC space, where $\varepsilon = 10^{-4}$.

| Panel | Description | Sample Size |
|-------|-------------|-------------|
| WBS | White British ancestry subset | 408,626 |
| WBS-R | White British ancestry subset, CoB Removed | 256,582 |
| $5\varepsilon$ | Random sample within $5\varepsilon$ of PC centroid | 100,000 |
| $10\varepsilon$ | Random sample within $10\varepsilon$ of PC centroid | 100,000 |
| $50\varepsilon$ | Random sample within $50\varepsilon$ of PC centroid | 100,000 |
| $200\varepsilon$ | Random sample within $200\varepsilon$ of PC centroid | 100,000 |
| ALL-R | Full UK Biobank | 318,560 |
| ALL | Full UK Biobank | 486,413 |

### Visualizing the full UK Biobank

We first performed PCA on the full UK Biobank and colored individuals by self-identified ethnic background to visualize the overall diversity of the cohort.

```r
# Load PCA results and ethnic background labels
dfPCA <- fread("../plots/ukbb/whole_biobank.eigenvec.gz")
dfEB <- fread("../plots/ukbb/EthnicBackground_21000.txt")
df <- inner_join(dfPCA, dfEB)

# Assign continental labels based on UKB ethnic background codes
df <- df %>% mutate(continental = case_when(
  EthnicBackground_21000 %in% c(4, 4001, 4002, 4003) ~ "Africa",
  EthnicBackground_21000 %in% c(1, 1001, 1002, 1003) ~ "Europe",
  EthnicBackground_21000 %in% c(2, 2001, 2002, 2003, 2004) ~ "Mixed",
  EthnicBackground_21000 %in% c(3, 3001, 3002, 3003, 3004) ~ "Asia"))

# Plot PC1 vs PC2 colored by self-identified ethnic background
ggplot(data = df, aes(x = PC1, y = PC2, color = continental)) +
  geom_point() +
  theme_classic(base_size = 14) +
  scale_color_manual(values = c("darkred", "navy", "goldenrod4", "gray50"),
                     na.value = "gray70")
```

![PC1 vs PC2 of the full UK Biobank colored by self-identified ethnic background](../assets/images/PCA.png)

### Constructing the intermediate GWAS panels

We calculated each individual's Euclidean distance from the PC centroid of the full biobank and sampled 100,000 individuals falling within each distance threshold.

```r
# Calculate distance from the PC centroid using PC1 and PC2
medianPC2 <- apply(df[,3:4], 2, median)
df$distance2 <- sqrt(rowSums((df[,3:4] - medianPC2)^2))

# Assign panel membership based on distance thresholds
epsilon <- 1e-4
dfDS <- df %>% mutate(
  e5   = distance2 <= (5   * epsilon),
  e10  = distance2 <= (10  * epsilon),
  e50  = distance2 <= (50  * epsilon),
  e200 = distance2 <= (200 * epsilon))

# Sample 100,000 individuals from each panel and save
set.seed(12121212)
dfe5   <- dfDS %>% filter(e5   == TRUE) %>% sample_n(100000) %>% select("#FID", "IID", "POP")
dfe10  <- dfDS %>% filter(e10  == TRUE) %>% sample_n(100000) %>% select("#FID", "IID", "POP")
dfe50  <- dfDS %>% filter(e50  == TRUE) %>% sample_n(100000) %>% select("#FID", "IID", "POP")
dfe200 <- dfDS %>% filter(e200 == TRUE) %>% sample_n(100000) %>% select("#FID", "IID", "POP")
```

Each sampled panel is visualized below, with selected individuals highlighted against the full biobank:

```r
# Plot each panel against the full biobank
p1 <- ggplot(data = dfPlot, aes(x = PC1, y = PC2, color = Sample)) +
  geom_point() + theme_classic(base_size = 14) +
  scale_color_manual(values = c("grey70", "purple1")) +
  ggtitle(TeX("$5 \\epsilon$")) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

# (repeated for e10, e50, e200 with colors darkorange, skyblue2, seagreen)

p <- grid.arrange(p1, p2, p3, p4, nrow = 1)
```

![The four intermediate GWAS panels highlighted against the full UK Biobank in PC space](../assets/images/gwasEpanels.png)

### Removing country of birth overlap

Individuals used to construct the country of birth prediction panel were removed from the WBS and ALL GWAS panels to avoid overlap between the prediction and GWAS samples.

```r
# Remove withdrawn participants
dfWBS <- fread("../plots/ukbb/WBS.txt") %>% select("#FID", "IID", "POP")
dfALL <- fread("../plots/ukbb/ALL.txt") %>% select("#FID", "IID", "POP")

# Remove CoB individuals from WBS and ALL
dfCoB <- fread("../plots/ukbb/CountryOfBirthUK.txt") %>% filter(!IID %in% wd$IID)
dfWBS_removed <- dfWBS %>% filter(!IID %in% dfCoB$IID)
dfALL_removed <- dfALL %>% filter(!IID %in% dfCoB$IID)

fwrite(dfWBS_removed, "../plots/ukbb/WBS-R.txt", sep = "\t")
fwrite(dfALL_removed, "../plots/ukbb/ALL-R.txt", sep = "\t")
```

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
