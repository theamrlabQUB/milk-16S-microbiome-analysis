# Downstream R Analysis

This directory contains the R scripts used for downstream statistical analysis, visualisation, and generation of manuscript tables and figures for the milk microbiome study.

Raw sequencing data are not stored in this GitHub repository. Sequencing data generated for this study are deposited in the European Nucleotide Archive (ENA).

ENA project/study accession: attached in the manuscript

Users should update the input file paths at the beginning of each script before running the analyses.

## Repository structure

```text
02_downstream_R_analysis/
├── README.md
└── scripts/
    ├── 00_tables1_2_and_Figure1_milk_traits.R
    ├── 01_bacterial_community_composition_Figure2_SuppFigures2_3_4.R
    ├── 02_alpha_beta_diversity_Figure3.R
    └── 03_taxa_trait_correlations_Figures4_5_SuppFigures5_6.R
```

## Script 00: Tables 1–2 and Figure 1

`00_tables1_2_and_Figure1_milk_traits.R`

This script analyses milk nutritional and production traits and generates Tables 1 and 2 and Figure 1.

### Table 1

Summarises milk composition according to lactation stage, farm, and parity group.

The analysed traits include:

- Fat
- Protein
- Lactose
- Milk energy
- Milk production
- Minerals

Results are presented as mean ± SEM.

Statistical analyses include:

- Paired t-tests for Early- vs Late-lactation comparisons using Cow ID
- Welch independent t-tests for parity effects
- Welch independent t-tests for farm effects

### Table 2

Compares milk nutritional and production traits between samples included and excluded from the microbiome analysis.

For each group, Early- and Late-lactation values are summarised as mean ± SEM and compared using paired t-tests based on Cow ID.

### Figure 1

**Principal component analysis (PCA) biplot showing variation in milk traits between the two farms.**

The PCA includes:

- Fat
- Protein
- Lactose
- Milk energy
- Milk production
- Minerals

Before PCA, multivariate outliers are identified using Mahalanobis distance calculated from standardised trait values. Observations exceeding the chi-square threshold at P < 0.001 are excluded.

PCA is then performed using centred and scaled variables.

The figure displays individual samples, farm-specific symbols, 95% confidence ellipses, and loading vectors for the six milk traits.

## Script 01: Figure 2 and Supplementary Figures 2–4

`01_bacterial_community_composition_Figure2_SuppFigures2_3_4.R`

This script analyses bacterial community composition using the QIIME 2 feature table, taxonomy assignments, and sample metadata.

### Figure 2

**Bacterial community composition across lactation stages and farms.**

The figure includes:

- Mean relative abundance of the top 20 species in Early- and Late-lactation milk samples
- Species-level composition stratified by farm and lactation stage
- Phylum-level composition stratified by farm and lactation stage
- Genus-level composition stratified by farm and lactation stage

Taxa outside the top 20 are grouped as `Other`.

### Supplementary Figure 2

**Identified taxa from raw bovine milk samples.**

Features are classified according to the deepest resolved taxonomic level:

- Unclassified
- Phylum
- Class
- Order
- Family
- Genus
- Species unspecified
- Species identified

### Supplementary Figure 3

**Relative abundance of the top 20 dominant microbial taxa across cow samples during Early and Late lactation in Farms A and B.**

Stacked bar plots show community composition at the phylum, genus, and species levels.

Each bar represents an individual cow sample. Early- and Late-lactation samples are paired within each farm.

### Supplementary Figure 4

**Distribution of the top 20 rare bacterial taxa across milk samples.**

Rare taxa are defined using a relative-abundance threshold of ≤ 0.01%.

The figure shows rare species:

- By lactation stage
- By farm
- Across individual cow milk samples

## Script 02: Figure 3

`02_alpha_beta_diversity_Figure3.R`

This script performs alpha- and beta-diversity analyses of the milk microbiome.

### Alpha diversity

The following diversity indices are calculated:

- Observed richness
- Chao1 richness
- Shannon diversity
- Simpson diversity

Early- and Late-lactation samples from the same cows are compared using paired Wilcoxon signed-rank tests.

P values are adjusted using the Benjamini-Hochberg false-discovery-rate procedure.

### Beta diversity

The analysis includes:

- Bray-Curtis dissimilarity
- Binary Jaccard dissimilarity
- Principal coordinates analysis (PCoA)
- Farm-adjusted paired PERMANOVA
- ANOSIM
- PERMDISP

The script generates:

**Figure 3. Alpha and beta diversity of the milk microbiome according to lactation stage.**

## Script 03: Figures 4–5 and Supplementary Figures 5–6

`03_taxa_trait_correlations_Figures4_5_SuppFigures5_6.R`

This script performs lactation stage-specific Spearman rank correlations between bacterial taxa and milk nutritional and production traits.

The analysed traits include:

- Milk production
- Fat
- Energy
- Lactose
- Minerals
- Protein

Taxa are prevalence-filtered before selecting the top 50 species and genera.

Analyses are performed separately for:

- Farm A – Early lactation
- Farm A – Late lactation
- Farm B – Early lactation
- Farm B – Late lactation
- Combined farms – Early lactation
- Combined farms – Late lactation

### Figure 4

**Lactation stage-specific correlations between bacterial species and milk nutritional and production traits.**

Significance annotations are based on Benjamini-Hochberg FDR-adjusted P values.

### Figure 5

**Lactation stage-specific correlations between bacterial genera and milk nutritional and production traits.**

Significance annotations are based on Benjamini-Hochberg FDR-adjusted P values.

### Supplementary Figure 5

Species-level correlations before multiple-testing correction.

Significance annotations are based on raw P values.

### Supplementary Figure 6

Genus-level correlations before multiple-testing correction.

Significance annotations are based on raw P values.

For the correlation heatmaps:

- `*` indicates P or FDR < 0.05
- `**` indicates P or FDR < 0.01
- `***` indicates P or FDR < 0.001

Main manuscript figures use FDR-adjusted P values, whereas supplementary figures use unadjusted P values.

## Required R packages

The scripts use packages including:

- tidyverse
- ggplot2
- vegan
- patchwork
- stringr
- phyloseq
- qiime2R
- readxl
- openxlsx
- ggrepel
- scales
- grid

Individual scripts check their required packages before running.

## Data availability

Raw sequencing data and associated sample metadata generated in this study have been deposited in the European Nucleotide Archive (ENA) under project accession PRJEB127152.


Raw FASTQ files and large QIIME 2 artifacts are not stored in this GitHub repository.

Users wishing to reproduce the analyses should obtain the required data, update the input file paths in the scripts, and run the scripts as appropriate.

## Reproducibility

The scripts include:

- Documented input requirements
- Statistical analysis steps
- Output file generation
- R session information for package and software version tracking
- Fixed random seeds where permutation-based analyses are used

