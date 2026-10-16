# Milk 16S microbiome analysis

Code for the manuscript "Nutritional composition and microbial community
membership of bovine milk vary across lactation with *Lactobacillus helveticus*
enrichment".

Raw PacBio HiFi full-length 16S rRNA gene reads are available from the European
Nucleotide Archive (ENA) under project accession **PRJEB127152**.

## Repository structure

```text
milk-16S-microbiome-analysis/
├── 01_milk_16S_metataxonomic_analysis/   QIIME 2 workflow (HPC, SLURM)
│   ├── README.md
│   └── scripts/
│       ├── 01_make_manifest.sh
│       ├── 02_import_demux.sh
│       ├── 03_denoise_all.slurm
│       ├── 04_taxonomy_all.slurm
│       ├── 05_filter_low_depth.sh
│       ├── 06_filter_mitochondria_chloroplast.sh
│       └── 07_export_for_R.sh
└── 02_downstream_R_analysis/             statistics, tables and figures (R)
    ├── README.md
    ├── 00_tables1_2_and_Figure1_milk_traits.R
    ├── 01_bacterial_community_composition_Figure2_SuppFigures3_4.R
    ├── 02_alpha_beta_diversity_Figure3.R
    ├── 03_taxa_trait_correlations_Figures4_5_SuppFigures6_7.R
    ├── 04_differential_abundance_TableS8.R
    ├── 05_global_microbiome_nutrition_TableS10.R
    ├── 06_decontam_and_L_helveticus_TableS17.R
    ├── 07_taxa_per_rank_SupplementaryFigure2.R
    └── 08_method_agreement_TableS16_SupplementaryFigure8.R
```

## Workflow

1. Run the QIIME 2 scripts in `01_milk_16S_metataxonomic_analysis/scripts/` in
   numerical order to denoise the reads (DADA2), assign taxonomy (Greengenes2),
   remove mitochondrial and chloroplast features and export the feature table,
   taxonomy and representative sequences.
2. Place the exported files as described in `02_downstream_R_analysis/README.md`
   and run the R scripts to reproduce the tables and figures.

## Software

QIIME 2 2024.10; R 4.4.3 with tidyverse, phyloseq, qiime2R, vegan, lme4,
lmerTest, decontam and patchwork. Package versions are recorded in the
`sessionInfo()` output saved by each R script.
