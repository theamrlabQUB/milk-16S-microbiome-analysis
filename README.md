# Milk 16S microbiome analysis

Code for the manuscript "Nutritional composition and microbial community
membership of bovine milk vary across lactation with *Lactobacillus helveticus*
enrichment".

Raw PacBio HiFi full-length 16S rRNA gene reads are available from the European
Nucleotide Archive (ENA) under project accession **PRJEB127152**.

The repository contains two sets of scripts:

- **`.sh` / `.slurm` files**: QIIME 2 workflow for processing the raw reads
  (run on an HPC cluster with SLURM).
- **`.R` files**: downstream statistical analysis, tables and figures (R).

## 1. QIIME 2 workflow (`.sh` / `.slurm`)

Run in numerical order. Set the paths at the top of each script
(`RAW`, `PROJECT_DIR`, classifier path) before running.

| Script | Step |
|---|---|
| `01_make_manifest.sh` | Create the QIIME 2 manifest from the raw FASTQ files |
| `02_import_demux.sh` | Import HiFi reads and summarise demultiplexing |
| `03_denoise_all.slurm` | DADA2 `denoise-ccs` (primer removal, 1,000–1,600 bp, chimera removal) |
| `04_taxonomy_all.slurm` | Naive Bayes classification with the Greengenes2 2024.09 full-length classifier |
| `05_filter_low_depth.sh` | Remove samples listed in `samples_to_remove.tsv` |
| `06_filter_mitochondria_chloroplast.sh` | Remove mitochondrial and chloroplast features |
| `07_export_for_R.sh` | Export the feature table, taxonomy and representative sequences for R |

Software: QIIME 2 amplicon distribution 2024.10.

## 2. Downstream analysis in R (`.R`)

| Script | Manuscript output |
|---|---|
| `00_tables1_2_and_Figure1_milk_traits.R` | Tables 1–2, Figure 1 |
| `01_bacterial_community_composition_Figure2_SuppFigures3_4.R` | Figure 2, Supplementary Figures 3–4 |
| `02_alpha_beta_diversity_Figure3.R` | Figure 3 (alpha diversity tested within each farm; PERMANOVA, ANOSIM, PERMDISP) |
| `03_taxa_trait_correlations_Figures4_5_SuppFigures6_7.R` | Figures 4–5, Supplementary Figures 6–7 (Spearman correlations) |
| `04_differential_abundance_TableS8.R` | Supplementary Table S8 (CLR linear mixed models, top 30 species) |
| `05_global_microbiome_nutrition_TableS10.R` | Supplementary Table S10, Supplementary Figure 5 |
| `06_decontam_and_L_helveticus_TableS17.R` | decontam contaminant screen; *L. helveticus* ASVs for BLASTn |
| `07_taxa_per_rank_SupplementaryFigure2.R` | Supplementary Figure 2 |
| `08_method_agreement_TableS16_SupplementaryFigure8.R` | Supplementary Table S16, Supplementary Figure 8 |

Scripts 00–03 read their inputs from the paths set at the top of each script.
Scripts 04–08 read their inputs from a `data/` folder in the working directory:

```text
data/
├── processed/
│   ├── feature-table_milk_42samples.tsv   ASV counts, 42 milk samples (first line "#OTU ID ...")
│   ├── taxonomy.tsv                       Greengenes2 taxonomy
│   ├── sequences.fasta                    representative ASV sequences
│   └── species_rel_no30.tsv               species-level relative abundance table
├── metadata/
│   ├── metadata_clean.csv                 SampleID, Sample_Name, Lactation, Farm     (script 04)
│   ├── metadata_cow_name.csv              SampleID, CowID, Farm, Lactation, traits   (script 05)
│   └── metadata_cow.csv                   SampleID, CowID, Farm, Lactation           (script 06)
└── validation/
    └── validation_lactoscan_soxhlet_bradford.tsv                                      (script 08)
```

Run a script from the working directory that contains `data/`, for example:

```bash
Rscript 07_taxa_per_rank_SupplementaryFigure2.R
```

Results are written to `results/<analysis>/`. Each script saves its
`sessionInfo()` with the outputs.

R 4.4.3 with tidyverse, phyloseq, qiime2R, vegan, lme4, lmerTest, car,
decontam, patchwork and readxl.
