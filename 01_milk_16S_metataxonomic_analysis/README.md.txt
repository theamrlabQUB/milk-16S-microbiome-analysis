# Milk 16S Metataxonomic Analysis

This directory contains the scripts used for processing PacBio HiFi
full-length 16S rRNA gene amplicon sequencing data generated for this study.

The analysis was performed using QIIME 2 on a High Performance Computing
(HPC) cluster using the SLURM job scheduler.

## Analysis workflow

The workflow consists of the following steps:

1. Creation of a QIIME 2 manifest file from raw FASTQ files.
2. Import of PacBio HiFi reads into QIIME 2.
3. Generation of a demultiplexing summary.
4. Denoising and ASV generation using DADA2 `denoise-ccs`.
5. Taxonomic classification using the Greengenes2 reference database.
6. Filtering of low-depth samples.
7. Removal of mitochondrial and chloroplast-associated features.
8. Export of the feature table, taxonomy, representative sequences,
   and metadata for downstream analysis in R.

## Scripts

The scripts are located in the `scripts/` directory and are organised
according to the order of the analysis:

- `01_make_manifest.sh` – creates the QIIME 2 manifest file.
- `02_import_demux.sh` – imports sequencing reads into QIIME 2 and generates a demultiplexing summary.
- `03_denoise_all.slurm` – performs DADA2 denoising and ASV generation.
- `04_taxonomy_all.slurm` – assigns taxonomy using Greengenes2.
- `05_filter_low_depth.sh` – removes low-depth samples.
- `06_filter_mitochondria_chloroplast.sh` – removes mitochondrial and chloroplast-associated features.
- `07_export_for_R.sh` – exports processed data for downstream statistical analysis in R.

## Computing environment

The analysis was performed using:

- QIIME 2 amplicon distribution: `qiime2-amplicon-2024.10`
- Sequencing platform: PacBio HiFi
- Target region: full-length 16S rRNA gene
- Computing environment: High Performance Computing (HPC) cluster
- Job scheduler: SLURM

The QIIME 2 environment was activated using:

```bash
module load apps/anaconda3/2024.10
source $(conda info --base)/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10