#!/bin/bash
#SBATCH --job-name=tax_all
#SBATCH --partition=k2-medpri
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=logs/tax_all_%j.out
#SBATCH --error=logs/tax_all_%j.err

set -euo pipefail

module load apps/anaconda3/2024.10
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate qiime2-amplicon-2024.10

# Set PROJECT_DIR to the QIIME 2 working directory
PROJECT_DIR=path/to/16s_milk/qiime2
cd "$PROJECT_DIR"

mkdir -p logs

qiime feature-classifier classify-sklearn \
    --i-classifier path/to/databases/gg2_2024.09_full_length_classifier.nb.qza \
    --i-reads repseq_all.qza \
    --p-n-jobs 8 \
    --o-classification taxonomy_gg2_full.qza