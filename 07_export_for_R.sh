#!/bin/bash
set -euo pipefail

mkdir -p export_for_R/table
mkdir -p export_for_R/taxonomy
mkdir -p export_for_R/repseqs

qiime tools export \
    --input-path filtered/table_noMitoChloro.qza \
    --output-path export_for_R/table

qiime tools export \
    --input-path taxonomy_gg2_full.qza \
    --output-path export_for_R/taxonomy

qiime tools export \
    --input-path filtered/repseq_noMitoChloro.qza \
    --output-path export_for_R/repseqs

cp metadata_qiime_no30.tsv export_for_R/metadata.tsv

cd export_for_R/table

biom convert \
    -i feature-table.biom \
    -o feature-table.tsv \
    --to-tsv

# biom convert writes "# Constructed from biom file" as the first line;
# remove it so that "#OTU ID" is the header line read by the R scripts.
sed -i '1{/^# Constructed from biom file/d}' feature-table.tsv
