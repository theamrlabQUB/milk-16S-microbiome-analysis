#!/bin/bash
set -euo pipefail

mkdir -p filtered

qiime taxa filter-table \
    --i-table filtered/table_filtered.qza \
    --i-taxonomy taxonomy_gg2_full.qza \
    --p-exclude mitochondria,chloroplast \
    --o-filtered-table filtered/table_noMitoChloro.qza

qiime feature-table filter-seqs \
    --i-data filtered/repseq_filtered.qza \
    --i-table filtered/table_noMitoChloro.qza \
    --o-filtered-data filtered/repseq_noMitoChloro.qza