#!/bin/bash
set -euo pipefail

mkdir -p filtered

qiime feature-table filter-samples \
    --i-table table_95_no30.qza \
    --m-metadata-file filtered/samples_to_remove.tsv \
    --p-exclude-ids \
    --o-filtered-table filtered/table_filtered.qza

qiime feature-table filter-seqs \
    --i-data repseq_95_no30.qza \
    --i-table filtered/table_filtered.qza \
    --o-filtered-data filtered/repseq_filtered.qza