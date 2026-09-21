#!/bin/bash
set -euo pipefail

qiime tools import \
    --type 'SampleData[SequencesWithQuality]' \
    --input-path manifests/manifest_all.tsv \
    --output-path demux_all.qza \
    --input-format SingleEndFastqManifestPhred33V2

qiime demux summarize \
    --i-data demux_all.qza \
    --o-visualization demux_all.qzv