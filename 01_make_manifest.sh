#!/bin/bash
set -euo pipefail

# Set RAW to the folder containing the demultiplexed HiFi FASTQ files
RAW=path/to/16s_milk/rawdata
OUT=manifests/manifest_all.tsv

mkdir -p manifests

echo -e "sample-id\tabsolute-filepath" > "$OUT"

for f in "$RAW"/*.fastq.gz; do
    base=$(basename "$f")
    sid=${base%%.fastq.gz}
    echo -e "${sid}\t${f}" >> "$OUT"
done

head "$OUT"
wc -l "$OUT"