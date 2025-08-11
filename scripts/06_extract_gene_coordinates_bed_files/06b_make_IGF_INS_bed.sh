#!/bin/bash
# Description: Convert gene coordinates to properly formatted BED file
# Author: aubsxb005

set -euo pipefail

INPUT="IGF_INS_gene_coords.txt"
OUTPUT="IGF_INS_targets.bed"

# Convert to BED format: 0-based start, include gene_feature as name
echo "[INFO] Creating BED file: $OUTPUT"
tail -n +2 "$INPUT" | awk -F'\t' 'BEGIN{OFS="\t"} {
    chrom = $1
    start = $2 - 1  # BED format uses 0-based start
    end = $3        # BED uses 1-based end
    feature = $4
    gene = $5
    strand = $6
    name = gene "_" feature
    print chrom, start, end, name, ".", strand
}' | sort -k1,1 -k2,2n > "$OUTPUT"

echo "[INFO] BED file successfully created at: $OUTPUT"