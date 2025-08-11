#!/bin/bash
# Script: 12a_build_missing_cds_matrix_FINAL.sh
# Purpose: Generate a per-gene per-sample missing data matrix for CDS.
# Output: Tab-separated matrix (1 = present, 0 = absent/low-quality)

set -o pipefail
set -x

# === Load ASC module environment ===
source /apps/profiles/modules_asax.sh.dyn
module load gcc/8.2.0

# === Configuration ===
# Input directory containing per-gene folders, each with sample FASTAs
INPUT_BASE="/scratch/aubsxb005/project1_data/project1_working_dataset/variants/01_final/cds_fasta_5/per_sample_per_gene_CDS_truthful_fastas_5"

# Output matrix path
OUTPUT="/scratch/aubsxb005/project1_data/project1_working_dataset/missing_CDS_data_matrix.tsv"

# Quality thresholds
MIN_REAL_BASES=5        # Minimum number of A/T/G/C bases required
MIN_BASE_RATIO=0.05     # Minimum ratio of A/T/G/C bases to total length

echo "[INFO] Building CDS missing data matrix..."

# === Step 1: Collect unique sample names ===
samples=()
for gene_dir in "$INPUT_BASE"/*; do
    for fasta in "$gene_dir"/*.fasta; do
        [[ -f "$fasta" ]] || continue
        sample=$(basename "$fasta" | sed -E 's/^([A-Za-z0-9]+_[A-Za-z0-9]+)_.*\.fasta$/\1/')
        samples+=("$sample")
    done
done

# Deduplicate and sort
samples=($(printf "%s\n" "${samples[@]}" | sort -u))

# === Step 2: Write header to output file ===
{
    echo -ne "Gene"
    for s in "${samples[@]}"; do echo -ne "\t$s"; done
    echo
} > "$OUTPUT"

# === Step 3: Process each gene directory ===
for gene_dir in "$INPUT_BASE"/*; do
    gene=$(basename "$gene_dir")  # e.g., IGF1
    row="$gene"

    for sample in "${samples[@]}"; do
        fasta_file=$(find "$gene_dir" -type f -name "${sample}_*.fasta" | head -n 1)

        if [[ -f "$fasta_file" ]]; then
            seq=$(grep -v "^>" "$fasta_file" | tr -d '\n' | tr 'a-z' 'A-Z')
            seq_len=${#seq}
            real_bases=$(echo "$seq" | grep -o '[ATGC]' | wc -l)
            base_ratio=$(awk -v r="$real_bases" -v l="$seq_len" 'BEGIN {if (l > 0) print r / l; else print 0}')

            if [[ "$real_bases" -lt $MIN_REAL_BASES ]] || (( $(echo "$base_ratio < $MIN_BASE_RATIO" | bc -l) )); then
                row+="\t0"
            else
                row+="\t1"
            fi
        else
            row+="\t0"
        fi
    done

    echo -e "$row" >> "$OUTPUT"
done

echo "[DONE] Missing CDS matrix written to: $OUTPUT"
