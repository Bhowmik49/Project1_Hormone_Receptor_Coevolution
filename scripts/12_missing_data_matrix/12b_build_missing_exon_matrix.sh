#!/bin/bash
# Script: 12b_build_missing_exon_matrix_FINAL.sh
# Purpose: Generate a per-gene per-sample missing data matrix for exon sequences.
# Output: Tab-separated matrix (1 = present, 0 = absent/low-quality)

set -o pipefail
set -x

# === Load environment ===
source /apps/profiles/modules_asax.sh.dyn
module load gcc/8.2.0

# === Configuration ===
INPUT_BASE="/scratch/aubsxb005/project1_data/project1_working_dataset/variants/01_final/exon_fasta_5/per_sample_per_gene_exon_truthful_fastas_5"
OUTPUT="/scratch/aubsxb005/project1_data/project1_working_dataset/missing_exon_data_matrix.tsv"

# Thresholds
MIN_REAL_BASES=5
MIN_BASE_RATIO=0.05

echo "[INFO] Building exon missing data matrix..."

# === Step 1: Extract unique sample IDs ===
samples=()
for gene_dir in "$INPUT_BASE"/*; do
    for fasta in "$gene_dir"/*.fasta; do
        [[ -f "$fasta" ]] || continue
        sample=$(basename "$fasta" | sed -E 's/^([A-Za-z0-9]+_[A-Za-z0-9]+)_.*$/\1/')
        samples+=("$sample")
    done
done

# Deduplicate and sort
samples=($(printf "%s\n" "${samples[@]}" | sort -u))

# === Step 2: Header ===
{
    echo -ne "Gene"
    for s in "${samples[@]}"; do echo -ne "\t$s"; done
    echo
} > "$OUTPUT"

# === Step 3: For each gene, create matrix row ===
for gene_dir in "$INPUT_BASE"/*; do
    gene=$(basename "$gene_dir")  # e.g., IGF1, IGF2R, etc.
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

echo "[DONE] Exon matrix written to: $OUTPUT"
