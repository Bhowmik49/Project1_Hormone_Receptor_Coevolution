#!/bin/bash
# Script: 11b_align_cds_by_gene_final.sh
# Purpose: Align CDS multi-FASTA per gene using MAFFT

source /apps/profiles/modules_asax.sh.dyn
module load mafft/7.505

set -o pipefail
set -x

# Genes to process
GENES=("IGF1" "IGF2" "INS" "IGF1R" "IGF2R" "INSR")

# Correct input and output directories
INPUT_ROOT="/scratch/aubsxb005/project1_data/project1_working_dataset/variants/01_final/cds_fasta_5/per_sample_per_gene_CDS_truthful_fastas_5"
ALIGN_DIR="/scratch/aubsxb005/project1_data/project1_working_dataset/alignments/final_aligned_CDS"
mkdir -p "$ALIGN_DIR"

# Main loop
for GENE in "${GENES[@]}"; do
    echo "[INFO] Aligning CDS for $GENE"

    GENE_DIR="${INPUT_ROOT}/${GENE}"
    TEMP_CAT="${ALIGN_DIR}/${GENE}_CDS_combined.fasta"
    OUTPUT_ALN="${ALIGN_DIR}/${GENE}_CDS_aligned.fasta"

    # Check and concatenate FASTAs
    if ls "${GENE_DIR}"/*.fasta 1> /dev/null 2>&1; then
        cat "${GENE_DIR}"/*.fasta > "$TEMP_CAT"
        mafft --auto "$TEMP_CAT" > "$OUTPUT_ALN"
    else
        echo "[WARN] No CDS FASTA files found for $GENE — skipping."
    fi
done

echo "[DONE] CDS alignments completed in $ALIGN_DIR"
