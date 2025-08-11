#!/bin/bash
# Script: 11a_align_exons_by_gene.sh
# Purpose: Concatenate all exon FASTAs per gene and run MAFFT alignment

set -o pipefail
set -x

# Load required module
source /apps/profiles/modules_asax.sh.dyn
module load mafft/7.505

# Gene list
GENES=("IGF1" "IGF2" "INS" "IGF1R" "IGF2R" "INSR")

# Input and output directories
INPUT_ROOT="/scratch/aubsxb005/project1_data/project1_working_dataset/variants/01_final/exon_fasta_5/per_sample_per_gene_exon_truthful_fastas_5"
ALIGN_DIR="/scratch/aubsxb005/project1_data/project1_working_dataset/alignments/final_aligned_exons"
mkdir -p "$ALIGN_DIR"

# Loop over each gene
for GENE in "${GENES[@]}"; do
    echo "[INFO] Processing $GENE"

    GENE_DIR="${INPUT_ROOT}/${GENE}"
    COMBINED_FASTA="${ALIGN_DIR}/${GENE}_exon_combined.fasta"
    OUTPUT_ALN="${ALIGN_DIR}/${GENE}_exon_aligned.fasta"

    # Combine all FASTA files for this gene
    cat "${GENE_DIR}"/*.fasta > "$COMBINED_FASTA"

    # Run MAFFT
    mafft --auto "$COMBINED_FASTA" > "$OUTPUT_ALN"
done

echo "[DONE] Exon alignments completed."
-