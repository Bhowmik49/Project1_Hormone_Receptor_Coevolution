!/bin/bash

#------------------------------------------------------------------
# Script: 10b_extract_concat_cds_per_sample_and_gene.sh
# Purpose:
#   1. Extract CDS regions from consensus FASTA files
#      - Per-sample, per-gene concatenated FASTA
#   2. Create multi-FASTA files per gene across all samples
#------------------------------------------------------------------

set -o pipefail
set -x

######################################
# Load required modules
######################################
source /apps/profiles/modules_asax.sh.dyn
module load bedtools/2.26.0
module load gcc/13.2.0

######################################
# Set directories and feature
######################################
FEATURE="CDS"
ROOT="/scratch/aubsxb005/project1_data/project1_working_dataset"

BED_DIR="$ROOT/bed_files_per_gene_by_feature"
CONS_DIR="$ROOT/variants/consensus_fastas_5"

# Output: individual sample-gene files (concatenated)
OUT_SAMPLE_DIR="$ROOT/variants/cds_fasta_5/per_sample_per_gene_${FEATURE}_truthful_fastas_5"

# Output: multi-FASTA per gene across all samples
OUT_GENE_CONCAT="$ROOT/variants/cds_fasta_5/per_gene_${FEATURE}_truthful_fastas_5_${FEATURE}"

mkdir -p "$OUT_SAMPLE_DIR" "$OUT_GENE_CONCAT"


######################################
# Process each gene-specific BED file
######################################
for BED in "$BED_DIR"/*_"$FEATURE".bed.sorted; do
    GENE=$(basename "$BED" _${FEATURE}.bed.sorted)
    echo "[INFO] Processing gene: $GENE"

    # Create gene subfolder for per-sample files
    GENE_DIR="${OUT_SAMPLE_DIR}/${GENE}"
    mkdir -p "$GENE_DIR"

    # File for per-gene multi-FASTA
    GENE_FASTA="${OUT_GENE_CONCAT}/${GENE}_${FEATURE}.fasta"
    echo "" > "$GENE_FASTA"

    ######################################
    # Process each sample consensus FASTA
    ######################################
    for CONS in "$CONS_DIR"/*.fasta; do
        SAMPLE=$(basename "$CONS" _truthful_consensus.fasta)

        # Output: per-sample per-gene FASTA (concatenated)
        SAMPLE_FASTA="${GENE_DIR}/${SAMPLE}_${GENE}_${FEATURE}.fasta"

        # Extract CDS regions using BED
        bedtools getfasta -fi "$CONS" -bed "$BED" -name -fo "temp_${SAMPLE}_${GENE}.fa"

        # Concatenate all CDS fragments into one line
        CONCAT_SEQ=$(grep -v "^>" "temp_${SAMPLE}_${GENE}.fa" | tr -d '\n')

        # Write per-sample output
        echo ">${SAMPLE}" > "$SAMPLE_FASTA"
        echo "$CONCAT_SEQ" >> "$SAMPLE_FASTA"

        # Append to gene-wide multi-FASTA
        echo ">${SAMPLE}" >> "$GENE_FASTA"
        echo "$CONCAT_SEQ" >> "$GENE_FASTA"

        # Clean up
        rm -f "temp_${SAMPLE}_${GENE}.fa"
    done
done

echo "[DONE] All per-sample and per-gene CDS FASTAs created successfully."

