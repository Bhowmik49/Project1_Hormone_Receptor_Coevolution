#!/bin/bash
# ---------------------------------------------------
# HISAT2 RNA-seq Mapping Script | Batch 1 of 3
# Author: aubsxb005 | Date: 2025-06-21
# ---------------------------------------------------

# Load modules
source /apps/profiles/modules_asax.sh.dyn
module load hisat2/2.2.0
module load samtools/1.18

# Set variables
MYID="aubsxb005"
WD="/scratch/${MYID}/project1_data"
IN_DIR="${WD}/cleaned_data/RNAseq_trimmed_lab/paired"
OUT_DIR="${WD}/mapping/rnaseq_bam"
INDEX="${WD}/reference/rAnoSag1_hisat2_index"
THREADS=8

# Create output directory
mkdir -p "$OUT_DIR"

# Sample list for Batch 2
SAMPLES=(Acy0101 Adis0011 Aeq0050)

# Loop through samples
for SAMPLE in "${SAMPLES[@]}"; do
  echo " Mapping sample: $SAMPLE"

  R1="${IN_DIR}/${SAMPLE}_1_paired.fastq.gz"
  R2="${IN_DIR}/${SAMPLE}_2_paired.fastq.gz"
  SAM="${OUT_DIR}/${SAMPLE}.sam"
  BAM="${OUT_DIR}/${SAMPLE}.bam"
  SORTED="${OUT_DIR}/${SAMPLE}.sorted.bam"
  LOG="${OUT_DIR}/${SAMPLE}_hisat2_summary.txt"

  # Run HISAT2 alignment
  hisat2 -p $THREADS -x "$INDEX" -1 "$R1" -2 "$R2" -S "$SAM" --summary-file "$LOG"

  # Convert SAM → BAM
  samtools view -@ $THREADS -bS "$SAM" > "$BAM"

  # Sort BAM
  samtools sort -@ $THREADS -o "$SORTED" "$BAM"

  # Index BAM
  samtools index "$SORTED"

  # Clean up
  rm "$SAM" "$BAM"

  echo " Finished mapping for: $SAMPLE"
done

echo " Batch 2 RNA-seq mapping complete!"
