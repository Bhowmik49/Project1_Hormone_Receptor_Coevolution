#!/bin/bash
# ---------------------------------------------------
# HISAT2 RNA-seq Mapping Script | Brown Anole Genome
# Author: aubsxb005 | Last updated: 2025-06-11
# Description: Align paired-end RNA-seq reads to rAnoSag1 genome using HISAT2 and sort/index BAM with SAMtools
# ---------------------------------------------------

# Load environments
source /apps/profiles/modules_asax.sh.dyn
module load hisat2/2.2.0
module load samtools/1.18

# Variables
MYID="aubsxb005"
WD="/scratch/$MYID/project1_data"
IN_DIR="${WD}/cleaned_data/RNAseq_trimmed/paired"
OUT_DIR="${WD}/mapping/rnaseq_bam"
INDEX="${WD}/reference/rAnoSag1_hisat2_index"  # built from GCF_037176765.1_rAnoSag1.mat_genomic.fna
THREADS=8

# Create output directory
mkdir -p "$OUT_DIR"

# Move to input directory
cd "$IN_DIR" || { echo "Error: Input directory not found!"; exit 1; }

# Generate sample list
ls *_1_paired.fastq | sed 's/_1_paired.fastq//' > sample_list.txt

# Alignment loop
while read -r SAMPLE; do
  echo ">> Mapping sample: $SAMPLE"

  R1="${IN_DIR}/${SAMPLE}_1_paired.fastq"
  R2="${IN_DIR}/${SAMPLE}_2_paired.fastq"
  SAM="${OUT_DIR}/${SAMPLE}.sam"
  BAM="${OUT_DIR}/${SAMPLE}.bam"
  SORTED="${OUT_DIR}/${SAMPLE}.sorted.bam"
  LOG="${OUT_DIR}/${SAMPLE}_hisat2_summary.txt"

  # HISAT2 alignment
  hisat2 -p $THREADS -x "$INDEX" -1 "$R1" -2 "$R2" -S "$SAM" --summary-file "$LOG"

  # SAM → BAM
  samtools view -@ $THREADS -bS "$SAM" -o "$BAM"

  # Sort BAM
  samtools sort -@ $THREADS -o "$SORTED" "$BAM"

  # Index BAM
  samtools index "$SORTED"

  # Clean up
  rm "$SAM" "$BAM"

  echo " Finished sample: $SAMPLE"
done < sample_list.txt

echo " RNA-seq alignment and sorting complete for all samples."
