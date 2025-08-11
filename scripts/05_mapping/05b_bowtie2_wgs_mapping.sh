#!/bin/bash
# ---------------------------------------------------
# Bowtie2 WGS Mapping Script | Brown Anole Genome
# Author: aubsxb005 | Last updated: 2025-06-11
# Description: Align paired-end WGS reads to rAnoSag1 genome using Bowtie2 and sort BAM with SAMtools
# ---------------------------------------------------

# Load environments
source /apps/profiles/modules_asax.sh.dyn
module load bowtie2/2.5.1
module load samtools/1.18

# Variables
MYID="aubsxb005"
WD="/scratch/$MYID/project1_data"
IN_DIR="${WD}/cleaned_data/WGS_trimmed/paired"
OUT_DIR="${WD}/mapping/wgs_bam"
INDEX="${WD}/reference/rAnoSag1_bowtie2_index"
THREADS=8

# Create output directory
mkdir -p "$OUT_DIR"

# Move to input directory
cd "$IN_DIR" || { echo "Input directory not found"; exit 1; }

# Generate sample list
ls *_1_paired.fastq | sed 's/_1_paired.fastq//' > sample_list.txt

# Loop for mapping
while read -r SAMPLE; do
  echo "Processing sample: $SAMPLE"

  R1="${IN_DIR}/${SAMPLE}_1_paired.fastq"
  R2="${IN_DIR}/${SAMPLE}_2_paired.fastq"
  SAM="${OUT_DIR}/${SAMPLE}.sam"
  BAM="${OUT_DIR}/${SAMPLE}.bam"
  SORTED="${OUTDIR}/${SAMPLE}.sorted.bam"

  # Align using Bowtie2
  bowtie2 -p "$THREADS" -x "$INDEX" -1 "$R1" -2 "$R2" -S "$SAM" 2> "${OUT_DIR}/${SAMPLE}_bowtie2_summary.txt"

  # Convert SAM to BAM, then sort and index
  samtools view -@ "$THREADS" -bS "$SAM" > "$BAM"
  samtools sort -@ "$THREADS" -o "$SORTED" "$BAM"
  samtools index "$SORTED"

  # Clean up intermediate files
  rm "$SAM" "$BAM"

  echo "Finished sample: $SAMPLE"
done < sample_list.txt

echo "WGS alignment and sorting complete."
