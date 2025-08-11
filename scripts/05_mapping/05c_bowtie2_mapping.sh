#!/bin/bash
# ---------------------------------------------------
# Bowtie2 WGS Mapping Script | Brown Anole Genome
# Author: aubsxb005 | Revised: 2025-06-16
# Description: Align paired-end WGS reads to rAnoSag1 genome using Bowtie2
# ---------------------------------------------------

# Load environment modules
source /apps/profiles/modules_asax.sh.dyn
module load bowtie2/2.5.1
module load samtools/1.18

# Define paths and variables
MYID="aubsxb005"
WD="/scratch/$MYID/project1_data"
IN_DIR="${WD}/cleaned_data/WGS_trimmed/paired"
OUT_DIR="${WD}/mapping/wgs_bam"
INDEX="${WD}/reference/rAnoSag1_bowtie2_index"
THREADS=8

# Create output directory if not present
mkdir -p "$OUT_DIR"

# Move to input directory
cd "$IN_DIR" || { echo " Input directory not found!"; exit 1; }

# Generate sample list
ls *_1_paired.fastq | sed 's/_1_paired.fastq//' > sample_list.txt

# Loop through each sample for mapping 
while read -r SAMPLE; do
  echo " Processing: $SAMPLE"
  
  R1="${IN_DIR}/${SAMPLE}_1_paired.fastq"
  R2="${IN_DIR}/${SAMPLE}_2_paired.fastq"
  SAM="${OUT_DIR}/${SAMPLE}.sam"
  BAM="${OUT_DIR}/${SAMPLE}.bam"
  SORTED="${OUT_DIR}/${SAMPLE}.sorted.bam"
  SUMMARY="${OUT_DIR}/${SAMPLE}_bowtie2_summary.txt"

  # Skip if sorted BAM already exists
  if [[ -f "$SORTED" ]]; then
    echo " $SAMPLE already mapped. Skipping..."
    continue
  fi

 # Align using Bowtie2
  echo " Mapping $SAMPLE..."
  bowtie2 -p "$THREADS" -x "$INDEX" -1 "$R1" -2 "$R2" -S "$SAM" 2> "$SUMMARY"

  # Convert SAM to BAM, then sort and index
  samtools view -@ "$THREADS" -bS "$SAM" > "$BAM"
  samtools sort -@ "$THREADS" -o "$SORTED" "$BAM"
  samtools index "$SORTED"

  # Clean up intermediate files 
  rm -f "$SAM" "$BAM"

  echo " Finished: $SAMPLE"
done < sample_list.txt

echo " All WGS mapping complete."




### (Choose 32 hours, 48 GB, and 16 cores just like before.)





