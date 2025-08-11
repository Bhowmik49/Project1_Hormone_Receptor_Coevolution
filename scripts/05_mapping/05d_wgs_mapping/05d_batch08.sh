#!/bin/bash
# Bowtie2 WGS Mapping Batch Script
# Author: aubsxb005
# Description: Maps a batch of 3 WGS samples using Bowtie2

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

# Mapping function for a single sample
map_sample() {
  SAMPLE="$1"
  echo "[$(date)] Processing: $SAMPLE"

  R1="${IN_DIR}/${SAMPLE}_1_paired.fastq"
  R2="${IN_DIR}/${SAMPLE}_2_paired.fastq"
  SAM="${OUT_DIR}/${SAMPLE}.sam"
  BAM="${OUT_DIR}/${SAMPLE}.bam"
  SORTED="${OUT_DIR}/${SAMPLE}.sorted.bam"
  BAI="${SORTED}.bai"
  SUMMARY="${OUT_DIR}/${SAMPLE}_bowtie2_summary.txt"
  LOG="${OUT_DIR}/${SAMPLE}_run.log"

  start_time=$(date +%s)

  if [[ -f "$SORTED" && -f "$BAI" ]]; then
    echo "[$(date)] $SAMPLE already mapped and indexed. Skipping..." | tee -a "$LOG"
    return
  fi

  echo "[$(date)] Running Bowtie2 for $SAMPLE..." | tee -a "$LOG"
  bowtie2 -p "$THREADS" -x "$INDEX" -1 "$R1" -2 "$R2" -S "$SAM" 2> "$SUMMARY"

  if [[ $? -ne 0 ]]; then
    echo "[$(date)] ERROR: Bowtie2 failed for $SAMPLE" | tee -a "$LOG"
    exit 1
  fi

  echo "[$(date)] Converting SAM to BAM..." | tee -a "$LOG"
  samtools view -@ "$THREADS" -bS "$SAM" > "$BAM"

  echo "[$(date)] Sorting BAM..." | tee -a "$LOG"
  samtools sort -@ "$THREADS" -o "$SORTED" "$BAM"

  echo "[$(date)] Indexing sorted BAM..." | tee -a "$LOG"
  samtools index "$SORTED"

  echo "[$(date)] Cleaning up intermediate files..." | tee -a "$LOG"
  rm -f "$SAM" "$BAM"

  end_time=$(date +%s)
  runtime=$((end_time - start_time))

  echo "[$(date)] Finished $SAMPLE in ${runtime}s" | tee -a "$LOG"
}


map_sample "SRR4238388"
map_sample "SRR4238389"
map_sample "SRR8148313"

echo "[$(date)] BATCH08 mapping complete."
