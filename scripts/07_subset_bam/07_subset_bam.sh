#!/bin/bash
# -------------------------------------------------------------------------
# Script: 05c_subset_bams_test_batch.sh
# Purpose: Subset test batch of RNA-seq and WGS BAMs using IGF/INS BED regions
# Platform: Alabama Supercomputer (ASC)
# Author: Sagar Bhowmik Ocean
# -------------------------------------------------------------------------

set -o pipefail
set -x

# Load ASC environment
source /apps/profiles/modules_asax.sh.dyn
module load samtools/1.18

# Define paths
BED="/scratch/aubsxb005/project1_data/IGF_INS_targets.bed"
OUT_DIR="/scratch/aubsxb005/project1_data/mapping/IGF_INS_bam_subsets_test"
LOG="/scratch/aubsxb005/project1_data/logs/bam_subsetting_test.log"

# Create output and log directories
mkdir -p "$OUT_DIR"
mkdir -p "$(dirname "$LOG")"

echo "[INFO] Test BAM Subsetting Started: $(date)" > "$LOG"

# Define test BAMs (paths verified)
WGS_BAMS=(
  "/scratch/aubsxb005/project1_data/mapping/wgs_bam/DRR360734.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/wgs_bam/SRR24030175.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/wgs_bam/SRR25113491.sorted.bam"
)

RNA_BAMS=(
  "/scratch/aubsxb005/project1_data/mapping/rnaseq_bam/Acar0321.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/rnaseq_bam/Acr0023.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/rnaseq_bam/Adis0011.sorted.bam"
)

# Function to process BAMs
subset_bams() {
  local TYPE="$1"
  shift
  local BAM_LIST=("$@")

  for BAM in "${BAM_LIST[@]}"; do
    if [[ ! -f "$BAM" ]]; then
      echo "[WARNING] File not found: $BAM" | tee -a "$LOG"
      continue
    fi

    BASENAME=$(basename "$BAM" .sorted.bam)
    OUT_BAM="${OUT_DIR}/${BASENAME}_${TYPE}_IGF_INS.bam"

    echo "[INFO] Processing $BASENAME ($TYPE)" | tee -a "$LOG"
    echo "[CMD] samtools view -@ 4 -b -L $BED $BAM -o $OUT_BAM" >> "$LOG"

    if samtools view -@ 4 -b -L "$BED" "$BAM" -o "$OUT_BAM"; then
      samtools index "$OUT_BAM"
      COUNT=$(samtools view "$OUT_BAM" | wc -l)
      echo "[SUCCESS] $BASENAME: $COUNT reads extracted" >> "$LOG"
    else
      echo "[ERROR] Failed to subset $BASENAME" >> "$LOG"
      rm -f "$OUT_BAM"
    fi
    echo "---" >> "$LOG"
  done
}

# Run subsetting for both BAM groups
subset_bams "WGS" "${WGS_BAMS[@]}"
subset_bams "RNAseq" "${RNA_BAMS[@]}"

echo "[INFO] Test BAM Subsetting Completed: $(date)" >> "$LOG"
