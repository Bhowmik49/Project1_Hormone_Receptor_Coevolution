#!/bin/bash
# -------------------------------------------------------------------------
# WGS Batch-6 (SRR4238385 SRR4238386 SRR4238387 SRR4238388 SRR4238389
#              SRR8148313)
# -------------------------------------------------------------------------

set -o pipefail
set -x
source /apps/profiles/modules_asax.sh.dyn
module load samtools/1.18

BED="/scratch/aubsxb005/project1_data/IGF_INS_targets.bed"
OUT="/scratch/aubsxb005/project1_data/mapping/IGF_INS_bam_subsets"
LOG="/scratch/aubsxb005/project1_data/logs/wgs_batch3.log"
mkdir -p "$OUT" "$(dirname "$LOG")"

echo "[START] WGS Batch-3  $(date)" >"$LOG"

declare -a BAMS=(
  "/scratch/aubsxb005/project1_data/mapping/wgs_bam/SRR4238385.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/wgs_bam/SRR4238386.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/wgs_bam/SRR4238387.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/wgs_bam/SRR4238388.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/wgs_bam/SRR4238389.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/wgs_bam/SRR8148313.sorted.bam"
)

for BAM in "${BAMS[@]}"; do
  [[ -f $BAM ]] || { echo "[WARN] missing $BAM" | tee -a "$LOG"; continue; }
  base=$(basename "$BAM" .sorted.bam)
  out="$OUT/${base}_WGS_IGF_INS.bam"
  echo "[INFO] $base" >>"$LOG"
  samtools view -@4 -b -L "$BED" "$BAM" -o "$out"
  samtools index "$out"
  samtools idxstats "$out" | awk '{s+=$3} END{printf("[DONE] %s %d reads\n", "'$base'", s)}' >>"$LOG"
done

echo "[END]  WGS Batch-3  $(date)" >>"$LOG"
