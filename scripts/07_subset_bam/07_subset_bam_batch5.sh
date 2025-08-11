#!/bin/bash
# -------------------------------------------------------------------------
# WGS Batch-5 (SRR30361216 SRR30501130 SRR31186253 SRR31190375 SRR31253438
#              SRR4238380 SRR4238381 SRR4238382 SRR4238383 SRR4238384)
# -------------------------------------------------------------------------

set -o pipefail
set -x
source /apps/profiles/modules_asax.sh.dyn
module load samtools/1.18

BED="/scratch/aubsxb005/project1_data/IGF_INS_targets.bed"
OUT="/scratch/aubsxb005/project1_data/mapping/IGF_INS_bam_subsets"
LOG="/scratch/aubsxb005/project1_data/logs/wgs_batch2.log"
mkdir -p "$OUT" "$(dirname "$LOG")"

echo "[START] WGS Batch-2  $(date)" >"$LOG"

declare -a BAMS=(
  "/scratch/aubsxb005/project1_data/mapping/wgs_bam/SRR30361216.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/wgs_bam/SRR30501130.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/wgs_bam/SRR31186253.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/wgs_bam/SRR31190375.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/wgs_bam/SRR31253438.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/wgs_bam/SRR4238380.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/wgs_bam/SRR4238381.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/wgs_bam/SRR4238382.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/wgs_bam/SRR4238383.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/wgs_bam/SRR4238384.sorted.bam"
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

echo "[END]  WGS Batch-2  $(date)" >>"$LOG"
