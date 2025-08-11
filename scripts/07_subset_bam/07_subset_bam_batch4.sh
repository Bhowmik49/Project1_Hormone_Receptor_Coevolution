#!/bin/bash
# -------------------------------------------------------------------------
# WGS Batch-4 (DRR360734 DRR360735 SRR21197535 SRR24030174 SRR24030175
#              SRR24030177 SRR24030180 SRR24030181 SRR25113491 SRR2651655)
# -------------------------------------------------------------------------

set -o pipefail
set -x
source /apps/profiles/modules_asax.sh.dyn
module load samtools/1.18

BED="/scratch/aubsxb005/project1_data/IGF_INS_targets.bed"
OUT="/scratch/aubsxb005/project1_data/mapping/IGF_INS_bam_subsets"
LOG="/scratch/aubsxb005/project1_data/logs/wgs_batch1.log"
mkdir -p "$OUT" "$(dirname "$LOG")"

echo "[START] WGS Batch-1  $(date)" >"$LOG"

declare -a BAMS=(
  "/scratch/aubsxb005/project1_data/mapping/wgs_bam/DRR360734.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/wgs_bam/DRR360735.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/wgs_bam/SRR21197535.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/wgs_bam/SRR24030174.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/wgs_bam/SRR24030175.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/wgs_bam/SRR24030177.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/wgs_bam/SRR24030180.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/wgs_bam/SRR24030181.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/wgs_bam/SRR25113491.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/wgs_bam/SRR2651655.sorted.bam"
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

echo "[END]  WGS Batch-1  $(date)" >>"$LOG"
