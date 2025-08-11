#!/bin/bash
# -------------------------------------------------------------------------
# RNA-seq Batch-3  (DRR232282 DRR232284 DRR232285 DRR232286 SRR389085)
# -------------------------------------------------------------------------

set -o pipefail
set -x
source /apps/profiles/modules_asax.sh.dyn
module load samtools/1.18

BED="/scratch/aubsxb005/project1_data/IGF_INS_targets.bed"
OUT="/scratch/aubsxb005/project1_data/mapping/IGF_INS_bam_subsets"
LOG="/scratch/aubsxb005/project1_data/logs/rnaseq_batch3.log"
mkdir -p "$OUT" "$(dirname "$LOG")"

echo "[START] RNA-seq Batch-3  $(date)" >"$LOG"

declare -a BAMS=(
  "/scratch/aubsxb005/project1_data/mapping/rnaseq_bam/DRR232282.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/rnaseq_bam/DRR232284.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/rnaseq_bam/DRR232285.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/rnaseq_bam/DRR232286.sorted.bam"
  "/scratch/aubsxb005/project1_data/mapping/rnaseq_bam/SRR389085.sorted.bam"
)

for BAM in "${BAMS[@]}"; do
  [[ -f $BAM ]] || { echo "[WARN] missing $BAM" | tee -a "$LOG"; continue; }
  base=$(basename "$BAM" .sorted.bam)
  out="$OUT/${base}_RNAseq_IGF_INS.bam"
  echo "[INFO] $base" >>"$LOG"
  samtools view -@4 -b -L "$BED" "$BAM" -o "$out"
  samtools index "$out"
  samtools idxstats "$out" | awk '{s+=$3} END{printf("[DONE] %s %d reads\n", "'$base'", s)}' >>"$LOG"
done

echo "[END]  RNA-seq Batch-3  $(date)" >>"$LOG"
