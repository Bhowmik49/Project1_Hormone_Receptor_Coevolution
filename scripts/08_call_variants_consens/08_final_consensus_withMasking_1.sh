#!/bin/bash

# Exit on error and echo commands
set -o pipefail
set -x

# --------------------------------------
# Load necessary ASC modules
# --------------------------------------
source /apps/profiles/modules_asax.sh.dyn
module load bcftools/1.13
module load samtools/1.18
module load bedtools/2.26.0

# --------------------------------------
# Define paths
# --------------------------------------
ROOT="/scratch/aubsxb005/project1_data/project1_working_dataset"
REF="$ROOT/reference/GCF_037176765.1_rAnoSag1.mat_genomic.fna"
BAM_DIR="$ROOT/mapping/IGF_INS_bam_subsets"
VCF_DIR="$ROOT/variants/vcfs_truthful_consensus"
MASK_DIR="$ROOT/temp_metadata/depth_masks"
CONS_DIR="$ROOT/variants/consensus_fastas_5"
LOG="$ROOT/logs/consensus_gen_v5_$(date +%Y%m%d_%H%M%S).log"

mkdir -p "$VCF_DIR" "$MASK_DIR" "$CONS_DIR" "$(dirname "$LOG")"

echo "[START] Consensus FASTA generation: $(date)" > "$LOG"

# --------------------------------------
# Count input BAMs
# --------------------------------------
BAM_COUNT=$(ls "$BAM_DIR"/*.bam 2>/dev/null | wc -l)
echo "[INFO] Total BAM files found: $BAM_COUNT" | tee -a "$LOG"

# --------------------------------------
# Loop over each BAM file
# --------------------------------------
for BAM in "$BAM_DIR"/*.bam; do
    if [[ -f "$BAM" ]]; then
        echo "[INFO] Processing BAM: $BAM" | tee -a "$LOG"

        BASENAME=$(basename "$BAM" .bam)
        SAMPLE=${BASENAME%%_IGF_INS*}

        VCF_OUT="$VCF_DIR/${SAMPLE}.vcf.gz"
        BED_MASK="$MASK_DIR/${SAMPLE}_mask_depth.bed"
        CONS_OUT="$CONS_DIR/${SAMPLE}_truthful_consensus.fasta"

        echo "[INFO] Sample: $SAMPLE" | tee -a "$LOG"

        # --------------------------------------
        # Step 1: Call variants (required for consensus)
        # --------------------------------------
        bcftools mpileup -Ou -f "$REF" "$BAM" | \
        bcftools call -c -Oz -o "$VCF_OUT"
        bcftools index "$VCF_OUT"

        # --------------------------------------
        # Step 2: Generate BED for low/high coverage regions (DP <1 or >100)
        # --------------------------------------
        samtools depth "$BAM" | \
        awk '{ if ($3 < 1 || $3 > 100) print $1"\t"$2-1"\t"$2 }' > "$BED_MASK"

        # --------------------------------------
        # Step 3: Generate N-masked consensus FASTA
        # --------------------------------------
        bcftools consensus -f "$REF" -m "$BED_MASK" "$VCF_OUT" > "$CONS_OUT"

        # --------------------------------------
        # Step 4: Index the FASTA file (.fai)
        # --------------------------------------
        if [[ -s "$CONS_OUT" ]]; then
            samtools faidx "$CONS_OUT"
            echo "[DONE] $SAMPLE → FASTA + FAI created" | tee -a "$LOG"
        else
            echo "[WARNING] Empty FASTA generated for $SAMPLE!" | tee -a "$LOG"
            rm -f "$CONS_OUT"  # Optionally remove failed output
        fi
    fi
done

echo "[END] All consensus generation complete: $(date)" | tee -a "$LOG"
