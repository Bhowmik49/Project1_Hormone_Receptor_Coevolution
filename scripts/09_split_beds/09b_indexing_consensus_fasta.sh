#!/bin/bash
#=====================================================
# STEP 6: Index consensus FASTAs (*.fasta → *.fasta.fai)
# PURPOSE: Prepares for region-based sequence extraction
#=====================================================

set -o pipefail
set -x

#======== ASC ENVIRONMENT SETUP ======================
source /apps/profiles/modules_asax.sh.dyn
module load samtools/1.18

#======== DIRECTORY SETUP ============================
CONS_DIR="/scratch/aubsxb005/project1_data/project1_working_dataset/variants/consensus_Nmasked_truthful"
cd "$CONS_DIR"

echo "[INFO] Indexing FASTAs in: $CONS_DIR"

#======== INDEXING LOOP ==============================
for fasta in *_truthful_consensus.fasta; do
    if [[ -f "$fasta" ]]; then
        echo "[INFO] Indexing: $fasta"
        samtools faidx "$fasta"
    else
        echo "[WARNING] No matching FASTA files found!"
        break
    fi
done

echo "[DONE] All consensus FASTAs are indexed."
