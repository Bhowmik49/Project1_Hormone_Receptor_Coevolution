#!/bin/bash
# Purpose: Split full BED into per-gene files for CDS and exon features only, removing duplicates safely

set -o pipefail
set -x

source /apps/profiles/modules_asax.sh.dyn
module load bedtools/2.26.0

BED_IN="/scratch/aubsxb005/project1_data/IGF_INS_targets.bed"
OUT_DIR="/scratch/aubsxb005/project1_data/bed_files_per_gene_by_feature"
mkdir -p "$OUT_DIR"

# Extract unique gene names from the BED's fourth column
cut -f4 "$BED_IN" | cut -d"_" -f1 | sort | uniq | while read GENE; do
    for FEATURE in CDS exon; do
        OUT_FILE="${OUT_DIR}/${GENE}_${FEATURE}.bed"

        awk -v gene="$GENE" -v feature="$FEATURE" '
            $4 == gene"_"feature {
                print $1"\t"$2"\t"$3
            }
        ' "$BED_IN" | sort -k1,1 -k2,2n | uniq > "$OUT_FILE"

        if [[ -s "$OUT_FILE" ]]; then
            echo "[INFO] Wrote clean BED: $OUT_FILE"
        else
            echo "[WARNING] No entries found for ${GENE}_${FEATURE}, file skipped."
            rm -f "$OUT_FILE"
        fi
    done
done

echo "[DONE] All clean per-gene, per-feature BEDs saved in: $OUT_DIR/"
