#!/bin/bash
# Extract gene coordinates for IGF/INS family genes from GFF

source /apps/profiles/modules_asax.sh.dyn
set -euo pipefail

# Variables
GENES=("IGF1" "IGF2" "INS" "IGF1R" "IGF2R" "INSR")
GFF="reference/GCF_037176765.1_rAnoSag1.mat_genomic.gff"
OUT="IGF_INS_gene_coords.txt"

# Create output file
echo -e "chrom\tstart\tend\tfeature_type\tgene_name\tstrand" > "$OUT"

# Extract entries for each gene
for GENE in "${GENES[@]}"; do
    grep -iw "$GENE" "$GFF" | \
    awk -v gene="$GENE" '$3 ~ /CDS|exon|transcript/ {
        match($0, /[ \t]strand=([-+])/, m);
        strand = (m[1] != "") ? m[1] : ".";
        print $1, $4, $5, $3, gene, strand;
    }' OFS="\t" >> "$OUT"
done

echo "[INFO] Gene coordinate extraction completed: $OUT"

