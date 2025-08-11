#!/bin/bash
# Script: download_reference.sh
# Purpose: Download and unpack Anolis sagrei reference genome and annotation files

# Exit on any error
set -e

# Set up working directory
BASE_DIR="/scratch/aubsxb005/project1_data"
REF_DIR="$BASE_DIR/reference"
mkdir -p "$REF_DIR"
cd "$REF_DIR"

# Define base URL (NCBI FTP)
BASE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/037/176/765/GCF_037176765.1_rAnoSag1.mat"

# List of required files
FILES=(
    "GCF_037176765.1_rAnoSag1.mat_genomic.fna.gz"
    "GCF_037176765.1_rAnoSag1.mat_genomic.gtf.gz"
    "GCF_037176765.1_rAnoSag1.mat_genomic.gff.gz"
    "GCF_037176765.1_rAnoSag1.mat_rna.fna.gz"
    "GCF_037176765.1_rAnoSag1.mat_cds_from_genomic.fna.gz"
    "GCF_037176765.1_rAnoSag1.mat_protein.faa.gz"
    "GCF_037176765.1_rAnoSag1.mat_feature_table.txt.gz"
    "GCF_037176765.1_rAnoSag1.mat_assembly_report.txt"
    "GCF_037176765.1_rAnoSag1.mat_assembly_stats.txt"
)

# Download loop
for FILE in "${FILES[@]}"; do
    if [[ -f "$FILE" ]]; then
        echo "File $FILE already exists. Skipping download."
    else
        echo "Downloading $FILE ..."
        wget -q "${BASE_URL}/${FILE}" -O "$FILE"
    fi
done

# Unzip compressed files
echo "Unzipping downloaded .gz files ..."
gunzip -f *.gz

echo "[$(date)] Reference download and extraction complete. Files stored in $REF_DIR"
