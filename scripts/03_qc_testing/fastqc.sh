#!/bin/bash
# Run FastQC on all FASTQ files in the raw_data directory

# Load the appropriate FastQC module
module purge
module load fastqc/0.12.1

# Define paths
INPUT_DIR="/scratch/aubsxb005/project1_data/raw_data"
OUTPUT_DIR="/scratch/aubsxb005/project1_data/fastqc_reports"
LOG_DIR="/scratch/aubsxb005/project1_data/logs"

# Make sure output directories exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# Navigate to input directory
cd "$INPUT_DIR" || { echo " ERROR: Cannot access input directory $INPUT_DIR"; exit 1; }

echo " Running FastQC on all FASTQ files in $INPUT_DIR ..."
for file in *.fastq *.fastq.gz *.fq *.fq.gz; do
    [ -e "$file" ] || continue  # Skip if no matching file
    echo " Processing: $file"
    fastqc "$file" --outdir="$OUTPUT_DIR" --threads 4
done

echo " FastQC run complete! Reports available in: $OUTPUT_DIR"
