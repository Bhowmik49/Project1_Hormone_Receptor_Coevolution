#!/bin/bash
# Run FastQC on all FASTQ files in the raw_data directory

# Set paths
WORKDIR="/scratch/aubsxb005/project1_data"
RAWDIR="$WORKDIR/raw_data"
OUTDIR="$WORKDIR/fastqc_reports"

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Load FastQC module (adjust if your ASC module name differs)
module load fastqc/0.11.9

# Run FastQC on all FASTQ files
for file in "$RAWDIR"/*.fastq
do
    echo "Running FastQC on $file"
    fastqc "$file" -o "$OUTDIR"
done

echo "FastQC completed. Reports saved to $OUTDIR"
