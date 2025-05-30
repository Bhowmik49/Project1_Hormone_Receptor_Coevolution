#!/bin/bash
#PBS -N qc_batch6
#PBS -q medium
#PBS -l select=1:ncpus=4:mem=8gb
#PBS -l walltime=04:00:00
#PBS -j oe
#PBS -o /scratch/aubsxb005/project1_data/logs/fastqc_batch6.out

# Load necessary modules
module purge
module load fastqc/0.12.1
module load multiqc/1.15


# Define directories 
INPUT_DIR="/scratch/aubsxb005/project1_data/raw_data"
FASTQC_DIR="/scratch/aubsxb005/project1_data/fastqc_reports/batch6"
MULTIQC_DIR="/scratch/aubsxb005/project1_data/multiqc_reports/batch6"
LOG_DIR="/scratch/aubsxb005/project1_data/logs"

mkdir -p "$FASTQC_DIR" "$MULTIQC_DIR" "$LOG_DIR"

# Move to input directory 
cd "$INPUT_DIR" || exit 1

echo "[$(date)] Starting FastQC for Batch6..."

for SAMPLE in SRR4238386 SRR4238387 SRR4238388 SRR4238389 SRR4238390 SRR7884561 SRR7884574 SRR8148313; do
    if [[ -f "${SAMPLE}_1.fastq" && -f "${SAMPLE}_2.fastq" ]]; then
        echo "  Processing $SAMPLE..."
        fastqc "${SAMPLE}_1.fastq" "${SAMPLE}_2.fastq" --outdir="$FASTQC_DIR" --threads 4
    else
        echo "  [Warning] Files for $SAMPLE missing. Skipping."
    fi
done

# Run MultiQC to summarize FastQC results
echo "[$(date)] Running MultiQC for Batch6..."
multiqc "$FASTQC_DIR" -o "$MULTIQC_DIR" --filename "multiqc_batch6_report"

echo "[$(date)] Batch6 FastQC + MultiQC completed."


