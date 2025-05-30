#!/bin/bash
#PBS -N qc_batch3
#PBS -q medium
#PBS -l select=1:ncpus=4:mem=8gb
#PBS -l walltime=04:00:00
#PBS -j oe
#PBS -o /scratch/aubsxb005/project1_data/logs/fastqc_batch3.out


# Load necessary modules
module purge
module load fastqc/0.12.1
module load multiqc/1.15

# Define directories 
INPUT_DIR="/scratch/aubsxb005/project1_data/raw_data"
FASTQC_DIR="/scratch/aubsxb005/project1_data/fastqc_reports/batch3"
MULTIQC_DIR="/scratch/aubsxb005/project1_data/multiqc_reports/batch3"
LOG_DIR="/scratch/aubsxb005/project1_data/logs"

mkdir -p "$FASTQC_DIR" "$MULTIQC_DIR" "$LOG_DIR"

# Move to input directory 
cd "$INPUT_DIR" || exit 1

echo "[$(date)] Starting FastQC for Batch3..."

for SAMPLE in SRR21197535 SRR21197563 SRR23804257 SRR24030174 SRR24030175 SRR24030177 SRR24030180 SRR24030181; do
    if [[ -f "${SAMPLE}_1.fastq" && -f "${SAMPLE}_2.fastq" ]]; then
        echo "  Processing $SAMPLE..."
        fastqc "${SAMPLE}_1.fastq" "${SAMPLE}_2.fastq" --outdir="$FASTQC_DIR" --threads 4
    else
        echo "  [Warning] Files for $SAMPLE missing. Skipping."
    fi
done

# Run MultiQC to summarize FastQC results
echo "[$(date)] Running MultiQC for Batch3..."
multiqc "$FASTQC_DIR" -o "$MULTIQC_DIR" --filename "multiqc_batch3_report"

echo "[$(date)] Batch3 FastQC + MultiQC completed."
