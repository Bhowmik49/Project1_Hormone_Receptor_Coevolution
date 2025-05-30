#!/bin/bash
#PBS -N qc_batch2
#PBS -q medium
#PBS -l select=1:ncpus=4:mem=8gb
#PBS -l walltime=04:00:00
#PBS -j oe
#PBS -o /scratch/aubsxb005/project1_data/logs/fastqc_batch2.out

# Load necessary modules
module purge
module load fastqc/0.12.1
module load multiqc/1.15

# Define directories 
INPUT_DIR="/scratch/aubsxb005/project1_data/raw_data"
FASTQC_DIR="/scratch/aubsxb005/project1_data/fastqc_reports/batch2"
MULTIQC_DIR="/scratch/aubsxb005/project1_data/multiqc_reports/batch2"
LOG_DIR="/scratch/aubsxb005/project1_data/logs"

mkdir -p "$FASTQC_DIR" "$MULTIQC_DIR" "$LOG_DIR"


# Move to input directory 
cd "$INPUT_DIR" || exit 1

echo "[$(date)] Starting FastQC for Batch2..."

# Run FastQC on each paired sample, if files exist
for SAMPLE in SRR13990293 SRR13990308 SRR13990391 SRR13990414 SRR20731671 SRR20731672 SRR20731677 SRR20731678; do
    if [[ -f "${SAMPLE}_1.fastq" && -f "${SAMPLE}_2.fastq" ]]; then
        echo "  Processing $SAMPLE..."
        fastqc "${SAMPLE}_1.fastq" "${SAMPLE}_2.fastq" --outdir="$FASTQC_DIR" --threads 4
    else
        echo "  [Warning] Files for $SAMPLE missing. Skipping."
    fi
done

# Run MultiQC to summarize FastQC results
echo "[$(date)] Running MultiQC for Batch2..."
multiqc "$FASTQC_DIR" -o "$MULTIQC_DIR" --filename "multiqc_batch2_report"

echo "[$(date)] Batch2 FastQC + MultiQC completed."
