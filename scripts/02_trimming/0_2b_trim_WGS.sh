#!/bin/bash
# ----------------------------------------------
# WGS Trimming Script for Anolis Dataset (ASC)
# Author: aubsxb005 | Last updated: 2025-06-10
# Description: Uses Trimmomatic to trim adapters and low-quality bases
# ----------------------------------------------

# Load environment
source /apps/profiles/modules_asax.sh.dyn
module load trimmomatic/0.39
module load fastqc/0.12.1

# Variables
MYID="aubsxb005"
WD="/scratch/$MYID/project1_data"
RAWDIR="${WD}/Genomic_data/WGS_data"
CLEANDIR="${WD}/cleaned_data/WGS_trimmed"
QCPOST="${WD}/fastqc_reports/PostCleanQC_WGS"
ADAPTERS="${WD}/scripts/TruSeq3-PE.fa"
THREADS=6

# Check adapter file
if [[ ! -f "$ADAPTERS" ]]; then
  echo " Adapter file not found at $ADAPTERS"
  exit 1
fi

# Create output directories
mkdir -p "${CLEANDIR}/paired" "${CLEANDIR}/unpaired" "${QCPOST}"

# Move to input directory and generate sample list
cd "$RAWDIR" || { echo " WGS_data folder not found!"; exit 1; }
ls *_1.fastq | sed 's/_1.fastq//' > sample_list.txt

# Main loop
while read -r SAMPLE; do
    echo " Trimming WGS sample: ${SAMPLE}"

    R1="${RAWDIR}/${SAMPLE}_1.fastq"
    R2="${RAWDIR}/${SAMPLE}_2.fastq"
    P1="${CLEANDIR}/paired/${SAMPLE}_1_paired.fastq"
    U1="${CLEANDIR}/unpaired/${SAMPLE}_1_unpaired.fastq"
    P2="${CLEANDIR}/paired/${SAMPLE}_2_paired.fastq"
    U2="${CLEANDIR}/unpaired/${SAMPLE}_2_unpaired.fastq"

    # Trimming
    trimmomatic PE -threads $THREADS -phred33 \
      "$R1" "$R2" \
      "$P1" "$U1" "$P2" "$U2" \
      ILLUMINACLIP:"$ADAPTERS":2:30:10 \
      LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36

    # Post-trim QC
    fastqc "$P1" "$P2" --outdir="${QCPOST}"

    echo " Finished: ${SAMPLE}"
done < sample_list.txt

echo " WGS trimming complete for all samples."
