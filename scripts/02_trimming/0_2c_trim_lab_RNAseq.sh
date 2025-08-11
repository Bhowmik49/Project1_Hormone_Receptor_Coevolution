#!/bin/bash
# ----------------------------------------------
# RNA-seq Trimming Script for Lab-generated Dataset (ASC)
# Author: aubsxb005 | Date: 2025-06-21
# Description: Trim adapters and low-quality reads (Phred 20) from lab RNA-seq data
# ----------------------------------------------

# Load necessary modules
source /apps/profiles/modules_asax.sh.dyn
module load trimmomatic/0.39
module load fastqc/0.12.1

# Set paths and variables
MYID="aubsxb005"
WD="/scratch/${MYID}/project1_data"
RAWDIR="${WD}/lab_rna_raw_data"
CLEANDIR="${WD}/cleaned_data/RNAseq_trimmed_lab"
QCPOST="${WD}/fastqc_reports/PostCleanQC_RNAseq_Lab"
ADAPTERS="${WD}/scripts/TruSeq3-PE.fa"
THREADS=6

# Check adapter file
if [[ ! -f "$ADAPTERS" ]]; then
  echo " Adapter file not found at $ADAPTERS"
  exit 1
fi


# Create output directories
mkdir -p "${CLEANDIR}/paired" "${CLEANDIR}/unpaired" "${QCPOST}"

# Navigate to raw directory
cd "$RAWDIR" || { echo " ERROR: lab_rna_raw_data directory not found."; exit 1; }

# Generate sample list by detecting directories (e.g., Acar0321)
ls -d */ | sed 's#/##' > sample_list.txt

# Loop over each sample folder
while read -r SAMPLE; do
    echo " Trimming sample: $SAMPLE"

    R1="${RAWDIR}/${SAMPLE}/${SAMPLE}_1.fq.gz"
    R2="${RAWDIR}/${SAMPLE}/${SAMPLE}_2.fq.gz"

    P1="${CLEANDIR}/paired/${SAMPLE}_1_paired.fastq.gz"
    U1="${CLEANDIR}/unpaired/${SAMPLE}_1_unpaired.fastq.gz"
    P2="${CLEANDIR}/paired/${SAMPLE}_2_paired.fastq.gz"
    U2="${CLEANDIR}/unpaired/${SAMPLE}_2_unpaired.fastq.gz"

    # Run Trimmomatic
    trimmomatic PE -threads $THREADS -phred33 \
        "$R1" "$R2" \
        "$P1" "$U1" "$P2" "$U2" \
        ILLUMINACLIP:"$ADAPTERS":2:30:10 \
        LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36

    # Post-trim QC
    fastqc "$P1" "$P2" --outdir="$QCPOST"

    echo " Finished trimming and QC: $SAMPLE"
done < sample_list.txt

echo " All lab RNA-seq samples trimmed successfully."
