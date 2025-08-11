#!/bin/bash
# ----------------------------------------------
# MultiQC Report for Cleaned RNA-seq and WGS Data
# Author: aubsxb005 | Last updated: 2025-06-11
# Description: Summarize post-trimming FastQC outputs using MultiQC
# ----------------------------------------------

# Load environment
source /apps/profiles/modules_asax.sh.dyn
module load multiqc/1.15

# Define directories
MYID="aubsxb005"
WD="/scratch/$MYID/project1_data"
POSTQC_RNA="${WD}/fastqc_reports/PostCleanQC_RNAseq"
POSTQC_WGS="${WD}/fastqc_reports/PostCleanQC_WGS"
OUTDIR="${WD}/multiqc_reports"


# Create output subdirectories
mkdir -p "${OUTDIR}/RNAseq" "${OUTDIR}/WGS"

# Run MultiQC on RNA-seq QC outputs
multiqc "$POSTQC_RNA" --outdir "${OUTDIR}/RNAseq" --force --title "Post-Trim RNA-seq Quality Report"

# Run MultiQC on WGS QC outputs
multiqc "$POSTQC_WGS" --outdir "${OUTDIR}/WGS" --force --title "Post-Trim WGS Quality Report"

echo " MultiQC reports saved in: ${OUTDIR}/RNAseq and ${OUTDIR}/WGS"