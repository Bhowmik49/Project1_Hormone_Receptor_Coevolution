#!/bin/bash
# Download Batch 5: A. garridoi – A. allogus
# Author: aubsxb005
# Date: 2025-05-21

module load sra/3.0.0
cd /scratch/aubsxb005/project1_data/raw_data

# Batch 4b SRAs
fastq-dump --split-files DRR232282   # A. garridoi
fastq-dump --split-files SRR20731677 # A. isolepis
fastq-dump --split-files DRR232281   # A. alutaceus
fastq-dump --split-files DRR360734   # A. allogus

echo " Batch 5 complete."
