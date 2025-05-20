#!/bin/bash
# Download Batch 4: A. frenatus – A. krugi
# Author: aubsxb005
# Date: 2025-05-19

module load sra/3.0.0
cd /scratch/aubsxb005/project1_data/raw_data

# Batch 4 SRAs
fastq-dump --split-files SRR20731678 # A. frenatus
fastq-dump --split-files DRR232286   # A. mestrei
fastq-dump --split-files DRR232284   # A. allisoni
fastq-dump --split-files DRR232285   # A. porcatus
fastq-dump --split-files DRR232282   # A. garridoi
fastq-dump --split-files SRR20731677 # A. isolepis
fastq-dump --split-files DRR232281   # A. alutaceus
fastq-dump --split-files DRR360734   # A. allogus

echo " Batch 4 complete."
