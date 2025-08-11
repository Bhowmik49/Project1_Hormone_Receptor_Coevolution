#!/bin/bash
# Download Batch 4: A. frenatus – A. porcatus
# Author: aubsxb005
# Date: 2025-05-21

module load sra/3.0.0
cd /scratch/aubsxb005/project1_data/raw_data

# Batch 4a SRAs
fastq-dump --split-files SRR20731678 # A. frenatus
fastq-dump --split-files DRR232286   # A. mestrei
fastq-dump --split-files DRR232284   # A. allisoni
fastq-dump --split-files DRR232285   # A. porcatus

echo "Batch 4 complete."



