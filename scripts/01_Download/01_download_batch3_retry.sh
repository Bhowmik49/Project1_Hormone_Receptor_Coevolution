#!/bin/bash
# Retry Batch 3: Remaining samples from A. marmoratus – A. auratus
# Author: aubsxb005
# Date: 2025-05-21

module load sra/3.0.0
cd /scratch/aubsxb005/project1_data/raw_data

# Re-download failed ones
fastq-dump --split-files SRR25113491 # A. marmoratus
fastq-dump --split-files SRR20731672 # A. baleatus
fastq-dump --split-files SRR20731671 # A. semilineatus
fastq-dump --split-files SRR8148313  # A. auratus

echo "Retry of Batch 3 complete."
