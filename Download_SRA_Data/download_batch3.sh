#!/bin/bash
# Download Batch 3: A. rodriguezii – A. porcatus
# Author: aubsxb005
# Date: 2025-05-19

module load sra/3.0.0
cd /scratch/aubsxb005/project1_data/raw_data

# Batch 3 SRAs
fastq-dump --split-files SRR30361216 # A. rodriguezii
fastq-dump --split-files SRR31186253 # A. sericeus
fastq-dump --split-files SRR21197535 # A. sagrei ordinatus (1)
fastq-dump --split-files SRR21197563 # A. sagrei ordinatus (2)
fastq-dump --split-files SRR25113491 # A. marmoratus
fastq-dump --split-files SRR20731672 # A. baleatus
fastq-dump --split-files SRR20731671 # A. semilineatus
fastq-dump --split-files SRR8148313  # A. auratus

echo " Batch 3 complete."
