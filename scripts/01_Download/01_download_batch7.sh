#!/bin/bash
# Download Batch 7: A. krugi – A. stratulus
# Author: aubsxb005
# Date: 2025-05-21

module load sra/3.0.0
cd /scratch/aubsxb005/project1_data/raw_data

# Batch 5b SRAs
fastq-dump --split-files SRR24030181 # A. krugi
fastq-dump --split-files SRR2651655  # A. loysiana
fastq-dump --split-files SRR26324123 # A. pulchellus
fastq-dump --split-files SRR24030175 # A. stratulus

echo " Batch 7 complete."


