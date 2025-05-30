#!/bin/bash
# Retry Batch 3: SRR20731671 (A. semilineatus)
# Author: aubsxb005
# Date: 2025-05-21

module load sra/3.0.0
cd /scratch/aubsxb005/project1_data/raw_data

fastq-dump --split-files SRR20731671  # A. semilineatus

echo " SRR20731671 download complete."
