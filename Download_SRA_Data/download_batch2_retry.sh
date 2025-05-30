#!/bin/bash
# Retry Batch 2 missing samples
# Author: aubsxb005
# Date: 2025-05-20

module load sra/3.0.0
cd /scratch/aubsxb005/project1_data/raw_data

# Re-download only missing ones
fastq-dump --split-files SRR4238390  # A. chlorocyanus
fastq-dump --split-files SRR389085   # A. carolinensis

echo " Retry Batch Complete."
