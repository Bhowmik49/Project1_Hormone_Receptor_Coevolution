#!/bin/bash
# Download Batch 8: A. gundlachi – A. brasiliensis
# Author: aubsxb005
# Date: 2025-05-21

module load sra/3.0.0
cd /scratch/aubsxb005/project1_data/raw_data

# Batch 8a SRAs
fastq-dump --split-files SRR24030174 # A. gundlachi
fastq-dump --split-files SRR24030177 # A. poncensis
fastq-dump --split-files SRR13990308 # A. fuscoauratus
fastq-dump --split-files SRR13990414 # A. brasiliensis

echo " Batch 8 complete."
