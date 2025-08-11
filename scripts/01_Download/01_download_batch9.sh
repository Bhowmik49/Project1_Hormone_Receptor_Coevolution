#!/bin/bash
# Download Batch 9: A. ortonii – A. dissimilis
# Author: aubsxb005
# Date: 2025-05-21

module load sra/3.0.0
cd /scratch/aubsxb005/project1_data/raw_data

# Batch 8b SRAs
fastq-dump --split-files SRR13990293 # A. ortonii
fastq-dump --split-files SRR13990391 # A. chrysolepis
fastq-dump --split-files SRR7884574  # A. tandai
fastq-dump --split-files SRR7884561  # A. dissimilis

echo "Batch 9 complete."





