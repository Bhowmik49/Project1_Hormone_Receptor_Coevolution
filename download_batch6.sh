#!/bin/bash
# Download Batch 6: A. gundlachi – A. chrysolepis
# Author: aubsxb005
# Date: 2025-05-19

module load sra/3.0.0
cd /scratch/aubsxb005/project1_data/raw_data

# Batch 6 SRAs
fastq-dump --split-files SRR24030174 # A. gundlachi
fastq-dump --split-files SRR24030177 # A. poncensis
fastq-dump --split-files SRR13990308 # A. fuscoauratus
fastq-dump --split-files SRR13990414 # A. brasiliensis
fastq-dump --split-files SRR13990293 # A. ortonii
fastq-dump --split-files SRR13990391 # A. chrysolepis
fastq-dump --split-files SRR7884574  # A. tandai
fastq-dump --split-files SRR7884561  # A. dissimilis

echo " Batch 6 complete."
