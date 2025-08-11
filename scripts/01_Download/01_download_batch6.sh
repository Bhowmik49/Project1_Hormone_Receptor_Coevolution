#!/bin/bash
# Download Batch 6: A. distichus favillarum – A. homolechis
# Author: aubsxb005
# Date: 2025-05-21

module load sra/3.0.0
cd /scratch/aubsxb005/project1_data/raw_data

# Batch 5a SRAs
fastq-dump --split-files SRR32222532 # A. distichus favillarum
fastq-dump --split-files ERR3672845  # A. barbatus
fastq-dump --split-files SRR24030180 # A. cuvieri
fastq-dump --split-files DRR360735   # A. homolechis

echo "Batch 6 complete."
