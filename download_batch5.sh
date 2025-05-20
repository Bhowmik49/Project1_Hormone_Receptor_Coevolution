#!/bin/bash
# Download Batch 5: A. distichus favillarum – A. stratulus
# Author: aubsxb005
# Date: 2025-05-19

module load sra/3.0.0
cd /scratch/aubsxb005/project1_data/raw_data

# Batch 5 SRAs
fastq-dump --split-files SRR32222532 # A. distichus favillarum
fastq-dump --split-files ERR3672845  # A. barbatus
fastq-dump --split-files SRR24030180 # A. cuvieri
fastq-dump --split-files DRR360735   # A. homolechis
fastq-dump --split-files SRR24030181 # A. krugi
fastq-dump --split-files SRR2651655  # A. loysiana
fastq-dump --split-files SRR26324123 # A. pulchellus
fastq-dump --split-files SRR24030175 # A. stratulus

echo " Batch 5 complete."
