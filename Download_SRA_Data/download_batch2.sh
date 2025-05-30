#!/bin/bash
# Download Batch 2: A. evermanni – A. sericeus
# Author: aubsxb005
# Date: 2025-05-19

module load sra/3.0.0
cd /scratch/aubsxb005/project1_data/raw_data

# Batch 2 SRAs
fastq-dump --split-files SRR4238388  # A. evermanni
fastq-dump --split-files SRR4238389  # A. angusticeps
fastq-dump --split-files SRR4238390  # A. chlorocyanus
fastq-dump --split-files SRR31190375 # A. lemurinus
fastq-dump --split-files SRR389085   # A. carolinensis
fastq-dump --split-files SRR23804257 # A. apletophallus
fastq-dump --split-files SRR30501130 # A. tropidonotus
fastq-dump --split-files SRR31253438 # A. planiceps

echo "Batch 2 complete."
