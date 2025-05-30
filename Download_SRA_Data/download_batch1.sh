#!/bin/bash
# Download Batch 1: A. insolitus – A. cristatellus
# Author: aubsxb005
# Date: 2025-05-19

module load sra/3.0.0
cd /scratch/aubsxb005/project1_data/raw_data

# Batch 1 SRAs
fastq-dump --split-files SRR4238380  # A. insolitus
fastq-dump --split-files SRR4238381  # A. cybotes
fastq-dump --split-files SRR4238382  # A. occultus
fastq-dump --split-files SRR4238383  # A. lineatopus
fastq-dump --split-files SRR4238384  # A. grahami
fastq-dump --split-files SRR4238385  # A. valencienni
fastq-dump --split-files SRR4238386  # A. sagrei
fastq-dump --split-files SRR4238387  # A. cristatellus

echo "Batch 1 complete."
