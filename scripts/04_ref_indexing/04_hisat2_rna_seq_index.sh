#!/bin/bash

# Load HISAT2
module load hisat2/2.2.0

# Navigate to reference directory
cd /scratch/aubsxb005/project1_data/reference

# Run HISAT2 indexing
hisat2-build GCF_037176765.1_rAnoSag1.mat_genomic.fna rAnoSag1_hisat2_index
