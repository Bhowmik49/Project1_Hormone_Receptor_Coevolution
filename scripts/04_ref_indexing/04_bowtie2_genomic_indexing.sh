#!/bin/bash

# Load bowtie2
module load bowtie2/2.5.1


# Navigate to reference directory
cd /scratch/aubsxb005/project1_data/reference

# Run HISAT2 indexing
bowtie2-build GCF_037176765.1_rAnoSag1.mat_genomic.fna rAnoSag1_bowtie2_index


