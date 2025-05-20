#!/bin/bash

#==============================
# Script to Download SRA Data
# Project: Anolis IGFs
# Author: aubsxb005
# Date: 05-18-2025
#==============================

# Load SRA toolkit
module load sra/3.0.0

# Open vdb-config if needed
# Only once per session or new user
vdb-config --interactive

# Move to raw_data directory
cd /scratch/aubsxb005/project1_data/raw_data

#-----------------------
# Start Download Section
#-----------------------

# Each line downloads a sample and includes annotation for tracking
fastq-dump --split-files SRR4238380  # A. insolitus 
fastq-dump --split-files SRR4238381  # A. cybotes 
fastq-dump --split-files SRR4238382  # A. occultus 
fastq-dump --split-files SRR4238383  # A. lineatopus 
fastq-dump --split-files SRR4238384  # A. grahami 
fastq-dump --split-files SRR4238385  # A. valencienni 
fastq-dump --split-files SRR4238386  # A. sagrei 
fastq-dump --split-files SRR4238387  # A. cristatellus 
fastq-dump --split-files SRR4238388  # A. evermanni 
fastq-dump --split-files SRR4238389  # A. angusticeps 
fastq-dump --split-files SRR4238390  # A. chlorocyanus 
fastq-dump --split-files SRR31190375 # A. lemurinus 
fastq-dump --split-files SRR389085   # A. carolinensis 
fastq-dump --split-files SRR23804257 # A. apletophallus
fastq-dump --split-files SRR30501130 # A. tropidonotus 
fastq-dump --split-files SRR31253438 # A. planiceps
fastq-dump --split-files SRR30361216 # A. rodriguezii
fastq-dump --split-files SRR31186253 # A. sericeus
fastq-dump --split-files SRR21197535 # A. sagrei ordinatus (subspecies)
fastq-dump --split-files SRR21197563 # A. sagrei ordinatus (subspecies)
fastq-dump --split-files SRR25113491 # A. marmoratus 
fastq-dump --split-files SRR20731672 # A. baleatus 
fastq-dump --split-files SRR20731671 # A. semilineatus 
fastq-dump --split-files SRR8148313  # A. auratus 
fastq-dump --split-files SRR20731678 # A. frenatus 
fastq-dump --split-files DRR232286   # A. mestrei 
fastq-dump --split-files DRR232284   # A. allisoni 
fastq-dump --split-files DRR232285   # A. porcatus 
fastq-dump --split-files DRR232282   # A. garridoi 
fastq-dump --split-files SRR20731677 # A. isolepis 
fastq-dump --split-files DRR232281   # A. alutaceus 
fastq-dump --split-files DRR360734   # A. allogus 
fastq-dump --split-files SRR32222532 # A. distichus favillarum 
fastq-dump --split-files ERR3672845  # A. barbatus 
fastq-dump --split-files SRR24030180 # A. cuvieri 
fastq-dump --split-files DRR360735   # A. homolechis 
fastq-dump --split-files SRR24030181 # A. krugi 
fastq-dump --split-files SRR2651655  # A. loysiana
fastq-dump --split-files SRR26324123 # A. pulchellus 
fastq-dump --split-files SRR24030175 # A. stratulus 
fastq-dump --split-files SRR24030174 # A. gundlachi 
fastq-dump --split-files SRR24030177 # A. poncensis 
fastq-dump --split-files SRR13990308 # A. fuscoauratus 
fastq-dump --split-files SRR13990414 # A. brasiliensis 
fastq-dump --split-files SRR13990293 # A. ortonii 
fastq-dump --split-files SRR13990391 # A. chrysolepis
fastq-dump --split-files SRR7884574  # A. tandai 
fastq-dump --split-files SRR7884561  # A. dissimilis



echo "All downloads complete!"


