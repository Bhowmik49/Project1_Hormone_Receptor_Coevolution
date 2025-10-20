# 🦎 Project1: Hormone–Receptor Coevolution in *Anolis* Lizards

This repository contains the computational workflow, scripts, and dataset summaries used to process whole-genome and RNA-seq data for the study of **insulin/insulin-like signaling (IIS) pathway evolution** across *Anolis* lizards.  
The pipeline converts raw reads into **N-masked consensus sequences**, then builds **per-gene alignments** for comparative and evolutionary analysis.

---

## ⚙️ Overview

All analyses were performed on the **Alabama Supercomputer (ASC)** using a modular Bash-based workflow.  
Each step automatically records modules, runtime logs, and outputs for full reproducibility.

The pipeline includes **12 major stages**:

1. **Download data** (SRA Toolkit)
2. **Retrieve reference genome** (NCBI)
3. **Quality control & trimming** (FastQC, MultiQC, Trimmomatic)
4. **Reference indexing** (Bowtie2, HISAT2)
5. **Read mapping** (Bowtie2/HISAT2 + SAMtools)
6. **Extract IIS gene coordinates** (awk, GTF)
7. **Subset BAMs to target genes** (SAMtools, BEDTools)
8. **Variant calling & consensus FASTA** (bcftools, vcfutils.pl)
9. **Split BEDs by gene & feature**
10. **Extract & concatenate CDS FASTAs** (BEDTools)
11. **Multiple sequence alignment** (MAFFT; codon-aware via PAL2NAL)
12. **Missing-data matrix generation** (bash, awk, python)

---

## 🧬 Repository Structure
```
Project1_Hormone_Receptor_Coevolution/
├── scripts/
│   ├── 01_Download/
│   ├── 01b_reference/
│   ├── 02_trimming/
│   ├── 03_qc_testing/
│   ├── 04_ref_indexing/
│   ├── 05_mapping/
│   ├── 06_extract_gene_coordinates_bed_files/
│   ├── 07_subset_bam/
│   ├── 08_call_variants_consens/
│   ├── 09_split_beds/
│   ├── 10_extracted_consensus_fastas_by_gene/
│   ├── 11_alignment/
│   └── 12_missing_data_matrix/
├── accessions.txt
├── Data_coverage.xls
├── Coverage_table.docx
├── Overview_Available_Data (Working_Dataset2).csv
└── README.md

```
---


## 🧰 Software Environment

| Tool | Purpose |
|------|----------|
| **FastQC / MultiQC** | Quality control |
| **Trimmomatic** | Adapter & quality trimming |
| **Bowtie2 / HISAT2** | WGS and RNA-seq mapping |
| **SAMtools** | Sorting, indexing, depth masking |
| **bcftools / vcfutils.pl** | Variant calling & consensus generation |
| **BEDTools** | Gene-level extraction |
| **MAFFT / PAL2NAL** | Alignment (nucleotide & codon-aware) |
| **ASC Modules** | Managed HPC environment |

---

## 📊 Dataset Summary

The dataset spans **37 *Anolis* species**, including both reference genomes and new RNA-seq libraries.  
Sequencing sources include:

- **NCBI SRA:** 15 WGS datasets + 2 reference genomes  
- **Schwartz Lab:** 13 RNA-seq libraries  
- **Tollis et al. 2018:** 3 transcriptomes  

Data completeness and sequence quality are summarized in:

- `Coverage_table.docx`
- `Data_coverage.xls`
- `Overview_Available_Data (Working_Dataset2).csv`

| Code | Description |
|------|--------------|
| **C** | Complete (≥95 % CDS coverage, < 5 % Ns) |
| **P** | Partial (70–94 % coverage, 5–30 % Ns) |
| **I** | Incomplete (30–69 % coverage, > 30 % Ns) |
| **M** | Missing (<30 % coverage or no data) |

---

## 🧩 Output Summary

| Folder | Output | Description |
|---------|---------|-------------|
| `variants/consensus_fastas_5/` | `.fasta` | N-masked per-sample consensus sequences |
| `variants/cds_fasta_5/per_gene/` | `.fasta` | Multi-species CDS FASTAs per gene |
| `alignments/final_aligned_CDS/` | `.fasta` | MAFFT-aligned CDS sequences |
| `missing_data/` | `.tsv` | Completeness matrix |
| `qc_reports/` | `.html` | MultiQC summaries |

---

## ✅ Validation

- *A. sagrei* reference sequences recovered at **> 99 % identity**  
- Coverage-based masking ensures only high-confidence bases  
- Completeness matrix quantifies recovery for each gene and species  
- Reproducibility confirmed via logged ASC environment and version control  

---

## 🔗 Repository Link

**GitHub:**  
[https://github.com/Bhowmik49/Project1_Hormone_Receptor_Coevolution](https://github.com/Bhowmik49/Project1_Hormone_Receptor_Coevolution)

---

## 📘 Citation

If you use this workflow or data, please cite:

> **Bhowmik S. (2025).** *Hormone–Receptor Coevolution in the Insulin/IGF Signaling Pathway of Anolis Lizards.* Auburn University Dissertation Project.  
> Tools referenced: Li et al. 2009 (SAMtools), Langmead & Salzberg 2012 (Bowtie2), Kim et al. 2015 (HISAT2), Katoh & Standley 2013 (MAFFT).

---

🧬 **Maintainer:** Sagar Bhowmik (Ocean)  
📧 szb0232@auburn.edu  | Schwartz Lab, Department of Biological Sciences, Auburn University, Auburn, AL 36849. 
