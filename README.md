Bioinformatics pipeline for variant discovery and functional analysis of rice genomes.

# Rice Variant Functional Analysis Pipeline

This repository contains the computational workflow used for genomic variant discovery and functional annotation in traditional rice varieties.

The pipeline is demonstrated using **Seeraga Samba** whole-genome sequencing data.  
The same workflow was applied to additional rice varieties for comparative genomic analysis.

---

## Pipeline Overview

1. Quality control and trimming using fastp
2. Alignment to Nipponbare reference genome using BWA
3. BAM processing using SAMtools and Picard
4. Variant calling using GATK
5. Variant annotation using SnpEff
6. Functional annotation using EggNOG-mapper
7. KEGG enrichment analysis using TBtools
8. GO enrichment analysis using PantherDB
9. Candidate gene extraction
10. Protein truncation analysis
11. Comparative protein analysis with indica rice

---

## Data Sources

Raw sequencing reads:
https://www.ncbi.nlm.nih.gov/sra/SRR3625307

Reference genome:
Oryza sativa Nipponbare IRGSP-1.0  
NCBI Assembly: GCF_001433935.1

Annotation:
IRGSP-1.0 GFF3

---

## Tools Used

fastp  
BWA  
SAMtools  
Picard  
GATK  
SnpEff  
EggNOG-mapper  
TBtools  
R (ggplot2, dplyr)

---

## Repository Structure

pipeline/ variant_pipeline.sh
scripts/ GO_enrichment_plot.R KEGG_plot.R
figures/ GO_barplot.png KEGG_bubble_plot.png

---

## Notes

Large sequencing files (FASTQ, BAM, VCF) are not included in this repository.  
They can be downloaded from NCBI SRA using the accession number provided above.

The pipeline is demonstrated using **Seeraga Samba** whole-genome sequencing data.  
The same workflow was applied to two more traditional rice varieties for comparative genomic analysis.

## Citation

If you use this pipeline in your research, please cite the associated publication (in preparation).
