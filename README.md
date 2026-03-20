Computational pipeline developed for genomic variant discovery,
functional annotation, and candidate gene prioritization in traditional rice varieties.

# Functional Annotation and Comparison of Traditional Rice Varieties

Bioinformatics pipeline for genomic variant discovery, functional annotation, and candidate gene prioritization in traditional rice varieties.

This repository contains the computational workflow used to analyze whole-genome sequencing data from traditional rice varieties and identify potentially functional genetic variants associated with stress response and metabolic pathways.

The pipeline is demonstrated using **Seeraga Samba** sequencing data.  
The same workflow was applied to additional rice varieties such as Vellai Kowni and Madu Muzhungi for comparative genomic analysis which are the core parts of this study.

---

## Pipeline Overview

The workflow consists of the following major steps:

1. **Quality Control and Trimming**
   - fastp

2. **Read Alignment**
   - Mapping to the *Oryza sativa* Nipponbare reference genome using BWA

3. **BAM Processing**
   - SAMtools sorting and indexing  
   - Picard duplicate marking

4. **Variant Calling**
   - SNP discovery using GATK HaplotypeCaller

5. **Variant Annotation**
   - Functional impact prediction using SnpEff

6. **Functional Annotation**
   - Orthology and pathway mapping using EggNOG-mapper

7. **Pathway and GO Enrichment**
   - KEGG enrichment using TBtools  
   - GO enrichment using PantherDB

8. **Candidate Gene Identification**
   - Extraction of genes with **HIGH impact variants**
   - Filtering for **loss-of-function mutations** including:
     - stop_gained
     - frameshift_variant
     - splice_donor_variant
     - splice_acceptor_variant
     - start_lost

9. **Protein Impact Analysis**
   - Prediction of truncation positions
   - Calculation of percent protein loss

10. **Comparative Analysis**
   - Comparison of affected proteins with *Oryza sativa indica* reference proteins using BLAST

---

## Data Sources

Raw sequencing data  
NCBI Sequence Read Archive (SRA)

Example dataset used in this repository:

SRR3625307  
https://www.ncbi.nlm.nih.gov/sra/SRR3625307

Reference genome:

*Oryza sativa* Nipponbare IRGSP-1.0  
NCBI Assembly: GCF_001433935.1

Genome FASTA and annotation files were downloaded from the NCBI genome database.

---

## Repository Structure

pipeline/ variant_pipeline.sh
scripts/ GO_enrichment_plot.R KEGG_plot.R
figures/ GO_barplot.png KEGG_bubble_plot.png
results/ (intermediate and processed outputs)

---

## Notes

Large sequencing files (FASTQ, BAM, VCF) are not included in this repository.

All raw data can be downloaded directly from the NCBI SRA using the accession numbers provided above.

---

## Citation

If you use this pipeline or repository in your work, please cite the associated publication (in preparation).
