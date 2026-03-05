##Functional Annotation and Comparison of Traditional Rice Varieties
Bioinformatics pipeline for genomic variant discovery, functional annotation, and candidate gene prioritization in traditional rice varieties.
This repository contains the computational workflow used to analyze whole-genome sequencing data from traditional rice varieties and identify potentially functional genetic variants associated with stress response and metabolic pathways.
The pipeline is demonstrated using Seeraga Samba sequencing data.
The same workflow was applied to additional rice varieties for comparative genomic analysis.

##Pipeline Overview
The workflow consists of the following major steps:


#Quality Control and Trimming

fastp

#Read Alignment

Mapping to the Oryza sativa Nipponbare reference genome using BWA


#BAM Processing

SAMtools sorting and indexing
Picard duplicate marking


#Variant Calling

SNP discovery using GATK HaplotypeCaller

#Variant Annotation

Functional impact prediction using SnpEff


#Functional Annotation

Orthology and pathway mapping using EggNOG-mapper


#Pathway and GO Enrichment

KEGG enrichment using TBtools
GO enrichment using PantherDB


#Candidate Gene Identification

Extraction of genes with HIGH impact variants
Filtering for loss-of-function mutations including:

stop_gained
frameshift_variant
splice_donor_variant
splice_acceptor_variant
start_lost



#Protein Impact Analysis

Prediction of truncation positions
Calculation of percent protein loss



#Comparative Analysis



Comparison of affected proteins with Oryza sativa indica reference proteins using BLAST.


##Data Sources
Raw sequencing data
NCBI Sequence Read Archive (SRA)
Example dataset used in this repository:
SRR3625307
https://www.ncbi.nlm.nih.gov/sra/SRR3625307
Reference genome:
Oryza sativa Nipponbare IRGSP-1.0
NCBI Assembly: GCF_001433935.1
Genome FASTA and annotation files were downloaded from the NCBI genome database.

##Repository Structure
pipeline/
variant_pipeline.sh
scripts/
GO_enrichment_plot.R
KEGG_plot.R
figures/
GO_barplot.png
KEGG_bubble_plot.png
results/
(intermediate and processed outputs)

##Notes
Large sequencing files (FASTQ, BAM, VCF) are not included in this repository.
All raw data can be downloaded directly from the NCBI SRA using the accession numbers provided above.

##Citation
If you use this pipeline or repository in your work, please cite the associated publication (in preparation).

