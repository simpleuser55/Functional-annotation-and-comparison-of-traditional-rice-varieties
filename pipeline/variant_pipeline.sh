# Ubuntu
 
prefetch SRR3625307
fasterq-dump SRR3625307 -O 

# Result files:

# SRR3625307_1.fastq
# SRR3625307_2.fastq

# FASTP (trimming + QC):

fastp \
  -i SRR3625307_1.fastq \
  -I SRR3625307_2.fastq \
  -o SRR3625307_1.trim.fastq \
  -O SRR3625307_2.trim.fastq \
  -h ../results/fastp_SRR3625307.html \
  -j ../results/fastp_SRR3625307.json \
  -w 6

# NIPPONBARE GENOME AND ANNOTATION:
# Genome (FASTA):
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_genomic.fna.gz

# Annotation (GFF3):

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_genomic.gff.gz

# UNZIPPING:

gunzip GCF_001433935.1_IRGSP-1.0_genomic.fna.gz(Fasta)
gunzip GCF_001433935.1_IRGSP-1.0_genomic.gff.gz(GFF3)

# Indexing the genome :
bwa index nipponbare.fa

# Mapping:
bwa mem -t 8 nipponbare.fa \
  ../data/SRR3625307_1.trim.fastq \
  ../data/SRR3625307_2.trim.fastq \
  | samtools view -bS - > ../results/seeraga_vs_nipponbare.bam

# Then:
samtools sort -@6 -o ../results/seeraga_vs_nipponbare.sorted.bam ../results/seeraga_vs_nipponbare.bam
samtools index ../results/seeraga_vs_nipponbare.sorted.bam
samtools flagstat ../results/seeraga_vs_nipponbare.sorted.bam


# Add sample metadata. Does not change reads

java -jar picard.jar AddOrReplaceReadGroups \
  I=seeraga_vs_nipponbare.sorted.bam \
  O=seeraga.RG.bam \
  RGID=seeraga1 \
  RGLB=lib1 \
  RGPL=ILLUMINA \
  RGPU=unit1 \
  RGSM=seeraga

# Index the Read-Group BAM

samtools index seeraga.RG.bam

# Resulted files:
# seeraga.RG.bam
# seeraga.RG.bam.bai

# Picard: Mark Duplicates:

java -jar picard.jar MarkDuplicates \
  I=seeraga.RG.bam \
  O=seeraga.dedup.bam \
  M=seeraga.dup_metrics.txt

# Index the Deduplicated BAM

samtools index seeraga.dedup.bam

# Results files:
# seeraga.dedup.bam
# seeraga.dedup.bam.bai

# Create .dict file:
java -jar /mnt/e/rice_poc/picard.jar CreateSequenceDictionary \
  R=nipponbare.fa \
  O=nipponbare.dict

# GATK: Variant Calling (HaplotypeCaller):
/mnt/e/rice_poc/gatk-4.6.1.0/gatk HaplotypeCaller \
  -R /mnt/e/rice_poc/results/nipponbare.fa \
  -I /mnt/e/rice_poc/results/seeraga.dedup.bam \
  -O /mnt/e/rice_poc/results/seeraga.raw.vcf.gz

# Resulted file: seeraga.raw.vcf.gz
# FILTER VARIANTS (keep only high-confidence):
# Extract SNPs:
/mnt/e/rice_poc/gatk-4.6.1.0/gatk SelectVariants \
  -R nipponbare.fa \
  -V seeraga.raw.vcf.gz \
  --select-type-to-include SNP \
  -O seeraga.snps.vcf.gz
# Apply basic hard filters (safe defaults):
/mnt/e/rice_poc/gatk-4.6.1.0/gatk VariantFiltration \
  -R nipponbare.fa \
  -V seeraga.snps.vcf.gz \
  -O seeraga.snps.filtered.vcf.gz \
  --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
  --filter-name "BASIC_SNP_FILTER"

# Resulted file(RF): seeraga.snps.filtered.vcf.gz
# Configure SnpEff for Nipponbare (IRGSP-1.0):
nano snpEff/snpEff.config
# Add this at the very bottom:
IRGSP1.0.genome : Oryza_sativa_Nipponbare
# Build SnpEff database: 
java -jar snpEff/snpEff.jar build \
  -gff3 \
  -noCheckCds \
  -noCheckProtein \
  IRGSP1.0

# After build succeeds, annotate variants:
java -jar snpEff/snpEff.jar IRGSP1.0 \
  /mnt/e/rice_poc/results/seeraga.snps.filtered.vcf.gz \
  > /mnt/e/rice_poc/results/seeraga.snps.ann.vcf

# Get the summary:
grep "ANN=" /mnt/e/rice_poc/results/seeraga.snps.ann.vcf \
| grep -o "HIGH\|MODERATE\|LOW" \
| sort | uniq -c
 
# Separate HIGH and MODERATE genes:
grep -E "HIGH|MODERATE" seeraga.snps.ann.vcf \
| grep -o "Os[0-9]\+g[0-9]\+" \
| sort | uniq > genes_HIGH_MODERATE.txt

# Get protein FASTA for those genes:
# We already have : Oryza_sativa.IRGSP-1.0.pep.all.fa
grep -F -f genes_HIGH_MODERATE.txt Oryza_sativa.IRGSP-1.0.pep.all.fa -A 1 \
| grep -v "^--$" \
> seeraga_HIGH_MODERATE_proteins.fa

# genes: 28592
# Protein seq: 32707

Collapse to 1 protein per gene:
awk '
/^>/ {
    match($0, /Os[0-9]+g[0-9]+/, gene)
    if (!(gene[0] in seen)) {
        seen[gene[0]] = 1
        keep = 1
    } else {
        keep = 0
    }
}
keep { print }
' seeraga_HIGH_MODERATE_proteins.fa \
> seeraga_HIGH_MODERATE_proteins_unique.fa

New protein count: 27251

# Run EggNOG:
emapper.py \
-i seeraga_HIGH_MODERATE_proteins_unique.fa \
--itype proteins \
-m diamond \
--data_dir /mnt/e/eggnog_data \
--output seeraga_eggnog \
--cpu 4

# Resulting file: seeraga_eggnog.emapper.annotations
# Extract KEGG KO IDs:

awk -F "\t" '$12 != "-" {print $1 "\t" $12}' seeraga_eggnog.emapper.annotations > gene_KO_raw.txt

# Split multiple KO IDs:

awk -F "\t" '{split($2,a,","); for(i in a) print $1 "\t" a[i]}' gene_KO_raw.txt > gene_KO_for_TBtools.txt

# Prepare Gene List for Background:

awk -F "\t" '
{
    if (a[$1]=="") a[$1]=$2;
    else a[$1]=a[$1]","$2
}
END {
    for (i in a) print i "\t" a[i]
}
' gene_KO_for_TBtools.txt > TBtools_KO_background.txt

# Reformats:

awk -F "\t" '
{
    match($1, /Os[0-9]+t[0-9]+/, id)
    gene=id[0]
    sub("t","g",gene)
    if (a[gene]=="") a[gene]=$2;
    else a[gene]=a[gene]","$2
}
END {
    for (i in a) print i "\t" a[i]
}
' gene_KO_for_TBtools.txt > TBtools_KO_genelevel.txt

# Then:

awk -F "\t" '
NR>1 {
    gsub("ko:","",$2);
    print $1 "\t" $2
}
' TBtools_KO_genelevel.txt > TBtools_KO_genelevel_clean.txt


# Clean to Gene-Level IDs:
awk -F "\t" '
NR>1 {
    gsub("ko:","",$2);
    print $1 "\t" $2
}
' TBtools_KO_genelevel.txt > TBtools_KO_genelevel_clean.txt

(background)
#Generate HIGH-only gene list:
grep "HIGH" seeraga.snps.ann.vcf \
| grep -o "Os[0-9]\+g[0-9]\+" \
| sort | uniq > genes_HIGH_only.txt
#High only:4194
