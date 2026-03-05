# Get the gene list for each biological process from pantherdb GO enrichment.
# Separate numeric summaries for specific biological process(Immune response, innate immune response etc..):


# To check: 
tail -n +12 "Biologcal GO.txt" | grep -i -E "immune|innate|symbiont|oomyc|salt"

# Then: 
tail -n +8 "Biologcal GO.txt" | grep -i "immune response"
tail -n +8 "Biologcal GO.txt" | grep -i "innate immune"
tail -n +8 "Biologcal GO.txt" | grep -i "symbiont"
tail -n +8 "Biologcal GO.txt" | grep -i "oomycete"
tail -n +8 "Biologcal GO.txt" | grep -i "salt stress"

# And then converted fold enrichment for each and added that to this file and named as comparison.xlsx

# Extraction of candidate genes:

# Get gene list for each biological process from panther then run:

awk '{print $2}' innate_panther_raw.txt | grep "^Os" | sort | uniq > innate_panther_genes.txt
awk '{print $2}' symbiont_panther_raw.txt | grep "^Os" | sort | uniq > symbiont_panther_genes.txt
awk '{print $2}' oomycetes_panther_raw.txt | grep "^Os" | sort | uniq > oomycetes_panther_genes.txt
awk '{print $2}' salt_panther_raw.txt | grep "^Os" | sort | uniq > salt_panther_genes.txt
awk '{print $2}' negative_salt_panther_raw.txt | grep "^Os" | sort | uniq > negative_salt_panther_genes.txt

# Intersect with High-only( genes_HIGH_only.txt)

grep -xF -f genes_HIGH_only.txt innate_panther_genes.txt > innate_HIGH_candidates.txt
grep -xF -f genes_HIGH_only.txt symbiont_panther_genes.txt > symbiont_HIGH_candidates.txt
grep -xF -f genes_HIGH_only.txt oomycetes_panther_genes.txt > oomycetes_HIGH_candidates.txt
grep -xF -f genes_HIGH_only.txt salt_panther_genes.txt > salt_HIGH_candidates.txt
grep -xF -f genes_HIGH_only.txt negative_salt_panther_genes.txt > negative_salt_HIGH_candidates.txt

# Check unique:

cat *_HIGH_candidates.txt | sort | uniq | wc -l

# Then: 

cat *_HIGH_candidates.txt | sort | uniq > final_stress_HIGH_candidates.txt
wc -l final_stress_HIGH_candidates.txt

# Candidate genes: final_stress_HIGH_candidates.txt

# Get gene ids from TBtools.KEGG.enrichment.result.xls.final.xlsx and create and paste it in zeatin_genes.txt

# Intersect with HIGH-only:

grep -xF -f genes_HIGH_only.txt zeatin_genes.txt > zeatin_HIGH_candidates.txt
wc -l zeatin_HIGH_candidates.txt

# Then:

cat final_stress_HIGH_candidates.txt zeatin_HIGH_candidates.txt | sort | uniq > final_candidate_genes.txt
wc -l final_candidate_genes.txt
# And then:

cat final_stress_HIGH_candidates.txt zeatin_HIGH_candidates.txt | sort | uniq > final_candidate_genes.txt
wc -l final_candidate_genes.txt

# COMPARISION WITH OTHER GROUPS:
# Extract only true protein altering variants:
grep -E "stop_gained|frameshift_variant|splice_donor_variant|splice_acceptor_variant|start_lost" seeraga.snps.ann.vcf > strong_all.vcf

# Then: 
grep -o "gene:Os[0-9A-Za-z]*" strong_all.vcf | cut -d':' -f2 | sort | uniq > strong_all_genes.txt

# Intersect with the 69 genes:

grep -xF -f final_candidate_genes.txt strong_all_genes.txt > strong_candidate_genes.txt
wc -l strong_candidate_genes.txt


# Extract only variants in the strong_candidate_genes.txt
grep -F -f strong_candidate_genes.txt seeraga.snps.ann.vcf \
| grep -E "stop_gained|frameshift_variant|splice_donor_variant|splice_acceptor_variant|start_lost" \
> strong60_protein_altering.vcf

# Extract Clean Protein Change Table:

grep "ANN=" strong60_protein_altering.vcf \
| awk -F'ANN=' '{print $2}' \
| tr ',' '\n' \
| grep -E "stop_gained|frameshift_variant|splice_donor_variant|splice_acceptor_variant|start_lost" \
> strong60_annotations.txt

# Clean column:

awk -F'|' '{
    effect=$2;
    transcript=$4;
    protein_change=$11;
    split(transcript,a,"-");
    gene=a[1];
    print gene"\t"effect"\t"protein_change
}' strong60_annotations.txt | sort | uniq > strong60_protein_summary.txt

# Conversion: 

sed 's/CDS://g' strong60_protein_summary.txt \
| sed 's/t/g/' \
| sort | uniq > strong60_protein_summary_clean.txt

# Correct Extraction of Nipponbare Proteins for 60 Genes:

grep -A1 -F -f strong_candidate_genes.txt Oryza_sativa.IRGSP-1.0.pep.all.fa > nipponbare_strong60_proteins.fa

# Keep Only One Protein Per Gene:

awk '
/^>/ {
    split($0,a,"gene:");
    split(a[2],b," ");
    gene=b[1];
    if (!(gene in seen)) {
        seen[gene]=1;
        print $0;
        getline;
        print;
    }
    next;
}
' nipponbare_strong60_proteins.fa > nipponbare_strong60_unique.fa

# Calculate Protein Lengths:

awk '
/^>/ {
    if (seq) print gene"\t"length(seq);
    split($0,a,"gene:");
    split(a[2],b," ");
    gene=b[1];
    seq="";
    next
}
{seq=seq $0}
END {print gene"\t"length(seq)}
' nipponbare_strong60_unique.fa > nipponbare_protein_lengths.txt
# Extract Truncation Position:

awk '
{
    gene=$1;
    effect=$2;
    prot=$3;
    if (prot ~ /p\./) {
        gsub("p\\.","",prot);
        gsub("\\*","",prot);
        split(prot,a,"[A-Za-z]+");
        pos=a[1];
    } else {
        pos="NA";
    }
    print gene"\t"effect"\t"pos;
}
' strong60_protein_summary_clean.txt > predicted_truncation_positions.txt

# Merge Length + Truncation:

join -t $'\t' \
<(sort nipponbare_protein_lengths.txt) \
<(sort predicted_truncation_positions.txt) \
> protein_comparison_table.txt

# Correct Protein Extraction (Full Sequence)(this changes the previous step.):

# Create look up: 

awk '{print $1}' strong_candidate_genes.txt > gene_lookup.txt

# Extract full sequences;

awk '
BEGIN {
    while ((getline line < "gene_lookup.txt") > 0) {
        genes[line]=1
    }
    RS=">"; FS="\n"
}
NR>1 {
    header=$1
    seq=""
    for(i=2;i<=NF;i++) seq=seq $i
    split(header,a,"gene:")
    split(a[2],b," ")
    gene=b[1]
    if (gene in genes) {
        print ">"header
        print seq
    }
}
' Oryza_sativa.IRGSP-1.0.pep.all.fa > nipponbare_strong60_full.fa

# Then:

awk '
BEGIN {RS=">"; FS="\n"}
NR>1 {
    header=$1
    seq=""
    for(i=2;i<=NF;i++) seq=seq $i
    split(header,a,"gene:")
    split(a[2],b," ")
    gene=b[1]
    print gene "\t" length(seq)
}
' nipponbare_strong60_full.fa > nipponbare_protein_lengths.txt

# Extract truncation position cleanly:

awk '
$2=="stop_gained" {
    gene=$1
    effect=$2
    prot=$3

    if (match(prot, /[0-9]+/)) {
        pos=substr(prot, RSTART, RLENGTH)
    } else {
        pos="NA"
    }

    print gene "\t" effect "\t" pos
}
' strong60_protein_summary_clean.txt > stop_only_positions.txt

# Merge With Protein Lengths:

sort nipponbare_protein_lengths.txt > lengths_sorted.txt
sort stop_only_positions.txt > stop_sorted.txt

# Then:

join -t $'\t' lengths_sorted.txt stop_sorted.txt > merged_protein_data.txt

Calculate % Protein Loss:

awk '
{
    gene=$1
    full=$2
    effect=$3
    trunc=$4

    loss = full - trunc
    percent = (loss/full)*100

    printf "%s\t%d\t%s\t%d\t%.2f%%\n", gene, full, effect, trunc, percent
}
' merged_protein_data.txt > final_protein_impact_table.txt

# Clean the 55 Gene List:

cut -f1 final_protein_impact_table.txt | sort -u > strong55_genes.txt

# Get Indica Protein FASTA(Oryza_sativa_indica.pep.all.fa)

# Extract Proteins for 55 Genes From Indica:

awk '
BEGIN {
    while ((getline line < "strong55_genes.txt") > 0) {
        genes[line]=1
    }
    RS=">"; FS="\n"
}
NR>1 {
    header=$1
    seq=""
    for(i=2;i<=NF;i++) seq=seq $i
    split(header,a,"gene:")
    split(a[2],b," ")
    gene=b[1]
    if (gene in genes) {
        print ">"header
        print seq
    }
}
' Oryza_sativa_indica.pep.all.fa > indica_55_proteins.fa

# Create BLAST database from indica:

makeblastdb -in Oryza_sativa_indica.pep.all.fa -dbtype prot -out indica_db

BLAST your 36 Nipponbare proteins:

blastp \
  -query nipponbare_strong60_full.fa \
  -db indica_db \
  -out nip_vs_indica_blast.txt \
  -outfmt 6 \
  -max_target_seqs 1 \
  -evalue 1e-10


# Indica Protein Lengths:

cut -f2 nip_vs_indica_blast.txt | sort -u > indica_hits_ids.txt

# Then extract:

awk '
BEGIN {
    while ((getline line < "indica_hits_ids.txt") > 0) {
        ids[line]=1
    }
    RS=">"; FS="\n"
}
NR>1 {
    header=$1
    seq=""
    for(i=2;i<=NF;i++) seq=seq $i
    split(header,a," ")
    id=a[1]
    if (id in ids) {
        print ">"header
        print seq
    }
}
' Oryza_sativa_indica.pep.all.fa > indica_matched_proteins.fa

# Compute Indica Lengths:

awk '
BEGIN {RS=">"; FS="\n"}
NR>1 {
    header=$1
    seq=""
    for(i=2;i<=NF;i++) seq=seq $i
    split(header,a," ")
    id=a[1]
    print id "\t" length(seq)
}
' indica_matched_proteins.fa > indica_lengths.txt

# Attach Indica Length To Each Nipponbare Gene:

awk '
{
    split($1,a,"-")
    transcript=a[1]
    gsub("t","g",transcript)
    gene=transcript
    print gene "\t" $2
}
' nip_vs_indica_blast.txt > gene_to_indicaID.txt

# Merge With Indica Lengths:

sort gene_to_indicaID.txt > gene_ind_sorted.txt
sort indica_lengths.txt > indica_len_sorted.txt

# First rearrange gene_ind_sorted so ID is first:

awk '{print $2 "\t" $1}' gene_ind_sorted.txt > flipped.txt

# Sort flipped:

sort flipped.txt > flipped_sorted.txt

# Join:

join -t $'\t' flipped_sorted.txt indica_len_sorted.txt > indica_gene_lengths.txt

# flip back:

awk '{print $2 "\t" $3}' indica_gene_lengths.txt > gene_indica_length_table.txt

# Merge All Data Together:

sort final_protein_impact_table.txt > final_sorted.txt
sort gene_indica_length_table.txt > indica_gene_sorted.txt

# Join:

join -t $'\t' final_sorted.txt indica_gene_sorted.txt > final_comparative_table.txt

# Convert it to excel:

awk '
{
    gene=$1
    nip_len=$2
    trunc=$4
    indica=$6

    # Keep longest Nip length
    if (!(gene in max_nip) || nip_len > max_nip[gene]) {
        max_nip[gene]=nip_len
    }

    # Keep earliest truncation (smallest number)
    if (!(gene in min_trunc) || trunc < min_trunc[gene]) {
        min_trunc[gene]=trunc
    }

    # Keep indica length (just store one)
    indica_len[gene]=indica
}
END {
    for (g in max_nip) {
        full=max_nip[g]
        trunc=min_trunc[g]
        loss=full-trunc
        percent=(loss/full)*100
        printf "%s\t%d\t%d\t%d\t%.2f%%\n", g, full, indica_len[g], trunc, percent
    }
}
' final_comparative_table.txt > final_clean_table.txt


##Notes:  The biological processes selected above are just for temporary comparision. The most enriched pathways will be analysed.
