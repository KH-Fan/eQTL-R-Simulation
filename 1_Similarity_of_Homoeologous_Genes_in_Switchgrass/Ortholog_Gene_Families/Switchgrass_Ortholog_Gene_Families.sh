#!/bin/bash

# Exclude genes form other species and save all the switchgrass gene families in a file "sw_gene_families.txt"
awk '/Pavirv*/{print $1 ; for(i=1;i<=NF;++i)if($i~/Pavirv*/)print $i}' OrthologousGroups_WorkingDirectory.txt > sw_gene_families.txt

# Counting numbers of genes within a gene family 
awk '/Pavirv*/{print $1 ; j=0 ; for(i=1;i<=NF;++i) if($i~/Pavirv*/) j++ ; print j }' OrthologousGroups_WorkingDirectory.txt > sw_gene_family_counts.txt

# Retreave the gene names and the gene family names for the switchgrass gene families with exactly 2 genes in it
awk '/Pavirv*/{j=0 ; b=$1 ; for(i=1;i<=NF;++i) if($i~/Pavirv*/) {j++ ; b=b " " $i} ; if(j==2) print b }' OrthologousGroups_WorkingDirectory.txt > sw_2_genes_per_family.txt

# Retreave the gene names and save as an index file
awk '!/axpr/{print $2 , $3}' sw_2_genes_per_family.txt > gene_index.txt

# Retreave the gene names in the switchgrass transcriptome fasta file
awk '/>Pav/{print $0}' Pvirgatumv2.1.primaryTrs.fa > fasta_gene_name_index.txt

# Make the sequence in one line form the switchgrass transcriptome fasta file 
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'  Pvirgatumv2.1.primaryTrs.fa > Pvirgatumv2.1_seq_in_one_line.fa

# Get the sequence data for each transcript; One line per transcript
awk '!/^>/ {print $0}' Pvirgatumv2.1_seq_in_one_line.fa > fasta_sequences.txt

# Get the length of the transcripts in the switchgrass transcriptome
awk '{print length($0)}' fasta_sequences.txt > len_of_transc.txt
