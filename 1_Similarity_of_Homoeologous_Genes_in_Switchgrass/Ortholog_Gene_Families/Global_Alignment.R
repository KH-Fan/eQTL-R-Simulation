
rm(list=ls())
options(stringsAsFactors=FALSE)
# source("http://bioconductor.org/biocLite.R")
# biocLite("Biostrings")
# biocLite("DECIPHER")
setwd("~/Desktop/Draft_of_first_manuscript/Supplemental Matrerial/Scripts/1_Similarity_of _Paralogous_Genes_in_Switchgrass/Ortholog_Gene_Families/")

## Load Biostrings package for doing sequences alignment ## 
library(Biostrings)
# library(DECIPHER)

#########################################################
#                                                       #
# Read in gene names within a gene family as the index  #
#   Each row has two gene names, which are belonging    #
#   to the same gene family identified by OrthoFinder.  #
#                                                       #
#########################################################
gf_indx <- read.csv("gene_index.txt",header=FALSE,sep=" ")
# str(gf_indx)

#########################################################
#                                                       #
# Read in gene names in the fasta file as the index     #
#   Each row is the head line of each transcript        #
#   in the fasta file.                                  #
#                                                       #
#########################################################
seq_indx <- read.csv("fasta_gene_name_index.txt",header=FALSE, sep=" ")
# dim(seq_indx)

## Remove all non-needed information except the gene name ##
sindx <- unname(sapply(seq_indx$V1,substring,2))
# length(sindx)

#########################################################
#                                                       #
# Read in gene sequences in the fasta file              #
#   Each row is the sequence of each transcript         #
#   in the fasta file.                                  #
#   This will take a while depends on the size          #
#                                                       #
#########################################################
sw_seqs <- read.csv("fasta_sequences.txt",header=FALSE)
# dim(sw_seqs)

#########################################################
#                                                       #
# Looping all 8913 pairs of sequences and do            #
#   1. Global alignment of 2 sequences                  #
#   2. Counting numbers of mismatches in the aligment   #
#   3. Counting numbers of opening gaps                 #
#   4. return the sum of 2. & 3. for each alignment     # 
#                                                       #
#########################################################

## Settings for nucleotide subsitution matrix and penality for gaps
match <- 4
mismatch <- -4
gapPenality <- -2
gapExPenality <- 0

nucleotideSubMatrix <- nucleotideSubstitutionMatrix(match=match, 
                                                    mismatch=mismatch, 
                                                    baseOnly = FALSE)

## Function for transfer global alignment result into a matrix
split_str <- function(x){
  unlist(strsplit(as.character(x),split=""))
}

## Strings for saving results form the for loop
num_mismatches <- num_gap_opening <- align_length <- rep(-1,nrow(gf_indx))

## Looping all 8913 gene families with exactly 2 switchgrass genes 
for (i in 1:nrow(gf_indx)){

  # process tracking
  if(i%%100==0) {print(i)}
  
  # Find the current gene names of 2 genes in the family
  cur_gs <- gf_indx[i,]

  # Retreave 2 DNA sequences for the current gene family
  seq1 <- DNAString(sw_seqs[which(sindx%in%cur_gs[1]),])
  seq2 <- DNAString(sw_seqs[which(sindx%in%cur_gs[2]),])

  # Align these two sequences globally
  align <- pairwiseAlignment(seq1,seq2,substitutionMatrix=nucleotideSubMatrix,
                             gapOpening=gapPenality, gapExtension=gapExPenality)

  # Save alignment result into a 2 by N matrix
  align_seqs <- matrix(c(split_str(align@pattern),
                         split_str(align@subject)),nrow=2,byrow=TRUE)
  
  # Find the length of the alignment
  align_length[i] <- ncol(align_seqs)
  
  # Counting numbers of mismatches
  mis_matches <- apply(align_seqs,2,function(x){
                       ifelse(all(x[1]!="-",x[2]!="-") & x[1]!=x[2], 1, 0)})
  num_mismatches[i] <- sum(mis_matches)

  # Finding all individual gaps 
  indv_gaps <- apply(align_seqs,2,function(x){
                     ifelse(any(x[1]=="-",x[2]=="-") & x[1]!=x[2],1,0)})
  
  # Counting numbers of indel as the numbers of opening gaps
  indels <- diff(c(0,indv_gaps))
  num_gap_opening[i] <- length(which(indels==1))  
  
}

# rm(list=setdiff(ls(),c("num_mismatches","num_gap_opening","align_length")))
# save.image("Similarity.RData")
