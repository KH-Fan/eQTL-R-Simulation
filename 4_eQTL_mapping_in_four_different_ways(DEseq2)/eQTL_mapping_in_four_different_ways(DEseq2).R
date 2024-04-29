####################################################################################
#
# The R script has 4 simulation functions:
# 
# The first function will
#   randomly generate gene expression levels (# of short reads)
#   for corresponding genotype from nagative binomial distribution
# 
# The second function will
#   return number of transcripts form either from A genome, 
#   from B genome, or from the unknown origion genome 
#
# The third function will
#   combines two previous functions (genotype_explevel & GSE_sim),
#   generate random data sets for assign unknown origion counts into A & B genomes
#   Note:
#     This function can only take care of 2 genomes so far
#
# The forth function will
#   iteratively call Gen_GSE_Data to create the matrixes 
#   ready for eQTL analysis by using Deseq2 package
#
# Finally, the last part of the script is applying different combination of 
# parameter settings into previous four funcitons to determine the lose of power 
# and increase of false postive corresponding to 3 different way 
# separating unknown origin counts (UOC). 
#
#
#     University of Georgia
#     Insititute of Bioinformatics
# 
#                             Kang-Hsien Fan
#                                 02/18/2016
####################################################################################

###################################
###   Simulation function # 1   ###
###################################

####################################################################################
#
# Function name: genotype_expression_level (genotype_explevel)
#
# Description: 
#   Randomly generate gene expression levels (# of short reads)
#   according to the genotype from negative binomial distribution
#   
# Parameters (default value):
# (1) mu_qq:            mean expression level for genotype qq (100)
# (2) fd:               flod changs between different genotypes (1)
#                           qq: mu_qq; Qq: fd*mu_qq; QQ: (2*fd-1)*mu_qq
# (3) phi:              overdispersion parameter (1)
# (4) ss:               sample size for the genotype qq (1)
#                           if there is only one element assigned, 
#                           then the function will amplify it with 
#                           the assigned genotypic ratio
# (5) genotypic_ratio:  the sample size ratio for the genotype (1:2:1)
#                           assum the hardy-weinberg equilibrium stands
# 
# Output:
#   A list with random negative binomial variables with corresponding genotype
#
# Dependance:
#   rnegbin function in MASS package
#
####################################################################################

library(MASS)
# rnegbin function in MASS package is needed for 
# generate random negative binomial variables.

genotype_explevel <- function(mu_qq=100,fd=1,phi=1,ss=1,genotypic_ratio=c(1,2,1)){ 
  
  if(length(ss)==1){ ss <- ss*genotypic_ratio }
  # check if there is only one element assign in ss,
  # if so, amplify it with the genotypic ratio.
  
  if(length(ss)!=3){
    break("There should be 3 sample saizes for each genotype")
  }
  # check if the sample sizes is properly assigned
  
  # corresponding mean value for each genotype
  mu_X <- c(mu_qq,mu_qq*fd,mu_qq*(2*fd-1))
  # generate ss random numbers of reads mapped to gene with genotype X
  
  n_AA <- rnegbin(ss[1],mu_X[1],phi)  
  # expression level for genotype AA
  n_AB <- rnegbin(ss[2],mu_X[2],phi)
  # expression level for genotype AB
  n_BB <- rnegbin(ss[3],mu_X[3],phi)
  # expression level for genotype BB
  
  return(list(AA=n_AA,AB=n_AB,BB=n_BB))
  # Return all the expression levels 
}

###################################
###   Simulation function # 2   ###
###################################

####################################################################################
#
# Function name:  genome specific counts simulator (GSE_sim)
#
# Description: 
#   simulate the number of RNA-seq counts form either from A genome, 
#   from B genome, or from the unknown origion genome 
#   
# Parameters (default value):
# (1) len_tans:       length of the transcript in bps (10kb)
# (2) n_diff_sites:   number of different site(s) between all the genomes (1)
# (3) n_reads:        number of reads mapped to each transcripts. (100)
#                         if there is only one element assigned, 
#                         then the function will assume equal read number 
#                         for all genome 
# (4) len_reads:      length of the sequencing short reads for each transcripts. (200) 
#                         if there is only one element assigned, 
#                         then the function will assume equal read length 
#                         for all genome
# (5) copy_genome:    copy of the genome (2)
#                         default is 2 for allotetraploid (A and B genome)
#
# Output:
#   A list of two components: 1) the random picked SNP(s) sites; 
#   and 2) the table of counts with information of genome origion 
#
# Dependance:
#   None
#
####################################################################################

GSE_sim <- function(len_trans=10000,n_diff_sites=1,n_reads=100,
                    len_reads=200,copy_genome=2){
  
  if(n_diff_sites>len_trans){
    break("Number of different sites between A and B genome is larger than the transcript itslef.")
  }
  # Check if the different sites are larger than transcript
  
  diff_sites <- sample(1:len_trans,n_diff_sites) 
  # randomly pick up the locations of the different sites on the A & B genome
  
  if(length(n_reads)==1){n_reads <- rep(n_reads,copy_genome)}
  if(length(len_reads)==1){len_reads <- rep(len_reads,copy_genome)}
  # if there is only one element in number of reads and/or the length of the read,
  # then we will assume they are all the same for all copy of the genomes
  
  if(length(n_reads)!=copy_genome){
    break("Please assign the number of reads for each transcript if they are different")
  }
  # check if number of reads for each genome copy are assigned properly 
  
  if(length(len_reads)!=copy_genome){
    break("Please assign short reads length for each transcript if they are different")
  }
  # check if the read length for each genome copy are assigned properly 
  
  reads_genome <- matrix(NA,nrow=copy_genome,ncol=2) 
  # an empty matrix for storing transcript counts for each genome copy
  
  # for each genome
  for(i in 1:copy_genome){
    
    if(n_reads[i]==0){
      
      reads_genome[i,] <- c(0,0)
      # if there is no reads mapped to the transcript
      
    }else{
      
      left_pos_reads <- sample(1:(len_trans-len_reads[i]),n_reads[i],replace=TRUE)
      # left end position of the reads
      
      reads_end_pos <- matrix(c(left_pos_reads,(left_pos_reads+len_reads[i])),
                              nrow=n_reads[i],byrow=FALSE)
      # a matrix of both ends of the reads 
      
      cover_diff <- do.call(cbind,lapply(diff_sites,function(x){
        ind_pos <- sign(reads_end_pos-x)
        # use sign to indicate if the different site is covered by a read
        return(ind_pos[,1]-ind_pos[,2]!=0)
        # find the reads that cover the different sites on the genome
      }))
      # create a TRUE/FASLE matrix to indicate if a reads covered each diffferent site
      
      reads_genome[i,] <- table(factor(apply(cover_diff,1,any),levels=c("FALSE","TRUE")))
      # return a count table of weather the reads can be determine form which trancsript
    }
  }
  
  reads_genome <- cbind(reads_genome,apply(reads_genome,1,sum))
  colnames(reads_genome) <- c("Reads_Uncovered","Reads_Covered","N_reads")
  rownames(reads_genome) <- LETTERS[1:copy_genome]
  
  return(list(diff_sites=diff_sites,reads_genome=reads_genome))  
  # return a count table of where the reads come from 
}

###################################
###   Simulation function # 3   ###
###################################

####################################################################################
#
# Function name:  Generate genome specific expression data (Gen_GSE_Data)
#
# Description: 
#   This function combines two previous functions (genotype_explevel & GSE_sim),
#   generate random data sets for assign unknown origion counts into A & B genomes
#
# Note:
#   This function can only take care of 2 genomes so far
#   
# Parameters (default value):
#     All the parameters have the same definitions as those two functions.
# (1) mu_A_qq:          mean expression level in genome A for genotype qq (100)
# (2) mu_B_qq:          mean expression level in genome B for genotype qq (100)
# (3) fd_A:             flod changs in genome A (1)
# (4) fd_B:             flod changs in genome B (3)
# (5) phi:              overdispersion parameter (1)
# (6) ss:               sample size for the genotype qq (25)
# (8) len_tans:         length of the transcript in bps (10kb)
# (9) n_diff_sites:     number of different site(s) between all the genomes (1)
# (10) len_reads:       length of the sequencing short reads for each transcripts. (200) 
# (11) genotypic_ratio: the sample size ratio for the genotype (1:2:1)
#
# Output:
#   A dataframe contain the genotype of genome A & B, as well as the true (simulated)
#   RNA-seq counts, genome specific counts for both genome, and the unknown origion 
#   gene counts, as well as 3 different ways to separate unknown-origion counts 
#   evenly, proportionally, and sums
#
# Dependance:
#   genotype_explevel & GSE_sim funcitons
#
####################################################################################

Gen_GSE_Data <- function(mu_A_qq=100,mu_B_qq=100,fd_A=1,fd_B=3,phi=1,ss=25,
                         len_trans=10000,n_diff_sites=1,len_reads=200,
                         genotypic_ratio=c(1,2,1)){
  
  ss_gen_A <- genotype_explevel(mu_A_qq,fd_A,phi,ss,genotypic_ratio)
  ss_gen_B <- genotype_explevel(mu_B_qq,fd_B,phi,ss,genotypic_ratio)
  # randomly generate counts for two genome transcripts using "genotype_explevel" function
  
  d_gen_A <- data.frame(Genotype_A=c(rep("qq",length(ss_gen_A[[1]])),
                                     rep("Qq",length(ss_gen_A[[2]])),
                                     rep("QQ",length(ss_gen_A[[3]]))),
                        count_true_A=c(ss_gen_A[[1]],ss_gen_A[[2]],ss_gen_A[[3]]))
  
  d_gen_B <- data.frame(Genotype_B=c(rep("qq",length(ss_gen_B[[1]])),
                                     rep("Qq",length(ss_gen_B[[2]])),
                                     rep("QQ",length(ss_gen_B[[3]]))),
                        count_true_B=c(ss_gen_B[[1]],ss_gen_B[[2]],ss_gen_B[[3]]))
  # expand the output of the "genotype_explevel" function 
  
  nr <- nrow(d_gen_A)
  
  ran_d_gen_A <- d_gen_A[sample(1:nr),]
  ran_d_gen_B <- d_gen_B[sample(1:nr),]
  # randomly shuffle the data 
  
  ran_d <- cbind(ran_d_gen_A,ran_d_gen_B)
  # combine two transcripts counts into one data frame
  
  counts <- ran_d[,c("count_true_A","count_true_B")]
  # extract transcript counts for "GSE_sim" function
  
  reads_genomes <- t(apply(counts,1,function(x){
    temp <- GSE_sim(n_reads=x,len_trans,n_diff_sites,len_reads)$reads_genome
    return(c(temp[,2],sum(temp[,1])))}
  ))
  # simulate reads could and coouldn't be identified from "GSE_sim" function
  
  d <- cbind(ran_d,reads_genomes)
  names(d)[5:7] <- c("GSE_A","GSE_B","Un_GSE")
  # combine all the information into a single data frame
  
  d$even_A <- round(d$GSE_A + d$Un_GSE/2,0)
  d$even_B <- round(d$GSE_B + d$Un_GSE/2,0)
  # calculate the genome expression by evenly separating unknown-origin counts (UOC)
  
  d$prop_for_A <- ifelse(d$GSE_A==0 & d$GSE_B==0, 0.5,
                         ifelse(d$GSE_A!=0 & d$GSE_B==0, d$GSE_A/(d$GSE_A+1),
                                ifelse(d$GSE_A==0 & d$GSE_B!=0, 1/(1+d$GSE_B),
                                       d$GSE_A/(d$GSE_A+d$GSE_B))))
  # calculate the GSE_A/GSE_B ratio
  
  d$prop_A <- round(d$GSE_A + (d$prop_for_A * d$Un_GSE), 0)
  d$prop_B <- round(d$GSE_B + ((1-d$prop_for_A) * d$Un_GSE), 0)
  # calculate the genome expression by proportionally separating unknown-origin counts 
  # according to the GSE_A/GSE_B raito  
  
  d$sum_count <- d$GSE_A + d$GSE_B + d$Un_GSE
  # summing all the homoe-allele gene counts together
  
  return(d)
}

###################################
###   Simulation function # 4   ###
###################################

####################################################################################
#
# Function name: Generate genome specific expression DEseq2 matrix (Gen_GSE_DEseq2_mat)
#
# Description: 
#   This function will iteratively call Gen_GSE_Data to create 
#   the matrixes ready for eQTL analysis by using Deseq2 package
#   
# Parameters (default value):
#   All the parameters have the same definitions as previous function.
# (1) mu_A_qq:          mean expression level in genome A for genotype qq (100)
# (2) mu_B_qq:          mean expression level in genome B for genotype qq (100)
# (3) fd_A:             flod changs in genome A (1)
# (4) fd_B:             flod changs in genome B (3)
# (5) phi:              overdispersion parameter (1)
# (6) ss:               sample size in for the genotype qq (25)
# (7) len_tans:         length of the transcript in bps (10kb)
# (8) n_diff_sites:     number of different site(s) between all the genomes (1)
# (9) len_reads:        length of the sequencing short reads (200) 
# (10) genotypic_ratio: the sample size ratio for the genotype (1:2:1)
# (11) n_rep:           number of repeats (genes) to generate for DEseq2 analysis
#
# Output:
#
#   A dataframe of both simulated (true) counts and 3 different ways for assigning  
#   unknown origion counts for A and B genome with corresponding genotypes.  
#   Note that the genotypes are independent between A and B genome. 
#
# Dependance:
#   Gen_GSE_Data funciton
#
####################################################################################

Gen_GSE_DEseq2_mat <- function(mu_A_qq=100, mu_B_qq=100, fd_A=1, fd_B=3, phi=1,
                               ss=25, len_trans=10000, n_diff_sites=1,
                               len_reads=200, genotypic_ratio=c(1,2,1), n_rep=1000){
  
  matrix_true_A <- matrix_true_B <- 
    matrix_GSE_A_GTA <- matrix_GSE_A_GTB <- 
    matrix_GSE_B_GTA <- matrix_GSE_B_GTB <- 
    matrix_even_A_GTA <- matrix_even_A_GTB <-
    matrix_even_B_GTA <- matrix_even_B_GTB <-
    matrix_prop_A_GTA <- matrix_prop_A_GTB <-
    matrix_prop_B_GTA <- matrix_prop_B_GTB <- 
    matrix_summ_GTA <- matrix_summ_GTB <- 
    matrix(NA,ncol=sum(ss*genotypic_ratio),nrow=n_rep)
  
  # create 16 matrixes with number of rows are total sample size, 
  # and the numbers of columns are number of replicates to store 
  # a. ture counts in genome A
  # b. true counts in genome B
  # c. genome A specific count 
  # d. genome B specific count 
  # e. evenly assign UOC in genome A
  # h. evenly assign UOC in genome B
  # i. proportionally assign UOC in genome A 
  # l. proportionally assign UOC in genome B 
  # m. summation of homoe-allele expression counts according to genotype A
  # n. summation of homoe-allele expression counts according to genotype B
  
  for (i in 1:n_rep){
    # for each replicates, we will do:
    
    d <- Gen_GSE_Data(mu_A_qq=mu_A_qq, mu_B_qq=mu_B_qq, fd_A=fd_A, fd_B=fd_B,
                      phi=phi, ss=ss, len_trans=len_trans, n_diff_sites=n_diff_sites,
                      len_reads=len_reads, genotypic_ratio=genotypic_ratio)
    # call the Gen_GSE_Data function to generate the simulated read counts
    
    ind_GT_A <- order(d$Genotype_A)
    # get the genotype in genome A by its order. 
    # The Gen_GSE_Data will shuffle the genotypes, 
    # but we need store the counts by genotypes.
    
    matrix_true_A[i,] <- d[ind_GT_A,"count_true_A"]
    matrix_even_A_GTA[i,] <- d[ind_GT_A,"even_A"]
    matrix_even_B_GTA[i,] <- d[ind_GT_A,"even_B"]
    matrix_prop_A_GTA[i,] <- d[ind_GT_A,"prop_A"]
    matrix_prop_B_GTA[i,] <- d[ind_GT_A,"prop_B"]
    matrix_summ_GTA[i,] <- d[ind_GT_A,"sum_count"]
    matrix_GSE_A_GTA[i,] <- d[ind_GT_A,"GSE_A"]
    matrix_GSE_B_GTA[i,] <- d[ind_GT_A,"GSE_B"]
    # save the true counts and swapped counts in corresponding matrix by genotype
  
    ind_GT_B <- order(d$Genotype_B)
    # get the genotype in genome B by its order. 
    
    matrix_true_B[i,] <- d[ind_GT_B,"count_true_B"]
    matrix_even_A_GTB[i,] <- d[ind_GT_B,"even_A"]
    matrix_even_B_GTB[i,] <- d[ind_GT_B,"even_B"]
    matrix_prop_A_GTB[i,] <- d[ind_GT_B,"prop_A"]
    matrix_prop_B_GTB[i,] <- d[ind_GT_B,"prop_B"]
    matrix_summ_GTB[i,] <- d[ind_GT_B,"sum_count"]
    matrix_GSE_A_GTB[i,] <- d[ind_GT_B,"GSE_A"]
    matrix_GSE_B_GTB[i,] <- d[ind_GT_B,"GSE_B"]
    # save the true counts and swapped counts in corresponding matrix by genotype
    
  }
  
  rownames(matrix_true_A) <- paste("True_A",1:n_rep,sep="_")   
  rownames(matrix_true_B) <- paste("True_B",1:n_rep,sep="_")
  rownames(matrix_GSE_A_GTA) <- paste("GSE_A_GTA",1:n_rep,sep="_")   
  rownames(matrix_GSE_A_GTB) <- paste("GSE_A_GTB",1:n_rep,sep="_")
  rownames(matrix_GSE_B_GTA) <- paste("GSE_B_GTA",1:n_rep,sep="_")   
  rownames(matrix_GSE_B_GTB) <- paste("GSE_B_GTB",1:n_rep,sep="_")
  rownames(matrix_even_A_GTA) <- paste("Even_A_GTA",1:n_rep,sep="_")
  rownames(matrix_even_A_GTB) <- paste("Even_A_GTB",1:n_rep,sep="_")
  rownames(matrix_even_B_GTA) <- paste("Even_B_GTA",1:n_rep,sep="_")
  rownames(matrix_even_B_GTB) <- paste("Even_B_GTB",1:n_rep,sep="_")
  rownames(matrix_prop_A_GTA) <- paste("Prop_A_GTA",1:n_rep,sep="_")   
  rownames(matrix_prop_A_GTB) <- paste("Prop_A_GTB",1:n_rep,sep="_")   
  rownames(matrix_prop_B_GTA) <- paste("Prop_B_GTA",1:n_rep,sep="_")
  rownames(matrix_prop_B_GTB) <- paste("Prop_B_GTB",1:n_rep,sep="_")
  rownames(matrix_summ_GTA) <- paste("Summ_GTA",1:n_rep,sep="_")
  rownames(matrix_summ_GTB) <- paste("Summ_GTB",1:n_rep,sep="_")
  # assign the row names for each matrix 
  # True_A: true counts for genome A
  # True_B: true counts for genome B
  # GSE_A: genome A specific counts
  # GSE_B: genome B specific counts
  # Even_A: evenly assign UOC for genome A
  # Even_B: evenly assign UOC for genome B
  # Prop_A: proportionally assign UOC for genome A
  # Prop_B: proportionally assign UOC for genome B
  # Summ_A: summation UOC for genome A
  # Summ_B: summation UOC for genome B
  
  gene_counts <- as.matrix(rbind(matrix_true_A, matrix_true_B, 
                                 matrix_GSE_A_GTA, matrix_GSE_A_GTB,
                                 matrix_GSE_B_GTA, matrix_GSE_B_GTB,
                                 matrix_even_A_GTA, matrix_even_A_GTB,
                                 matrix_even_B_GTA, matrix_even_B_GTB,
                                 matrix_prop_A_GTA, matrix_prop_A_GTB,
                                 matrix_prop_B_GTA, matrix_prop_B_GTB,
                                 matrix_summ_GTA, matrix_summ_GTB))
  # create a merged matrix with all the counts data ready for DEseq2 analysis
  
  GT_names <- paste("GT",1:sum(ss*genotypic_ratio),sep="")
  # create a systematic names for each sample
  colnames(gene_counts) <- GT_names
  # assign the column names for gene count table 
  
  col_data <- data.frame(Genotype=factor(c(rep("qq",ss*genotypic_ratio[1]),
                                           rep("Qq",ss*genotypic_ratio[2]),
                                           rep("QQ",ss*genotypic_ratio[3])),
                                         levels=c("qq","Qq","QQ")))
  # create the required column data frame for DEseq2 analysis
  rownames(col_data) <- GT_names
  # assign row names for column data frame
  
  return(list(gene_counts=gene_counts,
              colData=col_data))
  # return 2 matrixes both ready for DEseq2 to do analysis
}

#####################################################################################
###                                                                               ###
### determine the lose of power and increase of false postive                     ###
### corresponding to 3 different way separating unknown origin counts (UOC).      ###
###                                                                               ###
#####################################################################################

library(DESeq2)
library(MASS)
# load required packages for the simulation

n_rep <- 1000
# set numbers of replications 
ss <- 25
# numbers of sample size for genotype qq
fd <- 1.75
# set the fold change in genome B to 2

mu_A_qq <- mu_B_qq <- 25
# assign mean expression value for both A and B genome from 20 to 200
n_diff_sites <- c(1,3,5)
# number of SNPs between A & B genome 
len_trans <- 5000

paras <- expand.grid(mu_A_qq=mu_A_qq, mu_B_qq=mu_B_qq, n_rep=n_rep,
                     n_diff_sites=n_diff_sites, fd=fd, ss=ss,
                     len_trans=len_trans)
# create a matrix with all different combinations of parameters
nrow(paras)

FDR_True_A <- Power_True_B <- 
  FDR_GSE_A_GTA <- FDR_GSE_A_GTB <- Power_GSE_B_GTA <- Power_GSE_B_GTB <- 
  FDR_Even_A_GTA <- FDR_Even_A_GTB <- Power_Even_B_GTA <- Power_Even_B_GTB <-
  FDR_Prop_A_GTA <- FDR_Prop_A_GTB <- Power_Prop_B_GTA <- Power_Prop_B_GTB <- 
  FDR_Summ_A <-  Power_Summ_B <- rep(-1,nrow(paras))
# create 16 vectors to store the results of 
# false discovery rate and power before FDR adjustment

adj_FDR_True_A <- adj_Power_True_B <- 
  adj_FDR_GSE_A_GTA <- adj_FDR_GSE_A_GTB <- 
  adj_Power_GSE_B_GTA <- adj_Power_GSE_B_GTB <- 
  adj_FDR_Even_A_GTA <- adj_FDR_Even_A_GTB <- 
  adj_Power_Even_B_GTA <- adj_Power_Even_B_GTB <-
  adj_FDR_Prop_A_GTA <- adj_FDR_Prop_A_GTB <- 
  adj_Power_Prop_B_GTA <- adj_Power_Prop_B_GTB <- 
  adj_FDR_Summ_A <-  adj_Power_Summ_B <- rep(-1,nrow(paras))
# create 16 vectors to store the results of 
# false discovery rate and power after FDR adjustment

ind <- matrix(1:(16*n_rep),ncol=16,byrow=FALSE)

for(i in 1:nrow(paras)){
  # for each parameter setting
  
  print(i)
  # print current iteration to track process
  
  temp <- Gen_GSE_DEseq2_mat(n_rep =  paras[i,"n_rep"], ss = paras[i,"ss"], 
                             mu_A_qq = paras[i,"mu_A_qq"],mu_B_qq = paras[i,"mu_B_qq"], 
                             fd_B = paras[i,"fd"],n_diff_sites = paras[i,"n_diff_sites"],
                             len_trans = paras[i,"len_trans"])
  
  # call Gen_GSE_DEseq2_mat to generate 2 matrixes for DEseq2 analysis
  
  dds <- DESeqDataSetFromMatrix(countData=temp$gene_counts,
                                colData=temp$colData,
                                design= ~ Genotype)
  # DESeqDataSetFromMatrix function from DEseq2 to set experiment design
  
  normalFactors <- matrix(1,ncol=ncol(dds),nrow=nrow(dds),
                          dimnames = list(1:nrow(dds),1:ncol(dds)))
  # make an idendity matrix as the nomalization factors
  
  normalizationFactors(dds) <- normalFactors
  # assign the nomalization factor to be 1 
  
  DE_dds <- DESeq(dds, test="LRT", reduced = ~ 1, fitType="mean") 
  # find the differencial expressed genes by design (genotype)
  
  res_DE_dds <- results(DE_dds)
  # get the results form the model 
  
  ng <- nrow(res_DE_dds)
  # get number of rows of the result matrix
  
  FDR_True_A[i] <- sum(res_DE_dds$pvalue[ind[,1]]<0.05,na.rm = TRUE)/n_rep
  Power_True_B[i] <- sum(res_DE_dds$pvalue[ind[,2]]<0.05,na.rm = TRUE)/n_rep
  FDR_GSE_A_GTA[i] <- sum(res_DE_dds$pvalue[ind[,3]]<0.05,na.rm = TRUE)/n_rep
  FDR_GSE_A_GTB[i] <- sum(res_DE_dds$pvalue[ind[,4]]<0.05,na.rm = TRUE)/n_rep
  Power_GSE_B_GTA[i] <- sum(res_DE_dds$pvalue[ind[,5]]<0.05,na.rm = TRUE)/n_rep
  Power_GSE_B_GTB[i] <- sum(res_DE_dds$pvalue[ind[,6]]<0.05,na.rm = TRUE)/n_rep
  FDR_Even_A_GTA[i] <- sum(res_DE_dds$pvalue[ind[,7]]<0.05,na.rm = TRUE)/n_rep
  FDR_Even_A_GTB[i] <- sum(res_DE_dds$pvalue[ind[,8]]<0.05,na.rm = TRUE)/n_rep
  Power_Even_B_GTA[i] <- sum(res_DE_dds$pvalue[ind[,9]]<0.05,na.rm = TRUE)/n_rep
  Power_Even_B_GTB[i] <- sum(res_DE_dds$pvalue[ind[,10]]<0.05,na.rm = TRUE)/n_rep
  FDR_Prop_A_GTA[i] <- sum(res_DE_dds$pvalue[ind[,11]]<0.05,na.rm = TRUE)/n_rep
  FDR_Prop_A_GTB[i] <- sum(res_DE_dds$pvalue[ind[,12]]<0.05,na.rm = TRUE)/n_rep
  Power_Prop_B_GTA[i] <- sum(res_DE_dds$pvalue[ind[,13]]<0.05,na.rm = TRUE)/n_rep
  Power_Prop_B_GTB[i] <- sum(res_DE_dds$pvalue[ind[,14]]<0.05,na.rm = TRUE)/n_rep
  FDR_Summ_A[i] <- sum(res_DE_dds$pvalue[ind[,15]]<0.05,na.rm = TRUE)/n_rep
  Power_Summ_B[i] <- sum(res_DE_dds$pvalue[ind[,16]]<0.05,na.rm = TRUE)/n_rep
  # save the results of false discovery rate and power before the FDR adjustment

  adj_FDR_True_A[i] <- sum(res_DE_dds$padj[ind[,1]]<0.05,na.rm = TRUE)/n_rep
  adj_Power_True_B[i] <- sum(res_DE_dds$padj[ind[,2]]<0.05,na.rm = TRUE)/n_rep
  adj_FDR_GSE_A_GTA[i] <- sum(res_DE_dds$padj[ind[,3]]<0.05,na.rm = TRUE)/n_rep
  adj_FDR_GSE_A_GTB[i] <- sum(res_DE_dds$padj[ind[,4]]<0.05,na.rm = TRUE)/n_rep
  adj_Power_GSE_B_GTA[i] <- sum(res_DE_dds$padj[ind[,5]]<0.05,na.rm = TRUE)/n_rep
  adj_Power_GSE_B_GTB[i] <- sum(res_DE_dds$padj[ind[,6]]<0.05,na.rm = TRUE)/n_rep
  adj_FDR_Even_A_GTA[i] <- sum(res_DE_dds$padj[ind[,7]]<0.05,na.rm = TRUE)/n_rep
  adj_FDR_Even_A_GTB[i] <- sum(res_DE_dds$padj[ind[,8]]<0.05,na.rm = TRUE)/n_rep
  adj_Power_Even_B_GTA[i] <- sum(res_DE_dds$padj[ind[,9]]<0.05,na.rm = TRUE)/n_rep
  adj_Power_Even_B_GTB[i] <- sum(res_DE_dds$padj[ind[,10]]<0.05,na.rm = TRUE)/n_rep
  adj_FDR_Prop_A_GTA[i] <- sum(res_DE_dds$padj[ind[,11]]<0.05,na.rm = TRUE)/n_rep
  adj_FDR_Prop_A_GTB[i] <- sum(res_DE_dds$padj[ind[,12]]<0.05,na.rm = TRUE)/n_rep
  adj_Power_Prop_B_GTA[i] <- sum(res_DE_dds$padj[ind[,13]]<0.05,na.rm = TRUE)/n_rep
  adj_Power_Prop_B_GTB[i] <- sum(res_DE_dds$padj[ind[,14]]<0.05,na.rm = TRUE)/n_rep
  adj_FDR_Summ_A[i] <- sum(res_DE_dds$padj[ind[,15]]<0.05,na.rm = TRUE)/n_rep
  adj_Power_Summ_B[i] <- sum(res_DE_dds$padj[ind[,16]]<0.05,na.rm = TRUE)/n_rep
  # save the results of false discovery rate and power after the FDR adjustment
  
}

paras$MEV_ratio <- paras$mu_A_qq/paras$mu_B_qq

paras$col_mu <- ifelse(paras$mu_A==paras$mu_B, 1,ifelse(paras$mu_A>paras$mu_B,2,4))
# set the color to represent the ratio of the mu_A and mu_B
# black means mu_A = mu_B
# red means mu_A > mu_B
# blue means mu_A < mu_B
paras$col_fd <- sapply(paras$fd, function(x){which(x==unique(paras$fd))})
# set the color to represent the fold changes for genome B
# set black to fd = 1.5; red to fd = 1.75; green to fd = 2; 
# blue to fd = 2.25; light blue fd = 2.5

Four_way_table <- cbind(paras,FDR_True_A, Power_True_B,
                         FDR_GSE_A_GTA, FDR_GSE_A_GTB, 
                         Power_GSE_B_GTA, Power_GSE_B_GTB,
                         FDR_Even_A_GTA, FDR_Even_A_GTB,
                         Power_Even_B_GTA, Power_Even_B_GTB, 
                         FDR_Prop_A_GTA, FDR_Prop_A_GTB,
                         Power_Prop_B_GTA, Power_Prop_B_GTB, FDR_Summ_A, Power_Summ_B,
                         adj_FDR_True_A, adj_Power_True_B, 
                         adj_FDR_GSE_A_GTA, adj_FDR_GSE_A_GTB, 
                         adj_Power_GSE_B_GTA, adj_Power_GSE_B_GTB,
                         adj_FDR_Even_A_GTA, adj_FDR_Even_A_GTB,
                         adj_Power_Even_B_GTA, adj_Power_Even_B_GTB,
                         adj_FDR_Prop_A_GTA, adj_FDR_Prop_A_GTB,
                         adj_Power_Prop_B_GTA, adj_Power_Prop_B_GTB,
                         adj_FDR_Summ_A, adj_Power_Summ_B)

# setwd("~/Desktop/Manuscript_1/Appendix_Material/4_Expected_Number_of_ASE_Reads_by_Given_Number_of_SNPs/")

rm(list=setdiff(ls(),"Four_way_table"))

write.table(Four_way_table, "~/Desktop/Depth=1_Simulation.csv", quote = FALSE, row.names= FALSE, sep=",")

# save.image("~/Desktop/temp/F175_1kb.RData")

save.image("Simulate_Power_and_False_Discovery_by_assign_UOC_four_ways.RData")

