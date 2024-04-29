
rm(list=ls())

library(MASS)
# rnegbin function in MASS package is needed for generate random negative binomial variables.

########################
## A function to randomly generate gene expression levels (# of short reads)
## for corresponding genotype from nagative binomial distribution by given 
# (1) the mean number of reads for genotype AA
# (2) the fold change (effect size)
# (3) over-dispersion parameter for negative binomial
# (4) sample size
#     if there is only one element assigned, then the function will amplify it with  
#     the assigned genotypic ratio 
# (5) genotypic_ratio 
#     the sample size ratio for the genotype. The default is 1:2:1,
#     which follows Hardy Weinberg genotypic ratio
########################

genotype_explevel <- function(mu_AA=100,fd=1,phi=1,ss=1,genotypic_ratio=c(1,2,1)){ 
  
  if(length(ss)==1){ ss <- ss*genotypic_ratio }
  # check if there is only one element assign in ss,
  # if so, amplify it with the genotypic ratio.
  
  if(length(ss)!=3){
    break("There should be 3 sample saizes for each genotype")
  }
  # check if the sample sizes is properly assigned
  
  # corresponding mean value for each genotype
  mu_X <- c(mu_AA,mu_AA*fd,mu_AA*(2*fd-1))
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

########################
## A function to return number of transcripts form either from A genome, 
## from B genome, or from the genome that we can't tell by given
# (1) length of the transcript, default is 10kb
# (2) number of different sites between A and B genome, default is 10
# (3) number of reads mapped to each transcripts. 
#      if there is only one element assigned, 
#      then the function will assume equal read number for all genome 
# (4) length of the sequencing short reads for each transcripts. 
#      if there is only one element assigned, then the function will 
#     assume equal read length for all genome
# (5) copy of the genome, default is 2 for allotetraploid (A and B genome)
########################

Transcript_Counts_freq <- function(len_trans=10000,
                                   n_diff_sites=10,
                                   n_reads=100,
                                   len_reads=200,
                                   copy_genome=2)
{ 
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

########################
### A function generate random data set for estimating Zab & Zba
### This function combines two previous functions (genotype_explevel & Transcript_Counts_freq)
### Therefore all the parameters have the same definitions as those two functions.
########################

Gen_Sim_Data <- function(mu_1=100,mu_2=100,fd_A=1,fd_B=3,phi=1,ss_1=25,ss_2=25,
                         len_trans=10000,n_diff_sites=10,len_reads=200,
                         genotypic_ratio=c(1,2,1)){

  ss_gen_A <- genotype_explevel(mu_1,fd_A,phi,ss_1,genotypic_ratio)
  ss_gen_B <- genotype_explevel(mu_2,fd_B,phi,ss_2,genotypic_ratio)
  # randomly generate counts for two genome transcripts using "genotype_explevel" function
  
  d_gen_A <- data.frame(Genotype_A=c(rep("AA",length(ss_gen_A[[1]])),
                                     rep("AB",length(ss_gen_A[[2]])),
                                     rep("BB",length(ss_gen_A[[3]]))),
                        count_true_A=c(ss_gen_A[[1]],ss_gen_A[[2]],ss_gen_A[[3]]))
  
  d_gen_B <- data.frame(Genotype_B=c(rep("AA",length(ss_gen_B[[1]])),
                                     rep("AB",length(ss_gen_B[[2]])),
                                     rep("BB",length(ss_gen_B[[3]]))),
                        count_true_B=c(ss_gen_B[[1]],ss_gen_B[[2]],ss_gen_B[[3]]))
  # expand the output of the "genotype_explevel" function 
  
  nr <- nrow(d_gen_A)
  
  ran_d_gen_A <- d_gen_A[sample(1:nr),]
  ran_d_gen_B <- d_gen_B[sample(1:nr),]
  # randomly shuffle the data 
  
  ran_d <- cbind(ran_d_gen_A,ran_d_gen_B)
  # combine two transcripts counts into one data frame
  
  counts <- ran_d[,c("count_true_A","count_true_B")]
  # extract transcript counts for "Transcript_Counts_freq" function
  
  reads_genomes <- t(apply(counts,1,function(x){
    temp <- Transcript_Counts_freq(n_reads=x,len_trans,n_diff_sites,
                                   len_reads)$reads_genome
    return(c(temp[,2],sum(temp[,1])))}
  ))
  # simulate reads could and coouldn't be identified from "Transcript_Counts_freq" function
  
  d <- cbind(ran_d,reads_genomes)
  names(d)[5:7] <- c("GS_counts_A","GS_counts_B","Uncertain_counts")
  # combine all the information into a single data frame
  
  return(d)
}

########################
###
### Estimate Zab & Zba
###
########################

### a small function to avoid zero while calculating proportions ###

av_ze <- function(x){
  sapply(x, function(y){ if(y==0) {y <- 1} else {y <- y} })
}

### a small function to calculate Pa & Pb

Pa_Pb <- function(TA, TB, SA, SB){
  if(TA==0 & TB==0){
    Pa <- Pb <- 0
  }else if(TA==0){
    Pa <- 0
    Pb <- (SB-TB)/TB
  }else if(TB==0){
    Pa <- (SA-TA)/TA
    Pb <- 0
  }else{
    Pa <- (SA-TA)/TA
    Pb <- (SB-TB)/TB
  }
  return(c(Pa=Pa,Pb=Pb))
}

### Calculate Zab & Pb according to different parameters settings

mu_1 <- c(100, 200)
mu_2 <- c(100, 200)
fd_B <- 2
n_diff_sites <- c(1,5,10)
ss <- 1600

paras <- expand.grid(mu_1=mu_1, mu_2=mu_2, fd_B=fd_B, 
                     ss=ss, n_diff_sites=n_diff_sites)

# nrow(paras)

paras$mean_Pa_prop <- paras$mean_Pb_prop <- paras$mean_Pa_even <- paras$mean_Pb_even <- -1

PaPbs_prop <- PaPbs_even <- vector("list",nrow(paras))

for (i in 1:nrow(paras)){

  print(i)
  
  d <- Gen_Sim_Data(mu_1=paras[i,1], mu_2=paras[i,2], fd_B=paras[i,3],
                    ss_1=paras[i,4], ss_2=paras[i,4], 
                    n_diff_sites=paras[i,5],len_trans=1000)

  d$prop_A <- av_ze(d$GS_counts_A)/(av_ze(d$GS_counts_A)+av_ze(d$GS_counts_B))
  d$prop_B <- 1-d$prop_A
    
  d$Sim_prop_counts_A <- d$GS_counts_A+round(d$prop_A*d$Uncertain_counts)
  d$Sim_prop_counts_B <- d$GS_counts_B+round(d$prop_B*d$Uncertain_counts)

  d$Sim_even_counts_A <- d$GS_counts_A+round(d$Uncertain_counts/2,0)
  d$Sim_even_counts_B <- d$GS_counts_B+round(d$Uncertain_counts/2,0)
  
  temp <- d[,c("count_true_A","count_true_B",
               "Sim_prop_counts_A","Sim_prop_counts_B",
               "Sim_even_counts_A","Sim_even_counts_B")]
  
  Zabba_prop <- do.call(rbind,lapply(as.data.frame(t(temp)),
                                     function(x){Pa_Pb(x[1],x[2],x[3],x[4])}))
  
  Zabba_even <- do.call(rbind,lapply(as.data.frame(t(temp)),
                                     function(x){Pa_Pb(x[1],x[2],x[5],x[6])}))
  
  PaPbs_prop[[i]] <- data.frame(Zabba_prop,GB=d$Genotype_B)
  PaPbs_even[[i]] <- data.frame(Zabba_even,GB=d$Genotype_B)
  
  paras$mean_Pa_prop[i] <- mean(Zabba_prop[,"Pa"])  
  paras$mean_Pb_prop[i] <- mean(Zabba_prop[,"Pb"])

  paras$mean_Pa_even[i] <- mean(Zabba_even[,"Pa"])  
  paras$mean_Pb_even[i] <- mean(Zabba_even[,"Pb"])
  
}

est_Pa_Pb <- paras

est_Pa_Pb$mu_ratio <- est_Pa_Pb$mu_1/est_Pa_Pb$mu_2
  
rm(list=setdiff(ls(),c("PaPbs_prop","PaPbs_even","est_Pa_Pb")))

setwd("~/Desktop/Draft_of_first_manuscript/Figures/Scripts_&_data_for_Figures/")
save.image("Estimating_Pa_and_Pb_by_given_number_of_SNPs.RData")



