
rm(list=ls())
setwd("~/Desktop/Draft_of_first_manuscript/Figures/Scripts_&_data_for_Figures/")

##########

len_align <- scan("align_length.txt")
gap_open <- scan("gap_openings.txt")
mis_matches <- scan("mismatches.txt")

##########

diff_per_1k_align <- (mis_matches+gap_open)/len_align*1000
mismatches_per_1k_align <- mis_matches/len_align*1000
gap_opening_per_1k_align <- gap_open/len_align*1000

x <- summary(diff_per_1k_align)

pdf("../Figure_3.pdf")

par(mfrow=c(3,1))

hist(mismatches_per_1k_align,breaks=seq(-0.5,222.5,1),
     xlab="Numbers of mismatches per 1K alignment",
     main="Next Best BLAST Hit in Switchgrss",las=1)
mtext("(a)",side=2, line=0, las=1, at=10000)
text(200,5000,"N = 74153")

hist(gap_opening_per_1k_align,breaks=seq(-0.5,55.5,1),
     xlab="Numbers of gap opening per 1K alignment",
     main="",las=1)
mtext("(b)",side=2, line=0, las=1, at=37000)

hist(diff_per_1k_align,breaks=seq(-0.5,222.5,1),
     xlab="Numbers of differences per 1K alignment",
     main="",las=1)
mtext("(c)",side=2, line=0, las=1, at=10000)
legend(100,7000,legend=paste(names(x),x),bty="n")


dev.off()

#####

sum(mismatches_per_1k_align==0)/length(mismatches_per_1k_align)

sum(diff_per_1k_align==0)/length(diff_per_1k_align)
sum(diff_per_1k_align>0 & diff_per_1k_align<11)/length(diff_per_1k_align)
sum(diff_per_1k_align>10 & diff_per_1k_align<21)/length(diff_per_1k_align)
sum(diff_per_1k_align>20)/length(diff_per_1k_align)

summary(diff_per_1k_align)
summary(mismatches_per_1k_align)
summary(gap_opening_per_1k_align)
