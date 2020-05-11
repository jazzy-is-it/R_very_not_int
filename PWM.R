setwd("~/Desktop/R_zettel_2")
library(plyr)
library(stringi)
library(stringr)
library(big.matrix)
library(biganalytics)

#read table with alignment, factor=F! WHY NORA?
seq <- read.table("seq.txt", stringsAsFactors = F)
#data frames are easy
seq <- as.data.frame(seq)
seq




#vector which will be multiplied by 4 by the loop and will save all the count from the different letters
count <- vector()

#a vector where the 4 dna letters are, so the loop works
nucleot <- c("A","T","G","C")

#make vectors where save all the counts..the vector count will be multilied and A T G C will be added
for(i in nucleot) {
  assign(paste("count",i,sep="_"), count)
}

#loop for counting all the As Ts Cs Gs
for (c in 1:dim(seq)[2]){
#in As the counts of As from every column will be saved
#As <- str_count(seq[ ,c], "A")
#in numA the sum of the counts of As will be saved 
numA <- sum(str_count(seq[seq[,c] == "A", c]))
#Ts <- str_count(seq[ ,c], "T")
numT <- sum(str_count(seq[seq[,c] == "T", c]))
#Gs <- str_count(seq[ ,c], "G")
numG <- sum(str_count(seq[seq[,c] == "G", c]))
#Cs <- str_count(seq[ ,c], "C")
numC <- sum(str_count(seq[seq[,c] == "C", c]))
#the empty vector count_A will be bound to numA, where the sum of the counts from A is
count_A <- cbind(count_A, c(numA))
count_T <- cbind(count_T, c(numT))
count_G <- cbind(count_G, c(numG))
count_C <- cbind(count_C, c(numC))
}

#bind together all sum of all counts 
table <- rbind(count_A, count_T, count_C, count_G)
#add names to columns
colnames(table) <- c("1", "2", "3", "4", "5", "6", "7")
#add names to rows
rownames(table) <- c("A", "T", "C", "G")
#check the class of table
class("table")
#make table a matrix and give it the name PCM Position count matrix
pcm <- as.matrix(table)
pcm

#add pseudo count of 1 to every string in the matrix and call the new matrix pseudo
pseudo <- pcm + 1

#make it a PFM, where u divied it by 19(not 15 because of the pseudo)
pfm <- pseudo/19

#make the pfm in percent
perc_pfm <- pfm*100

#based on the GC-percentage in the organism, divide the observed frequncy, which is the pfm matrix, by the biological frequency
#likelihiid ratio, meaning how many times the sequence differs from a rnd sequence
p_log_ratioA <- pfm[1,]/0.3
p_log_ratioT <- pfm[2,]/0.3
p_log_ratioC <- pfm[3,]/0.2
p_log_ratioG <- pfm[4,]/0.2

#make a vecrtor where to save the ration between observed and real frequency (likelihood ration)
all_p_likeli_ratio <- vector()
all_p_likeli_ratio <- rbind(p_log_ratioA, p_log_ratioT, p_log_ratioC, p_log_ratioG)

#make the pwm, by taking the log2 of the likelihood ratio
pwm <- log2(all_p_likeli_ratio)
rownames(pwm) <-  c("A", "T", "C", "G")

#find the max number in every column
max <-as.big.matrix(pwm)
max
colmax(max)

#reversse column and rows 
test <- t(pwm)


