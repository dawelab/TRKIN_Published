library(tidyverse)
library(ggplot2)
library(reshape2)

setwd("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/Paper/TRKIN_Paper")

#Load the data
COUNTS <- read.csv("/Volumes/Transcend/FeatureCounts_counts.txt", sep="")
COUNTS$GeneID <- rownames(COUNTS)
ANN <- read.csv("/Volumes/Transcend/FeatureCounts_annotation.txt", sep="")
OF_ANN <- read.delim("/Volumes/Transcend/Zm-B73_AB10-REFERENCE-NAM-1.0_Zm00043a.1.noheader.gff3", header=FALSE, comment.char="#")

#merge the Counts with the gene length
COUNT_LEN <- merge(COUNTS, ANN[c("GeneID", "Length")], by="GeneID")

COUNT_LEN$Length_kb <- COUNT_LEN$Length/1000

tpm3 <- function(counts,len) {x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

colnames(COUNT_LEN)

TPM <- tpm3(COUNT_LEN[,2:15],COUNT_LEN$Length_kb)
TPM_DF <- as.data.frame(TPM)
TPM_DF$GeneID <- COUNT_LEN$GeneID

write.csv(TPM_DF, "K10L2_Ab10I_N10_TranscriptsPerMillion.csv")

TRKIN_TPM <- subset(TPM_DF, GeneID=="gene:Zm00043a049192")

TRKIN_MELT <- melt(TRKIN_TPM)

TRKIN_MELT$Gene <- "trkin 1"

TRKIN_MELT$Name <- c("Ab10\n1", "Ab10\n2", "Ab10\n3", "K10L2\n1", "K10L2\n2", "K10L2\n3", "K10L2\n4", "K10L2\n5", "N10\n1", "N10\n2", "N10\n3", "N10\n4", "N10\n5", "N10\n6")

PTRKIN_TPM <- subset(TPM_DF, GeneID=="gene:Zm00043a049309")

PTRKIN_MELT <- melt(PTRKIN_TPM)
PTRKIN_MELT$Gene <- "trkin 2"
PTRKIN_MELT$Name <- c("Ab10\n1", "Ab10\n2", "Ab10\n3", "K10L2\n1", "K10L2\n2", "K10L2\n3", "K10L2\n4", "K10L2\n5", "N10\n1", "N10\n2", "N10\n3", "N10\n4", "N10\n5", "N10\n6")

#This adds a sum dataframe
SUM <- PTRKIN_MELT
SUM$GeneID <- "gene:Zm00043a049192 	gene:Zm00043a049309"
SUM$value <- TRKIN_MELT$value+PTRKIN_MELT$value
SUM$Gene <- "sum"

#This brings together the trkin 1 and trkin 2 datasources
ALL_trkin <- rbind(TRKIN_MELT[,c(2,3,4,5)], PTRKIN_MELT[,c(2,3,4,5)], SUM[,c(2,3,4,5)])

pdf("Ab10trkinsTPM.pdf", height = 5, width=11)
ggplot(ALL_trkin, aes(x=Name, y=value, fill=Gene)) +
  geom_bar(position='dodge', stat="identity") +
  labs(y="TPM", x="Chromosome 10 Haplotype") +
  theme(axis.text.x = element_text()) +
  scale_y_continuous(breaks=c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3)) +
  scale_fill_viridis_d(option = "mako", begin = 0.2, end = 0.7) +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 15), legend.title = element_text(size = 18), legend.text = element_text(size = 15), legend.position = "bottom")
dev.off()
