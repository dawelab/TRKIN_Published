library(ggplot2)
library(reshape2)

setwd("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/Paper/TRKIN_Published/trkin1and2_Expression")

Counts <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/Differential_Expression/FeatureCounts_counts.txt", sep="")
Anno <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/Differential_Expression/FeatureCounts_annotation.txt", sep="")

Counts$chr <- Anno$Chr
Counts$start <- Anno$Start
Counts$end <- Anno$End
Counts$strand <- Anno$Strand
Counts$Length <- Anno$Length 

Counts$LengthKB <- Counts$Length/1000

tpm3 <- function(counts,len) {x <- counts/len
return(t(t(x)*1e6/colSums(x)))
}
TPM <- tpm3(Counts[,1:20],Counts$LengthKB)
TPM_DF <- as.data.frame(TPM)
TPM_DF$gene_ID <- rownames(Counts)
TPM_DF <- separate(TPM_DF, gene_ID, into = c("gene_ID", "exon"), sep=".exon")
colnames(TPM_DF) <- c("embryo_1", "embryo_2", "endosperm_1", "endosperm_2", "root_1", "root_2", "shoot_1", "shoot_2", "anther_1", "anther_2", "base_1", "base_2", "middle_1", "middle_2", "tip_1", "tip_2", "ear_1", "ear_2", "tassel_1", "tassel_2", "gene_ID", "exon")


TRKIN_TPM <- subset(TPM_DF, gene_ID=="g5486.t1")
TRKIN_MELT <- melt(TRKIN_TPM, id=c("gene_ID", "exon"))
TRKIN_MELT$Tissue <- c(rep("embryo", 19*2), rep("endosperm", 19*2), rep("root", 19*2), rep("shoot", 19*2), rep("anther", 19*2), rep("leaf base", 19*2), rep("leaf middle", 19*2), rep("leaf tip", 19*2), rep("ear", 19*2), rep("tassel", 19*2))

TRKIN_MELT_SUB <- subset(TRKIN_MELT, exon == 7 | exon == 8)

PTRKIN_TPM <- subset(TPM_DF, gene_ID=="g5491.t1")
PTRKIN_MELT <- melt(PTRKIN_TPM, id=c("gene_ID", "exon"))
PTRKIN_MELT$Tissue <- c(rep("embryo", 19*2), rep("endosperm", 19*2), rep("root", 19*2), rep("shoot", 19*2), rep("anther", 19*2), rep("leaf base", 19*2), rep("leaf middle", 19*2), rep("leaf tip", 19*2), rep("ear", 19*2), rep("tassel", 19*2))

PTRKIN_MELT_SUB <- subset(PTRKIN_MELT, exon == 7 | exon == 8)

#T test tooking fir difference in values
t.test(TRKIN_MELT_SUB$value, PTRKIN_MELT_SUB$value)

#Calculating reduction
((mean(TRKIN_MELT_SUB$value) - mean(PTRKIN_MELT_SUB$value))  / mean(TRKIN_MELT_SUB$value))*100


TRKIN_MELT_SUB$gene_ID <- "trkin1"
PTRKIN_MELT_SUB$gene_ID <- "trkin2"
ALL <- rbind(TRKIN_MELT_SUB, PTRKIN_MELT_SUB)

pdf(file="All.pdf", height=6, width=4)
ggplot(ALL, aes(x=exon, y=value, color=gene_ID)) +
  geom_point(alpha=0.5) +
  scale_colour_manual(values=c("#0072b2", "orange")) +
  facet_wrap(~Tissue, ncol=2) +
  ylab("TPM") +
  ggtitle("") +
  scale_y_continuous(limits=c(0,5), breaks=c(0,5)) +
  theme(plot.title = element_text(size =20, hjust = 0.5), axis.text.y = element_text(size =15),axis.text.x = element_text(size =20, angle = 90, vjust = 0.5, hjust=1), axis.title = element_text(size=15), strip.text.x = element_text(size=13), legend.position = "right")
dev.off()

