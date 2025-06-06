library(ggplot2)
library(reshape2)

setwd("")

#From step4
Counts <- read.csv("FeatureCounts_counts.txt", sep="")
Anno <- read.csv("FeatureCounts_annotation.txt", sep="")

Counts$chr <- Anno$Chr
Counts$start <- Anno$Start
Counts$end <- Anno$End
Counts$strand <- Anno$Strand
Counts$Length <- Anno$Length 

Counts$LengthKB <- Counts$Length/1000

#Calculate tpm
tpm3 <- function(counts,len) {x <- counts/len
return(t(t(x)*1e6/colSums(x)))
}
TPM <- tpm3(Counts[,1:20],Counts$LengthKB)
TPM_DF <- as.data.frame(TPM)
TPM_DF$gene_ID <- rownames(Counts)
TPM_DF <- separate(TPM_DF, gene_ID, into = c("gene_ID", "exon"), sep=".exon")
colnames(TPM_DF) <- c("embryo_1", "embryo_2", "endosperm_1", "endosperm_2", "root_1", "root_2", "shoot_1", "shoot_2", "anther_1", "anther_2", "base_1", "base_2", "middle_1", "middle_2", "tip_1", "tip_2", "ear_1", "ear_2", "tassel_1", "tassel_2", "gene_ID", "exon")


TRKIN1_TPM <- subset(TPM_DF, gene_ID=="g5486.t1")
TRKIN1_MELT <- melt(TRKIN1_TPM, id=c("gene_ID", "exon"))
TRKIN1_MELT$Tissue <- c(rep("embryo", 19*2), rep("endosperm", 19*2), rep("root", 19*2), rep("shoot", 19*2), rep("anther", 19*2), rep("leaf base", 19*2), rep("leaf middle", 19*2), rep("leaf tip", 19*2), rep("ear", 19*2), rep("tassel", 19*2))

TRKIN1_MELT_SUB <- subset(TRKIN1_MELT, exon == 7 | exon == 8)

TRKIN2_TPM <- subset(TPM_DF, gene_ID=="g5491.t1")
TRKIN2_MELT <- melt(TRKIN2_TPM, id=c("gene_ID", "exon"))
TRKIN2_MELT$Tissue <- c(rep("embryo", 19*2), rep("endosperm", 19*2), rep("root", 19*2), rep("shoot", 19*2), rep("anther", 19*2), rep("leaf base", 19*2), rep("leaf middle", 19*2), rep("leaf tip", 19*2), rep("ear", 19*2), rep("tassel", 19*2))

TRKIN2_MELT_SUB <- subset(TRKIN2_MELT, exon == 7 | exon == 8)

#T test tooking fir difference in values
t.test(TRKIN1_MELT_SUB$value, TRKIN1_MELT_SUB$value)

#Calculating reduction
((mean(TRKIN1_MELT_SUB$value) - mean(TRKIN2_MELT_SUB$value))  / mean(TRKIN1_MELT_SUB$value))*100

#Label the genes for plotting
TRKIN1_MELT_SUB$gene_ID <- "trkin1"
TRKIN2_MELT_SUB$gene_ID <- "trkin2"
ALL <- rbind(TRKIN1_MELT_SUB, TRKIN2_MELT_SUB)

#Plot
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

