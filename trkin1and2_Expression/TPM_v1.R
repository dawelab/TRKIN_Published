library(ggplot2)
library(reshape2)

Counts <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/Differential_Expression/FeatureCounts.csv")
Anno <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/Differential_Expression/FeatureCounts_Annot.csv")

Counts$chr <- Anno$Chr
Counts$start <- Anno$Start
Counts$end <- Anno$End
Counts$strand <- Anno$Strand
Counts$Length <- Anno$Length 

Counts$LengthKB <- Counts$Length/1000

tpm3 <- function(counts,len) {x <- counts/len
return(t(t(x)*1e6/colSums(x)))
}

TPM <- tpm3(Counts[,2:21],Counts$LengthKB)
TPM_DF <- as.data.frame(TPM)
TPM_DF$gene_ID <- Counts$X
colnames(TPM_DF) <- c("embryo_1", "embryo_2", "endosperm_1", "endosperm_2", "root_1", "root_2", "shoot_1", "shoot_2", "anther_1", "anther_2", "base_1", "base_2", "middle_1", "middle_2", "tip_1", "tip_2", "ear_1", "ear_2", "tassel_1", "tassel_2", "gene_ID")

TRKIN_TPM <- subset(TPM_DF, gene_ID=="gene:Zm00043a049192")
TRKIN_TPM$Exon <- c(1:19)
TRKIN_MELT <- melt(TRKIN_TPM, id=c("gene_ID", "Exon"))
TRKIN_MELT$Tissue <- c(rep("embryo", 19*2), rep("endosperm", 19*2), rep("root", 19*2), rep("shoot", 19*2), rep("anther", 19*2), rep("leaf base", 19*2), rep("leaf middle", 19*2), rep("leaf tip", 19*2), rep("ear", 19*2), rep("tassel", 19*2))

TRKIN_MELT_SUB <- subset(TRKIN_MELT, Exon == 7 | Exon == 8)

pdf(file="Ab10_TRKIN_1_7_8.pdf", height=8, width=2)
ggplot(TRKIN_MELT_SUB, aes(x=Exon, y=value, color = Tissue)) +
  geom_point(alpha=0.5) +
  facet_wrap(~Tissue, ncol=1) +
  ylab("TPM") +
  ggtitle("Ab10 trkin 1") +
  ylim(0,7) +
  scale_x_continuous(breaks=c(7:8), minor_breaks = c(0)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
dev.off()


PTRKIN_TPM <- subset(TPM_DF, gene_ID=="gene:Zm00043a049309")
PTRKIN_TPM$Exon <- c(19:1)
PTRKIN_MELT <- melt(PTRKIN_TPM, id=c("gene_ID", "Exon"))
PTRKIN_MELT$Tissue <- c(rep("embryo", 19*2), rep("endosperm", 19*2), rep("root", 19*2), rep("shoot", 19*2), rep("anther", 19*2), rep("leaf base", 19*2), rep("leaf middle", 19*2), rep("leaf tip", 19*2), rep("ear", 19*2), rep("tassel", 19*2))

PTRKIN_MELT_SUB <- subset(PTRKIN_MELT, Exon == 7 | Exon == 8)

pdf(file="Ab10_TRKIN_2_7_8.pdf", height=8, width=2)
ggplot(PTRKIN_MELT_SUB, aes(x=Exon, y=value, color=Tissue)) +
  geom_point(alpha=0.5) +
  facet_wrap(~Tissue, ncol=1) +
  ylab("TPM") +
  ggtitle("Ab10 trkin 2") +
  ylim(0,7) +
  scale_x_continuous(breaks=c(7:8), minor_breaks = c(0)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
dev.off()

t.test(TRKIN_MELT_SUB$value, PTRKIN_MELT_SUB$value)
