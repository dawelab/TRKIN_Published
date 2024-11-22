library(tidyverse)

setwd("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/Paper/TRKIN_Paper")

Ab10_TE <- read.delim("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/OrthoFinder/Ab10_HiFi_v2_corrected.TE.gff", header=FALSE)
colnames(Ab10_TE) <- c("seqnames", "method", "repeattype", "start", "end", "something1", "strand", "dot", "description")

Ab10TE_trkin1 <- subset(Ab10_TE, seqnames == "chr10" & start >= 148968414 & end <= 149079834 )
Ab10TE_trkin1 <- separate(Ab10TE_trkin1, description, into=c("ID", "Motif"), sep=" ")
Ab10TE_trkin1$Motif <- gsub("Motif:", "", Ab10TE_trkin1$Motif)
Uniq <- unique(Ab10TE_trkin1$Motif)

write.csv(Ab10TE_trkin1, "Ab10_trkin1_TE.csv", row.names = FALSE, quote = FALSE)


Ab10TE_trkin2 <- subset(Ab10_TE, seqnames == "chr10" & start >= 150358418 & end <= 150457939)
Ab10TE_trkin2 <- separate(Ab10TE_trkin2, description, into=c("ID", "Motif"), sep=" ")
Ab10TE_trkin2$Motif <- gsub("Motif:", "", Ab10TE_trkin2$Motif)
Uniq <- unique(Ab10TE_trkin2$Motif)

write.csv(Ab10TE_trkin2, "Ab10_trkin2_TE.csv", row.names = FALSE, quote = FALSE)


K10L2_TE <- read.delim("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/OrthoFinder/CI66_K10L2_v1.TE.gff", header=FALSE)
colnames(K10L2_TE) <- c("seqnames", "method", "repeattype", "start", "end", "something1", "strand", "dot", "description")


K10L2TE_trkin <- subset(K10L2_TE, seqnames == "K10L2" & start >= 16280429 & end <= 16357330)
K10L2TE_trkin <- separate(K10L2TE_trkin, description, into=c("ID", "Motif"), sep=" ")
K10L2TE_trkin$Motif <- gsub("Motif:", "", K10L2TE_trkin$Motif)
Uniq <- unique(K10L2TE_trkin$Motif)

write.csv(K10L2TE_trkin, "K10L2_trkin_TE.csv", row.names = FALSE, col.names = TRUE, quote = FALSE)

