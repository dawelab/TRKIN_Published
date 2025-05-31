library(tidyverse)

setwd("")

#Read in the Ab10 TE annotation
Ab10_TE <- read.delim("B73_Ab10_HiFi_v2.TE.gff", header=FALSE)
colnames(Ab10_TE) <- c("seqnames", "method", "repeattype", "start", "end", "something1", "strand", "dot", "description")

#Isolate TEs within Ab10 trkin 1 
Ab10TE_trkin1 <- subset(Ab10_TE, seqnames == "chr10" & start >= 148968414 & end <= 149079834 )
#Separate out the superfamily identity
Ab10TE_trkin1 <- separate(Ab10TE_trkin1, description, into=c("ID", "Motif"), sep=" ")
Ab10TE_trkin1$Motif <- gsub("Motif:", "", Ab10TE_trkin1$Motif)

#write the Ab10 trkin1 TE annotation
write.csv(Ab10TE_trkin1, "Ab10_trkin1_TE.csv", row.names = FALSE, quote = FALSE)

#Isolate TEs within Ab10 trkin 2
Ab10TE_trkin2 <- subset(Ab10_TE, seqnames == "chr10" & start >= 150358418 & end <= 150457939)
#Separate out the superfamily identity
Ab10TE_trkin2 <- separate(Ab10TE_trkin2, description, into=c("ID", "Motif"), sep=" ")
Ab10TE_trkin2$Motif <- gsub("Motif:", "", Ab10TE_trkin2$Motif)

#write the Ab10 trkin2 TE annotation
write.csv(Ab10TE_trkin2, "Ab10_trkin2_TE.csv", row.names = FALSE, quote = FALSE)

#Read in the K10L2 TE annotation
K10L2_TE <- read.delim("CI66_K10L2_v1.TE.gff", header=FALSE)
colnames(K10L2_TE) <- c("seqnames", "method", "repeattype", "start", "end", "something1", "strand", "dot", "description")

#Isolate TEs within K10l2 trkin
K10L2TE_trkin <- subset(K10L2_TE, seqnames == "K10L2" & start >= 16280429 & end <= 16357330)
#Separate out the superfamily identity
K10L2TE_trkin <- separate(K10L2TE_trkin, description, into=c("ID", "Motif"), sep=" ")
K10L2TE_trkin$Motif <- gsub("Motif:", "", K10L2TE_trkin$Motif)


#write the K10L2 trkin TE annotation
write.csv(K10L2TE_trkin, "K10L2_trkin_TE.csv", row.names = FALSE, col.names = TRUE, quote = FALSE)

